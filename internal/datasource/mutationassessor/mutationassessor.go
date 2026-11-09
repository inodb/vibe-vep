// Package mutationassessor serves Mutation Assessor's functional impact
// predictions for amino acid substitutions.
//
// The predictions are distributed as four gzipped TSVs keyed by UniProt
// accession and protein change, which is how they are looked up here: a
// variant's transcript gives a UniProt id, and its protein change gives the
// substitution. genome-nexus imports the same files into MongoDB, gunzipping
// roughly eight gigabytes to do it; this reads them as a stream and writes
// only the index.
//
// The Mutation Assessor web service is long gone, so these files are the only
// remaining source and will not change again.
package mutationassessor

import (
	"bufio"
	"compress/gzip"
	"database/sql"
	"fmt"
	"os"
	"strconv"
	"strings"

	_ "modernc.org/sqlite"
)

// Annotation is one Mutation Assessor prediction.
type Annotation struct {
	UniprotID                  string
	SV                         int
	HGVSpShort                 string
	FunctionalImpactScore      float64
	FunctionalImpactPrediction string
	MSA                        string
	MAV                        int
}

// Store answers predictions from a built index.
type Store struct {
	db       *sql.DB
	lookupPS *sql.Stmt
}

// Open opens a built index and prepares the lookup.
func Open(dbPath string) (*Store, error) {
	db, err := sql.Open("sqlite", dbPath+"?_pragma=mmap_size%3D1073741824&_pragma=journal_mode%3DWAL")
	if err != nil {
		return nil, fmt.Errorf("open sqlite: %w", err)
	}
	ps, err := db.Prepare(`SELECT sv, f_score, f_impact, msa, mav
		FROM mutation_assessor WHERE uniprot_id=? AND hgvsp_short=?`)
	if err != nil {
		db.Close()
		return nil, fmt.Errorf("prepare lookup: %w", err)
	}
	return &Store{db: db, lookupPS: ps}, nil
}

// Close releases the database.
func (s *Store) Close() error {
	if s == nil || s.db == nil {
		return nil
	}
	return s.db.Close()
}

// Lookup returns the prediction for a protein change on a UniProt entry.
// hgvspShort is the one-letter form, e.g. "p.L25A".
func (s *Store) Lookup(uniprotID, hgvspShort string) (Annotation, bool) {
	if s == nil || s.lookupPS == nil || uniprotID == "" || hgvspShort == "" {
		return Annotation{}, false
	}
	a := Annotation{UniprotID: uniprotID, HGVSpShort: hgvspShort}
	err := s.lookupPS.QueryRow(uniprotID, hgvspShort).Scan(
		&a.SV, &a.FunctionalImpactScore, &a.FunctionalImpactPrediction, &a.MSA, &a.MAV)
	if err != nil {
		return Annotation{}, false
	}
	return a, true
}

// Ready reports whether a usable index exists at dbPath.
func Ready(dbPath string) bool {
	info, err := os.Stat(dbPath)
	if err != nil || info.Size() == 0 {
		return false
	}
	db, err := sql.Open("sqlite", dbPath)
	if err != nil {
		return false
	}
	defer db.Close()
	var one int
	return db.QueryRow(`SELECT 1 FROM mutation_assessor LIMIT 1`).Scan(&one) == nil
}

// Build streams the distributed TSVs into an index at dbPath, replacing any
// existing one. The inputs stay gzipped throughout.
func Build(dbPath string, tsvPaths []string) error {
	if err := os.Remove(dbPath); err != nil && !os.IsNotExist(err) {
		return fmt.Errorf("remove existing index: %w", err)
	}

	db, err := sql.Open("sqlite", dbPath)
	if err != nil {
		return fmt.Errorf("open sqlite: %w", err)
	}
	defer db.Close()

	// Bulk-load tuning, as for the genomic index: the journal is pointless
	// when a crash just means rebuilding from the source files.
	for _, pragma := range []string{
		"PRAGMA journal_mode = OFF",
		"PRAGMA synchronous = OFF",
		"PRAGMA temp_store = MEMORY",
		"PRAGMA cache_size = -64000",
		"PRAGMA page_size = 8192",
	} {
		if _, err := db.Exec(pragma); err != nil {
			return fmt.Errorf("set pragma %q: %w", pragma, err)
		}
	}

	if _, err := db.Exec(`CREATE TABLE mutation_assessor (
		uniprot_id TEXT NOT NULL,
		hgvsp_short TEXT NOT NULL,
		sv INTEGER NOT NULL DEFAULT 0,
		f_score REAL NOT NULL DEFAULT 0,
		f_impact TEXT NOT NULL DEFAULT '',
		msa TEXT NOT NULL DEFAULT '',
		mav INTEGER NOT NULL DEFAULT 0,
		PRIMARY KEY (uniprot_id, hgvsp_short)
	) WITHOUT ROWID`); err != nil {
		return fmt.Errorf("create table: %w", err)
	}

	for _, path := range tsvPaths {
		if err := loadTSV(db, path); err != nil {
			return fmt.Errorf("load %s: %w", path, err)
		}
	}
	return nil
}

// loadTSV streams one gzipped TSV into the index in a single transaction.
func loadTSV(db *sql.DB, path string) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()

	gz, err := gzip.NewReader(f)
	if err != nil {
		return fmt.Errorf("gzip: %w", err)
	}
	defer gz.Close()

	tx, err := db.Begin()
	if err != nil {
		return err
	}
	defer tx.Rollback()

	// The same substitution can appear in more than one of the four files;
	// the first prediction wins rather than failing the load.
	stmt, err := tx.Prepare(`INSERT OR IGNORE INTO mutation_assessor
		(uniprot_id, hgvsp_short, sv, f_score, f_impact, msa, mav)
		VALUES (?, ?, ?, ?, ?, ?, ?)`)
	if err != nil {
		return err
	}
	defer stmt.Close()

	scanner := bufio.NewScanner(gz)
	scanner.Buffer(make([]byte, 0, 64*1024), 1024*1024)

	if !scanner.Scan() {
		return fmt.Errorf("empty file")
	}
	cols := indexColumns(strings.Split(scanner.Text(), "\t"),
		"uniprotId", "sv", "hgvspShort", "f_score", "f_impact", "msa", "mav")
	for _, name := range []string{"uniprotId", "hgvspShort"} {
		if cols[name] < 0 {
			return fmt.Errorf("missing required column %q", name)
		}
	}

	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		get := func(name string) string {
			i := cols[name]
			if i < 0 || i >= len(fields) {
				return ""
			}
			return fields[i]
		}
		uniprotID, hgvspShort := get("uniprotId"), get("hgvspShort")
		if uniprotID == "" || hgvspShort == "" {
			continue
		}
		sv, _ := strconv.Atoi(get("sv"))
		score, _ := strconv.ParseFloat(get("f_score"), 64)
		mav, _ := strconv.Atoi(get("mav"))
		if _, err := stmt.Exec(uniprotID, hgvspShort, sv, score, get("f_impact"), get("msa"), mav); err != nil {
			return fmt.Errorf("insert: %w", err)
		}
	}
	if err := scanner.Err(); err != nil {
		return fmt.Errorf("scan: %w", err)
	}
	return tx.Commit()
}

// indexColumns maps column names to their position, or -1 when absent.
func indexColumns(header []string, names ...string) map[string]int {
	idx := make(map[string]int, len(names))
	for _, name := range names {
		idx[name] = -1
	}
	for i, col := range header {
		col = strings.TrimSpace(col)
		if _, ok := idx[col]; ok {
			idx[col] = i
		}
	}
	return idx
}
