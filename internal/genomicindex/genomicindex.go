package genomicindex

import (
	"bufio"
	"compress/gzip"
	"database/sql"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/datasource/clinvar"
	"github.com/inodb/vibe-vep/internal/datasource/gnomad"
	"github.com/inodb/vibe-vep/internal/datasource/signal"
	_ "modernc.org/sqlite"
)

// Store provides point lookups against the unified genomic annotation SQLite database.
type Store struct {
	db       *sql.DB
	lookupPS *sql.Stmt
}

// Open opens an existing genomic annotation database and prepares the lookup statement.
func Open(dbPath string) (*Store, error) {
	db, err := sql.Open("sqlite", dbPath+"?_pragma=mmap_size%3D2147483648&_pragma=journal_mode%3DWAL")
	if err != nil {
		return nil, fmt.Errorf("open sqlite: %w", err)
	}

	ps, err := db.Prepare(`SELECT am_score, am_class, cv_clnsig, cv_revstat, cv_clndn,
		sig_mut_status, sig_count, sig_freq,
		gnomad_af, gnomad_ac, gnomad_an, gnomad_nhomalt, gnomad_version
		FROM genomic_annotations WHERE chrom=? AND pos=? AND ref=? AND alt=?`)
	if err != nil {
		db.Close()
		return nil, fmt.Errorf("prepare lookup: %w", err)
	}

	return &Store{db: db, lookupPS: ps}, nil
}

// Lookup performs a point lookup for a single variant.
func (s *Store) Lookup(chrom string, pos int64, ref, alt string) (Result, bool) {
	var r Result
	err := s.lookupPS.QueryRow(chrom, pos, ref, alt).Scan(
		&r.AMScore, &r.AMClass, &r.CVClnSig, &r.CVClnRevStat, &r.CVClnDN,
		&r.SigMutStatus, &r.SigCount, &r.SigFreq,
		&r.GnomadAF, &r.GnomadAC, &r.GnomadAN, &r.GnomadNhomalt, &r.GnomadVersion,
	)
	if err != nil {
		return Result{}, false
	}
	return r, true
}

// Close closes the prepared statement and database.
func (s *Store) Close() error {
	if s.lookupPS != nil {
		s.lookupPS.Close()
	}
	return s.db.Close()
}

// Ready returns true if dbPath exists, is newer than all source files, and
// passes a quick integrity check (table exists and is readable).
func Ready(dbPath string, sources BuildSources) bool {
	dbInfo, err := os.Stat(dbPath)
	if err != nil {
		return false
	}

	// Reject empty files (e.g. leftover from a failed build).
	if dbInfo.Size() == 0 {
		return false
	}

	dbMod := dbInfo.ModTime()

	for _, src := range []string{sources.AlphaMissenseTSV, sources.ClinVarVCF, sources.SignalTSV, sources.GnomadVCF} {
		if src == "" {
			continue
		}
		srcInfo, err := os.Stat(src)
		if err != nil {
			continue // source file missing — can't be stale
		}
		if srcInfo.ModTime().After(dbMod) {
			return false
		}
	}

	// Quick integrity check: open DB and verify the table is readable.
	if err := quickCheck(dbPath); err != nil {
		return false
	}
	return true
}

// quickCheck opens the database and verifies the genomic_annotations table
// exists and can return a row. This catches corruption, truncation, and
// schema mismatches without scanning the full table.
func quickCheck(dbPath string) error {
	db, err := sql.Open("sqlite", dbPath+"?mode=ro")
	if err != nil {
		return err
	}
	defer db.Close()

	// Verify table exists and is readable.
	var n int
	return db.QueryRow("SELECT 1 FROM genomic_annotations LIMIT 1").Scan(&n)
}

// Build creates the SQLite database from source files. Each source is loaded
// in a single transaction. The database is created at dbPath (overwriting any
// existing file).
func Build(dbPath string, sources BuildSources, logf func(string, ...any)) error {
	// Remove old DB so we start fresh.
	os.Remove(dbPath)
	os.Remove(dbPath + "-wal")
	os.Remove(dbPath + "-shm")

	db, err := sql.Open("sqlite", dbPath+"?_pragma=mmap_size%3D2147483648")
	if err != nil {
		return fmt.Errorf("create sqlite: %w", err)
	}
	defer db.Close()

	// Tuning for bulk load. journal_mode=OFF avoids a multi-GB WAL file
	// during the 71.7M-row AlphaMissense insert. Safe because we're building
	// from scratch — a crash just means we rebuild. WAL is set on Open() for
	// concurrent read access at runtime.
	for _, pragma := range []string{
		"PRAGMA journal_mode = OFF",
		"PRAGMA synchronous = OFF",
		"PRAGMA temp_store = MEMORY",
		"PRAGMA cache_size = -64000", // 64MB
		"PRAGMA page_size = 8192",
	} {
		if _, err := db.Exec(pragma); err != nil {
			return fmt.Errorf("set pragma %q: %w", pragma, err)
		}
	}

	if _, err := db.Exec(`CREATE TABLE genomic_annotations (
		chrom TEXT NOT NULL,
		pos INTEGER NOT NULL,
		ref TEXT NOT NULL,
		alt TEXT NOT NULL,
		am_score REAL NOT NULL DEFAULT 0,
		am_class TEXT NOT NULL DEFAULT '',
		cv_clnsig TEXT NOT NULL DEFAULT '',
		cv_revstat TEXT NOT NULL DEFAULT '',
		cv_clndn TEXT NOT NULL DEFAULT '',
		sig_mut_status TEXT NOT NULL DEFAULT '',
		sig_count TEXT NOT NULL DEFAULT '',
		sig_freq TEXT NOT NULL DEFAULT '',
		gnomad_af TEXT NOT NULL DEFAULT '',
		gnomad_ac TEXT NOT NULL DEFAULT '',
		gnomad_an TEXT NOT NULL DEFAULT '',
		gnomad_nhomalt TEXT NOT NULL DEFAULT '',
		gnomad_version TEXT NOT NULL DEFAULT '',
		PRIMARY KEY (chrom, pos, ref, alt)
	) WITHOUT ROWID`); err != nil {
		return fmt.Errorf("create table: %w", err)
	}

	// 1. AlphaMissense
	if sources.AlphaMissenseTSV != "" {
		if _, err := os.Stat(sources.AlphaMissenseTSV); err == nil {
			logf("loading AlphaMissense from %s", sources.AlphaMissenseTSV)
			n, err := loadAlphaMissense(db, sources.AlphaMissenseTSV)
			if err != nil {
				return fmt.Errorf("load AlphaMissense: %w", err)
			}
			logf("loaded %d AlphaMissense variants", n)
		}
	}

	// 2. ClinVar
	if sources.ClinVarVCF != "" {
		if _, err := os.Stat(sources.ClinVarVCF); err == nil {
			logf("loading ClinVar from %s", sources.ClinVarVCF)
			n, err := loadClinVar(db, sources.ClinVarVCF)
			if err != nil {
				return fmt.Errorf("load ClinVar: %w", err)
			}
			logf("loaded %d ClinVar variants", n)
		}
	}

	// 3. SIGNAL
	if sources.SignalTSV != "" {
		if _, err := os.Stat(sources.SignalTSV); err == nil {
			logf("loading SIGNAL from %s", sources.SignalTSV)
			n, err := loadSignal(db, sources.SignalTSV)
			if err != nil {
				return fmt.Errorf("load SIGNAL: %w", err)
			}
			logf("loaded %d SIGNAL variants", n)
		}
	}

	// 4. gnomAD
	if sources.GnomadVCF != "" {
		if _, err := os.Stat(sources.GnomadVCF); err == nil {
			logf("loading gnomAD from %s", sources.GnomadVCF)
			n, err := loadGnomad(db, sources.GnomadVCF, sources.GnomadVersion)
			if err != nil {
				return fmt.Errorf("load gnomAD: %w", err)
			}
			logf("loaded %d gnomAD variants", n)
		}
	}

	return nil
}

// loadAlphaMissense parses a gzipped AlphaMissense TSV and inserts into the DB.
// The TSV has 3 comment lines, then a header, then data lines:
//
//	#CHROM  POS  REF  ALT  genome  uniprot_id  transcript_id  protein_variant  am_pathogenicity  am_class
//
// Multiple rows per (chrom,pos,ref,alt) exist (one per transcript) with identical scores.
// INSERT OR IGNORE deduplicates.
func loadAlphaMissense(db *sql.DB, tsvPath string) (int64, error) {
	f, err := os.Open(tsvPath)
	if err != nil {
		return 0, err
	}
	defer f.Close()

	gz, err := gzip.NewReader(f)
	if err != nil {
		return 0, err
	}
	defer gz.Close()

	scanner := bufio.NewScanner(gz)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	// Skip comment lines and header (4 lines total: 3 comments + 1 header).
	for i := 0; i < 4 && scanner.Scan(); i++ {
	}

	tx, err := db.Begin()
	if err != nil {
		return 0, err
	}
	defer tx.Rollback()

	stmt, err := tx.Prepare(`INSERT OR IGNORE INTO genomic_annotations (chrom, pos, ref, alt, am_score, am_class) VALUES (?, ?, ?, ?, ?, ?)`)
	if err != nil {
		return 0, err
	}
	defer stmt.Close()

	var count int64
	for scanner.Scan() {
		line := scanner.Text()
		// Fields: CHROM POS REF ALT genome uniprot_id transcript_id protein_variant am_pathogenicity am_class
		// Indices: 0     1   2   3   4      5          6             7               8                9
		chrom, rest, ok := cutTab(line)
		if !ok {
			continue
		}
		posStr, rest, ok := cutTab(rest)
		if !ok {
			continue
		}
		ref, rest, ok := cutTab(rest)
		if !ok {
			continue
		}
		alt, rest, ok := cutTab(rest)
		if !ok {
			continue
		}
		// Skip fields 4-7 (genome, uniprot_id, transcript_id, protein_variant).
		for i := 0; i < 4; i++ {
			_, rest, ok = cutTab(rest)
			if !ok {
				break
			}
		}
		if !ok {
			continue
		}
		scoreStr, classStr, ok := cutTab(rest)
		if !ok {
			continue
		}

		// Only SNVs (single-base ref and alt).
		if len(ref) != 1 || len(alt) != 1 {
			continue
		}

		pos, err := strconv.ParseInt(posStr, 10, 64)
		if err != nil {
			continue
		}
		score, err := strconv.ParseFloat(scoreStr, 32)
		if err != nil {
			continue
		}

		chrom = strings.TrimPrefix(chrom, "chr")
		if _, err := stmt.Exec(chrom, pos, ref, alt, float32(score), classStr); err != nil {
			return 0, fmt.Errorf("insert AM row: %w", err)
		}
		count++
	}
	if err := scanner.Err(); err != nil {
		return 0, err
	}

	if err := tx.Commit(); err != nil {
		return 0, err
	}
	return count, nil
}

// loadClinVar parses a gzipped ClinVar VCF and upserts into the DB.
func loadClinVar(db *sql.DB, vcfPath string) (int64, error) {
	f, err := os.Open(vcfPath)
	if err != nil {
		return 0, err
	}
	defer f.Close()

	var scanner *bufio.Scanner
	if strings.HasSuffix(vcfPath, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return 0, err
		}
		defer gz.Close()
		scanner = bufio.NewScanner(gz)
	} else {
		scanner = bufio.NewScanner(f)
	}
	scanner.Buffer(make([]byte, 4*1024*1024), 4*1024*1024)

	tx, err := db.Begin()
	if err != nil {
		return 0, err
	}
	defer tx.Rollback()

	// Upsert: if variant already exists (from AlphaMissense), update ClinVar fields.
	stmt, err := tx.Prepare(`INSERT INTO genomic_annotations (chrom, pos, ref, alt, cv_clnsig, cv_revstat, cv_clndn)
		VALUES (?, ?, ?, ?, ?, ?, ?)
		ON CONFLICT(chrom, pos, ref, alt) DO UPDATE SET
			cv_clnsig=excluded.cv_clnsig, cv_revstat=excluded.cv_revstat, cv_clndn=excluded.cv_clndn`)
	if err != nil {
		return 0, err
	}
	defer stmt.Close()

	var count int64
	for scanner.Scan() {
		line := scanner.Text()
		if len(line) == 0 || line[0] == '#' {
			continue
		}

		entry, chrom, ok := clinvar.ParseVCFLine(line)
		if !ok {
			continue
		}

		// Normalize VCF-style alleles to canonical (MAF-style) form.
		nPos, nRef, nAlt := NormalizeAlleles(entry.Pos, entry.Ref, entry.Alt)
		if _, err := stmt.Exec(chrom, nPos, nRef, nAlt, entry.ClnSig, entry.RevStat, entry.ClnDN); err != nil {
			return 0, fmt.Errorf("upsert ClinVar row: %w", err)
		}
		count++
	}
	if err := scanner.Err(); err != nil {
		return 0, err
	}

	if err := tx.Commit(); err != nil {
		return 0, err
	}
	return count, nil
}

// loadSignal parses a SIGNAL TSV and upserts into the DB.
func loadSignal(db *sql.DB, tsvPath string) (int64, error) {
	f, err := os.Open(tsvPath)
	if err != nil {
		return 0, err
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	// Read header.
	if !scanner.Scan() {
		return 0, fmt.Errorf("empty signal file")
	}
	header := strings.Split(scanner.Text(), "\t")
	col := signal.IndexColumns(header,
		"Chromosome", "Start_Position",
		"Reference_Allele", "Alternate_Allele",
		"n_impact", "f_impact",
	)
	if col["Chromosome"] < 0 || col["Start_Position"] < 0 {
		return 0, fmt.Errorf("missing required SIGNAL columns")
	}

	tx, err := db.Begin()
	if err != nil {
		return 0, err
	}
	defer tx.Rollback()

	stmt, err := tx.Prepare(`INSERT INTO genomic_annotations (chrom, pos, ref, alt, sig_mut_status, sig_count, sig_freq)
		VALUES (?, ?, ?, ?, ?, ?, ?)
		ON CONFLICT(chrom, pos, ref, alt) DO UPDATE SET
			sig_mut_status=excluded.sig_mut_status, sig_count=excluded.sig_count, sig_freq=excluded.sig_freq`)
	if err != nil {
		return 0, err
	}
	defer stmt.Close()

	var count int64
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		chromIdx := col["Chromosome"]
		posIdx := col["Start_Position"]
		if chromIdx >= len(fields) || posIdx >= len(fields) {
			continue
		}

		chrom := strings.TrimPrefix(fields[chromIdx], "chr")
		pos, err := strconv.ParseInt(fields[posIdx], 10, 64)
		if err != nil {
			continue
		}

		ref := getField(fields, col["Reference_Allele"])
		alt := getField(fields, col["Alternate_Allele"])
		if ref == "-" {
			ref = ""
		}
		if alt == "-" {
			alt = ""
		}

		var countStr, freqStr string
		if idx := col["n_impact"]; idx >= 0 && idx < len(fields) {
			n, _ := strconv.Atoi(fields[idx])
			if n > 0 {
				countStr = strconv.Itoa(n)
			}
		}
		if idx := col["f_impact"]; idx >= 0 && idx < len(fields) {
			freq, _ := strconv.ParseFloat(fields[idx], 64)
			if freq > 0 {
				freqStr = signal.FormatFreq(freq)
			}
		}

		if _, err := stmt.Exec(chrom, pos, ref, alt, "germline", countStr, freqStr); err != nil {
			return 0, fmt.Errorf("upsert SIGNAL row: %w", err)
		}
		count++
	}
	if err := scanner.Err(); err != nil {
		return 0, err
	}

	if err := tx.Commit(); err != nil {
		return 0, err
	}
	return count, nil
}

// loadGnomad parses a gzipped gnomAD VCF and upserts into the DB.
func loadGnomad(db *sql.DB, vcfPath, version string) (int64, error) {
	f, err := os.Open(vcfPath)
	if err != nil {
		return 0, err
	}
	defer f.Close()

	var scanner *bufio.Scanner
	if strings.HasSuffix(vcfPath, ".gz") || strings.HasSuffix(vcfPath, ".bgz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return 0, err
		}
		defer gz.Close()
		scanner = bufio.NewScanner(gz)
	} else {
		scanner = bufio.NewScanner(f)
	}
	scanner.Buffer(make([]byte, 4*1024*1024), 4*1024*1024)

	tx, err := db.Begin()
	if err != nil {
		return 0, err
	}
	defer tx.Rollback()

	stmt, err := tx.Prepare(`INSERT INTO genomic_annotations (chrom, pos, ref, alt, gnomad_af, gnomad_ac, gnomad_an, gnomad_nhomalt, gnomad_version)
		VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
		ON CONFLICT(chrom, pos, ref, alt) DO UPDATE SET
			gnomad_af=excluded.gnomad_af, gnomad_ac=excluded.gnomad_ac,
			gnomad_an=excluded.gnomad_an, gnomad_nhomalt=excluded.gnomad_nhomalt,
			gnomad_version=excluded.gnomad_version`)
	if err != nil {
		return 0, err
	}
	defer stmt.Close()

	var count int64
	for scanner.Scan() {
		line := scanner.Text()
		if len(line) == 0 || line[0] == '#' {
			continue
		}

		entry, chrom, ok := gnomad.ParseVCFLine(line)
		if !ok {
			continue
		}

		nPos, nRef, nAlt := NormalizeAlleles(entry.Pos, entry.Ref, entry.Alt)
		afStr := gnomad.FormatAF(entry.AF)
		acStr := strconv.Itoa(entry.AC)
		anStr := strconv.Itoa(entry.AN)
		nhomaltStr := strconv.Itoa(entry.Nhomalt)

		if _, err := stmt.Exec(chrom, nPos, nRef, nAlt, afStr, acStr, anStr, nhomaltStr, version); err != nil {
			return 0, fmt.Errorf("upsert gnomAD row: %w", err)
		}
		count++
	}
	if err := scanner.Err(); err != nil {
		return 0, err
	}

	if err := tx.Commit(); err != nil {
		return 0, err
	}
	return count, nil
}

// cutTab splits s at the first tab, returning the part before, the part after, and ok.
func cutTab(s string) (string, string, bool) {
	i := strings.IndexByte(s, '\t')
	if i < 0 {
		return s, "", false
	}
	return s[:i], s[i+1:], true
}

func getField(fields []string, idx int) string {
	if idx < 0 || idx >= len(fields) {
		return ""
	}
	return fields[idx]
}

// NormalizeAlleles converts VCF-style variants to a canonical (MAF-style) form
// by stripping shared prefix/suffix bases and adjusting position.
//
// VCF indels include an anchor base:
//
//	deletion:  POS=99  REF=GAT  ALT=G    → (100, "AT", "")
//	insertion: POS=100 REF=A    ALT=ACG  → (100, "",   "CG")
//
// MAF indels have no anchor:
//
//	deletion:  Start=100 Ref=AT  Alt=-   → (100, "AT", "")
//	insertion: Start=100 Ref=-   Alt=CG  → (100, "",   "CG")
//
// After normalization both produce the same key. Safe to call on already-
// normalized (MAF-derived) variants — it is a no-op.
func NormalizeAlleles(pos int64, ref, alt string) (int64, string, string) {
	// Left-trim common prefix.
	trimmed := 0
	for len(ref) > 0 && len(alt) > 0 && ref[0] == alt[0] {
		ref = ref[1:]
		alt = alt[1:]
		trimmed++
	}
	pos += int64(trimmed)

	// Right-trim common suffix.
	for len(ref) > 0 && len(alt) > 0 && ref[len(ref)-1] == alt[len(alt)-1] {
		ref = ref[:len(ref)-1]
		alt = alt[:len(alt)-1]
	}

	// For insertions where we trimmed the VCF anchor base, adjust pos back.
	// MAF convention: insertion pos = base before the inserted sequence.
	if trimmed > 0 && len(ref) == 0 && len(alt) > 0 {
		pos--
	}

	return pos, ref, alt
}
