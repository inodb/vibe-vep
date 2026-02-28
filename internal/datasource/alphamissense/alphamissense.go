// Package alphamissense provides AlphaMissense pathogenicity score lookups
// backed by DuckDB. AlphaMissense data is loaded from the official TSV files
// (Cheng et al., Science 2023, CC BY 4.0).
package alphamissense

import (
	"database/sql"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	_ "github.com/marcboeker/go-duckdb"
)

// Store provides AlphaMissense score lookups backed by DuckDB.
type Store struct {
	db       *sql.DB
	lookupPS *sql.Stmt // prepared statement for Lookup, lazily initialized
}

// Open opens or creates a DuckDB database for AlphaMissense data at the given path.
func Open(dbPath string) (*Store, error) {
	if dbPath != "" {
		dir := filepath.Dir(dbPath)
		if err := os.MkdirAll(dir, 0755); err != nil {
			return nil, fmt.Errorf("create directory: %w", err)
		}
	}

	db, err := sql.Open("duckdb", dbPath)
	if err != nil {
		return nil, fmt.Errorf("open duckdb: %w", err)
	}

	s := &Store{db: db}
	if err := s.ensureSchema(); err != nil {
		db.Close()
		return nil, fmt.Errorf("ensure schema: %w", err)
	}

	return s, nil
}

func (s *Store) ensureSchema() error {
	if _, err := s.db.Exec(`CREATE TABLE IF NOT EXISTS alphamissense (
		chrom VARCHAR,
		pos BIGINT,
		ref VARCHAR,
		alt VARCHAR,
		am_pathogenicity FLOAT,
		am_class VARCHAR
	)`); err != nil {
		return err
	}
	// Index for fast point lookups
	s.db.Exec(`CREATE INDEX IF NOT EXISTS idx_am_lookup ON alphamissense (chrom, pos, ref, alt)`)
	return nil
}

// Loaded returns true if the AlphaMissense table has data.
func (s *Store) Loaded() bool {
	var count int64
	err := s.db.QueryRow("SELECT COUNT(*) FROM alphamissense").Scan(&count)
	return err == nil && count > 0
}

// Count returns the number of rows in the AlphaMissense table.
func (s *Store) Count() (int64, error) {
	var count int64
	err := s.db.QueryRow("SELECT COUNT(*) FROM alphamissense").Scan(&count)
	if err != nil {
		return 0, fmt.Errorf("count alphamissense rows: %w", err)
	}
	return count, nil
}

// Load bulk-loads AlphaMissense data from a gzipped TSV file using DuckDB's read_csv.
// The file has 3 comment lines, then a header:
//
//	#CHROM  POS  REF  ALT  genome  uniprot_id  transcript_id  protein_variant  am_pathogenicity  am_class
func (s *Store) Load(tsvPath string) error {
	// Clear any existing data first (idempotent reload)
	s.db.Exec(`DELETE FROM alphamissense`)

	// Stream directly from the gzipped TSV. The file has multiple rows per
	// (chrom,pos,ref,alt) — one per transcript — but scores are identical.
	// Duplicates are harmless; Lookup uses LIMIT 1.
	query := fmt.Sprintf(`INSERT INTO alphamissense
		SELECT column0, column1, column2, column3,
			CAST(column8 AS FLOAT), column9
		FROM read_csv('%s', delim='\t', header=false, skip=4,
			columns={
				'column0': 'VARCHAR',
				'column1': 'BIGINT',
				'column2': 'VARCHAR',
				'column3': 'VARCHAR',
				'column4': 'VARCHAR',
				'column5': 'VARCHAR',
				'column6': 'VARCHAR',
				'column7': 'VARCHAR',
				'column8': 'VARCHAR',
				'column9': 'VARCHAR'
			})`, tsvPath)

	if _, err := s.db.Exec(query); err != nil {
		return fmt.Errorf("loading AlphaMissense data: %w", err)
	}
	return nil
}

// Result holds a single AlphaMissense lookup result.
type Result struct {
	Score float64
	Class string
}

// Lookup queries the AlphaMissense score for a specific variant.
// It lazily prepares the lookup statement for efficient repeated calls.
func (s *Store) Lookup(chrom string, pos int64, ref, alt string) (Result, bool) {
	if s.lookupPS == nil {
		ps, err := s.db.Prepare(
			"SELECT am_pathogenicity, am_class FROM alphamissense WHERE chrom=? AND pos=? AND ref=? AND alt=? LIMIT 1",
		)
		if err != nil {
			return Result{}, false
		}
		s.lookupPS = ps
	}
	var r Result
	err := s.lookupPS.QueryRow(chrom, pos, ref, alt).Scan(&r.Score, &r.Class)
	if err != nil {
		return Result{}, false
	}
	return r, true
}

// LookupKey identifies a variant for batch lookup.
type LookupKey struct {
	Chrom string
	Pos   int64
	Ref   string
	Alt   string
}

// BatchLookup queries AlphaMissense scores for a batch of variants using a
// temporary table join for high throughput. Returns a map from LookupKey to Result.
func (s *Store) BatchLookup(keys []LookupKey) (map[LookupKey]Result, error) {
	if len(keys) == 0 {
		return nil, nil
	}

	// Create temp table and insert keys.
	if _, err := s.db.Exec(`CREATE TEMPORARY TABLE IF NOT EXISTS batch_keys (
		chrom VARCHAR, pos BIGINT, ref VARCHAR, alt VARCHAR
	)`); err != nil {
		return nil, fmt.Errorf("create temp table: %w", err)
	}
	defer s.db.Exec(`DROP TABLE IF EXISTS batch_keys`)
	s.db.Exec(`DELETE FROM batch_keys`)

	// Insert in chunks to avoid huge SQL statements.
	const chunkSize = 1000
	for i := 0; i < len(keys); i += chunkSize {
		end := i + chunkSize
		if end > len(keys) {
			end = len(keys)
		}
		chunk := keys[i:end]

		var sb strings.Builder
		sb.WriteString("INSERT INTO batch_keys VALUES ")
		for j, k := range chunk {
			if j > 0 {
				sb.WriteByte(',')
			}
			fmt.Fprintf(&sb, "('%s',%d,'%s','%s')", k.Chrom, k.Pos, k.Ref, k.Alt)
		}
		if _, err := s.db.Exec(sb.String()); err != nil {
			return nil, fmt.Errorf("insert batch keys: %w", err)
		}
	}

	// Join and collect results.
	rows, err := s.db.Query(`
		SELECT b.chrom, b.pos, b.ref, b.alt, a.am_pathogenicity, a.am_class
		FROM batch_keys b
		JOIN alphamissense a ON a.chrom=b.chrom AND a.pos=b.pos AND a.ref=b.ref AND a.alt=b.alt
	`)
	if err != nil {
		return nil, fmt.Errorf("batch lookup query: %w", err)
	}
	defer rows.Close()

	results := make(map[LookupKey]Result, len(keys))
	for rows.Next() {
		var k LookupKey
		var r Result
		if err := rows.Scan(&k.Chrom, &k.Pos, &k.Ref, &k.Alt, &r.Score, &r.Class); err != nil {
			return nil, fmt.Errorf("scan batch result: %w", err)
		}
		if _, exists := results[k]; !exists {
			results[k] = r
		}
	}
	if err := rows.Err(); err != nil {
		return nil, fmt.Errorf("batch lookup rows: %w", err)
	}
	return results, nil
}

// Close closes the database connection.
func (s *Store) Close() error {
	if s.lookupPS != nil {
		s.lookupPS.Close()
	}
	return s.db.Close()
}
