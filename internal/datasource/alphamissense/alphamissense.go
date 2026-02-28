// Package alphamissense provides AlphaMissense pathogenicity score lookups
// backed by DuckDB. AlphaMissense data is loaded from the official TSV files
// (Cheng et al., Science 2023, CC BY 4.0).
package alphamissense

import (
	"database/sql"
	"fmt"
	"os"
	"path/filepath"

	_ "github.com/marcboeker/go-duckdb"
)

// Store provides AlphaMissense score lookups backed by DuckDB.
type Store struct {
	db *sql.DB
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
	_, err := s.db.Exec(`CREATE TABLE IF NOT EXISTS alphamissense (
		chrom VARCHAR,
		pos BIGINT,
		ref VARCHAR,
		alt VARCHAR,
		am_pathogenicity FLOAT,
		am_class VARCHAR,
		PRIMARY KEY (chrom, pos, ref, alt)
	)`)
	return err
}

// Loaded returns true if the AlphaMissense table has data.
func (s *Store) Loaded() bool {
	var count int64
	err := s.db.QueryRow("SELECT COUNT(*) FROM alphamissense").Scan(&count)
	return err == nil && count > 0
}

// Load bulk-loads AlphaMissense data from a gzipped TSV file using DuckDB's read_csv.
func (s *Store) Load(tsvPath string) error {
	// The AlphaMissense TSV has 4 header comment lines starting with #,
	// then the actual header line: #CHROM  POS  REF  ALT  am_pathogenicity  am_class
	// DuckDB's read_csv with skip=3 and header=true handles this.
	query := fmt.Sprintf(`INSERT INTO alphamissense
		SELECT "#CHROM", "POS", "REF", "ALT",
			CAST(am_pathogenicity AS FLOAT), am_class
		FROM read_csv('%s', delim='\t', header=true, skip=3,
			columns={
				'#CHROM': 'VARCHAR',
				'POS': 'BIGINT',
				'REF': 'VARCHAR',
				'ALT': 'VARCHAR',
				'am_pathogenicity': 'VARCHAR',
				'am_class': 'VARCHAR'
			})`, tsvPath)

	_, err := s.db.Exec(query)
	if err != nil {
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
func (s *Store) Lookup(chrom string, pos int64, ref, alt string) (Result, bool) {
	var r Result
	err := s.db.QueryRow(
		"SELECT am_pathogenicity, am_class FROM alphamissense WHERE chrom=? AND pos=? AND ref=? AND alt=?",
		chrom, pos, ref, alt,
	).Scan(&r.Score, &r.Class)
	if err != nil {
		return Result{}, false
	}
	return r, true
}

// Close closes the database connection.
func (s *Store) Close() error {
	return s.db.Close()
}
