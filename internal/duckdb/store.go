// Package duckdb provides caching for transcript and variant annotation data.
// Transcripts are cached as gob files (fast, pure Go).
// Variant results are cached in DuckDB (queryable, append-only).
package duckdb

import (
	"database/sql"
	"fmt"
	"os"
	"path/filepath"

	_ "github.com/marcboeker/go-duckdb"
)

// Store manages a DuckDB connection for caching variant results.
type Store struct {
	db   *sql.DB
	path string
}

// Open opens or creates a DuckDB database at the given path.
// Use an empty string for an in-memory database.
func Open(path string) (*Store, error) {
	if path != "" {
		dir := filepath.Dir(path)
		if err := os.MkdirAll(dir, 0755); err != nil {
			return nil, fmt.Errorf("create cache directory: %w", err)
		}
	}

	db, err := sql.Open("duckdb", path)
	if err != nil {
		return nil, fmt.Errorf("open duckdb: %w", err)
	}

	s := &Store{db: db, path: path}
	if err := s.ensureSchema(); err != nil {
		db.Close()
		return nil, fmt.Errorf("ensure schema: %w", err)
	}

	return s, nil
}

// Close closes the database connection.
func (s *Store) Close() error {
	return s.db.Close()
}

// DB returns the underlying *sql.DB for direct access.
func (s *Store) DB() *sql.DB {
	return s.db
}

// ensureSchema creates tables if they don't exist.
func (s *Store) ensureSchema() error {
	_, err := s.db.Exec(`CREATE TABLE IF NOT EXISTS variant_results (
		chrom VARCHAR,
		pos BIGINT,
		ref VARCHAR,
		alt VARCHAR,
		transcript_id VARCHAR,
		gene_name VARCHAR,
		gene_id VARCHAR,
		consequence VARCHAR,
		impact VARCHAR,
		cds_position BIGINT,
		protein_position BIGINT,
		amino_acid_change VARCHAR,
		codon_change VARCHAR,
		is_canonical BOOLEAN,
		allele VARCHAR,
		biotype VARCHAR,
		exon_number VARCHAR,
		intron_number VARCHAR,
		cdna_position BIGINT,
		hgvsp VARCHAR,
		hgvsc VARCHAR,
		gene_type VARCHAR,
		PRIMARY KEY (chrom, pos, ref, alt, transcript_id)
	)`)
	return err
}
