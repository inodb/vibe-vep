// Package cache provides VEP cache loading functionality.
package cache

import (
	"database/sql"
	"encoding/json"
	"fmt"
	"strings"

	_ "github.com/marcboeker/go-duckdb"
)

// DuckDBLoader provides access to transcript data stored in a DuckDB database.
type DuckDBLoader struct {
	db   *sql.DB
	path string
}

// NewDuckDBLoader creates a new DuckDB-backed transcript loader.
// The path can be a local file path or an S3 URL (s3://bucket/path.duckdb).
func NewDuckDBLoader(path string) (*DuckDBLoader, error) {
	// Configure connection string for DuckDB
	connStr := path

	// For S3 URLs, we need to enable httpfs extension
	if strings.HasPrefix(path, "s3://") {
		// DuckDB will automatically use httpfs for S3 URLs
		connStr = path
	}

	db, err := sql.Open("duckdb", connStr)
	if err != nil {
		return nil, fmt.Errorf("open duckdb: %w", err)
	}

	// Enable httpfs extension for S3 support
	if strings.HasPrefix(path, "s3://") {
		if _, err := db.Exec("INSTALL httpfs; LOAD httpfs;"); err != nil {
			db.Close()
			return nil, fmt.Errorf("load httpfs extension: %w", err)
		}
	}

	return &DuckDBLoader{
		db:   db,
		path: path,
	}, nil
}

// Close closes the database connection.
func (l *DuckDBLoader) Close() error {
	return l.db.Close()
}

// Load loads all transcripts for a chromosome into the cache.
func (l *DuckDBLoader) Load(c *Cache, chrom string) error {
	transcripts, err := l.findTranscriptsByChrom(chrom)
	if err != nil {
		return err
	}

	for _, t := range transcripts {
		c.AddTranscript(t)
	}
	return nil
}

// LoadAll loads all transcripts into the cache.
func (l *DuckDBLoader) LoadAll(c *Cache) error {
	rows, err := l.db.Query(`
		SELECT id, gene_id, gene_name, chrom, start, end_, strand, biotype,
		       is_canonical, cds_start, cds_end, cds_sequence, protein_sequence
		FROM transcripts
		ORDER BY chrom, start
	`)
	if err != nil {
		return fmt.Errorf("query transcripts: %w", err)
	}
	defer rows.Close()

	for rows.Next() {
		t, err := l.scanTranscript(rows)
		if err != nil {
			return err
		}
		if err := l.loadExons(t); err != nil {
			return err
		}
		c.AddTranscript(t)
	}
	return rows.Err()
}

// FindTranscripts returns all transcripts overlapping a genomic position.
func (l *DuckDBLoader) FindTranscripts(chrom string, pos int64) ([]*Transcript, error) {
	rows, err := l.db.Query(`
		SELECT id, gene_id, gene_name, chrom, start, end_, strand, biotype,
		       is_canonical, cds_start, cds_end, cds_sequence, protein_sequence
		FROM transcripts
		WHERE chrom = ? AND start <= ? AND end_ >= ?
		ORDER BY start
	`, chrom, pos, pos)
	if err != nil {
		return nil, fmt.Errorf("query transcripts: %w", err)
	}
	defer rows.Close()

	var transcripts []*Transcript
	for rows.Next() {
		t, err := l.scanTranscript(rows)
		if err != nil {
			return nil, err
		}
		if err := l.loadExons(t); err != nil {
			return nil, err
		}
		transcripts = append(transcripts, t)
	}
	return transcripts, rows.Err()
}

// GetTranscript returns a specific transcript by ID.
func (l *DuckDBLoader) GetTranscript(id string) (*Transcript, error) {
	row := l.db.QueryRow(`
		SELECT id, gene_id, gene_name, chrom, start, end_, strand, biotype,
		       is_canonical, cds_start, cds_end, cds_sequence, protein_sequence
		FROM transcripts
		WHERE id = ?
	`, id)

	t := &Transcript{}
	var cdsSeq, protSeq sql.NullString
	err := row.Scan(
		&t.ID, &t.GeneID, &t.GeneName, &t.Chrom, &t.Start, &t.End,
		&t.Strand, &t.Biotype, &t.IsCanonical, &t.CDSStart, &t.CDSEnd,
		&cdsSeq, &protSeq,
	)
	if err == sql.ErrNoRows {
		return nil, nil
	}
	if err != nil {
		return nil, fmt.Errorf("scan transcript: %w", err)
	}

	t.CDSSequence = cdsSeq.String
	t.ProteinSequence = protSeq.String

	if err := l.loadExons(t); err != nil {
		return nil, err
	}
	return t, nil
}

// findTranscriptsByChrom returns all transcripts for a chromosome.
func (l *DuckDBLoader) findTranscriptsByChrom(chrom string) ([]*Transcript, error) {
	rows, err := l.db.Query(`
		SELECT id, gene_id, gene_name, chrom, start, end_, strand, biotype,
		       is_canonical, cds_start, cds_end, cds_sequence, protein_sequence
		FROM transcripts
		WHERE chrom = ?
		ORDER BY start
	`, chrom)
	if err != nil {
		return nil, fmt.Errorf("query transcripts: %w", err)
	}
	defer rows.Close()

	var transcripts []*Transcript
	for rows.Next() {
		t, err := l.scanTranscript(rows)
		if err != nil {
			return nil, err
		}
		if err := l.loadExons(t); err != nil {
			return nil, err
		}
		transcripts = append(transcripts, t)
	}
	return transcripts, rows.Err()
}

// scanTranscript scans a transcript row into a Transcript struct.
func (l *DuckDBLoader) scanTranscript(rows *sql.Rows) (*Transcript, error) {
	t := &Transcript{}
	var cdsSeq, protSeq sql.NullString
	err := rows.Scan(
		&t.ID, &t.GeneID, &t.GeneName, &t.Chrom, &t.Start, &t.End,
		&t.Strand, &t.Biotype, &t.IsCanonical, &t.CDSStart, &t.CDSEnd,
		&cdsSeq, &protSeq,
	)
	if err != nil {
		return nil, fmt.Errorf("scan transcript: %w", err)
	}
	t.CDSSequence = cdsSeq.String
	t.ProteinSequence = protSeq.String
	return t, nil
}

// loadExons loads exons for a transcript.
func (l *DuckDBLoader) loadExons(t *Transcript) error {
	rows, err := l.db.Query(`
		SELECT exon_number, start, end_, cds_start, cds_end, frame
		FROM exons
		WHERE transcript_id = ?
		ORDER BY exon_number
	`, t.ID)
	if err != nil {
		return fmt.Errorf("query exons: %w", err)
	}
	defer rows.Close()

	for rows.Next() {
		var e Exon
		var cdsStart, cdsEnd sql.NullInt64
		var frame sql.NullInt64
		err := rows.Scan(&e.Number, &e.Start, &e.End, &cdsStart, &cdsEnd, &frame)
		if err != nil {
			return fmt.Errorf("scan exon: %w", err)
		}
		e.CDSStart = cdsStart.Int64
		e.CDSEnd = cdsEnd.Int64
		if frame.Valid {
			e.Frame = int(frame.Int64)
		} else {
			e.Frame = -1
		}
		t.Exons = append(t.Exons, e)
	}
	return rows.Err()
}

// CreateSchema creates the database schema for storing transcripts.
func (l *DuckDBLoader) CreateSchema() error {
	schema := `
		CREATE TABLE IF NOT EXISTS transcripts (
			id VARCHAR PRIMARY KEY,
			gene_id VARCHAR,
			gene_name VARCHAR,
			chrom VARCHAR,
			start BIGINT,
			end_ BIGINT,
			strand TINYINT,
			biotype VARCHAR,
			is_canonical BOOLEAN,
			cds_start BIGINT,
			cds_end BIGINT,
			cds_sequence VARCHAR,
			protein_sequence VARCHAR
		);

		CREATE TABLE IF NOT EXISTS exons (
			transcript_id VARCHAR,
			exon_number INTEGER,
			start BIGINT,
			end_ BIGINT,
			cds_start BIGINT,
			cds_end BIGINT,
			frame TINYINT,
			PRIMARY KEY (transcript_id, exon_number)
		);

		CREATE INDEX IF NOT EXISTS idx_transcripts_pos ON transcripts(chrom, start, end_);
		CREATE INDEX IF NOT EXISTS idx_transcripts_gene ON transcripts(gene_name);
		CREATE INDEX IF NOT EXISTS idx_exons_transcript ON exons(transcript_id);
	`
	_, err := l.db.Exec(schema)
	return err
}

// InsertTranscript inserts a transcript and its exons into the database.
func (l *DuckDBLoader) InsertTranscript(t *Transcript) error {
	_, err := l.db.Exec(`
		INSERT INTO transcripts (id, gene_id, gene_name, chrom, start, end_, strand,
		                         biotype, is_canonical, cds_start, cds_end,
		                         cds_sequence, protein_sequence)
		VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
	`, t.ID, t.GeneID, t.GeneName, t.Chrom, t.Start, t.End, t.Strand,
		t.Biotype, t.IsCanonical, t.CDSStart, t.CDSEnd,
		nullString(t.CDSSequence), nullString(t.ProteinSequence))
	if err != nil {
		return fmt.Errorf("insert transcript: %w", err)
	}

	for _, e := range t.Exons {
		var frame interface{}
		if e.Frame >= 0 {
			frame = e.Frame
		}
		_, err := l.db.Exec(`
			INSERT INTO exons (transcript_id, exon_number, start, end_, cds_start, cds_end, frame)
			VALUES (?, ?, ?, ?, ?, ?, ?)
		`, t.ID, e.Number, e.Start, e.End,
			nullInt64(e.CDSStart), nullInt64(e.CDSEnd), frame)
		if err != nil {
			return fmt.Errorf("insert exon: %w", err)
		}
	}
	return nil
}

// TranscriptCount returns the total number of transcripts in the database.
func (l *DuckDBLoader) TranscriptCount() (int, error) {
	var count int
	err := l.db.QueryRow("SELECT COUNT(*) FROM transcripts").Scan(&count)
	return count, err
}

// Chromosomes returns a sorted list of chromosomes in the database.
func (l *DuckDBLoader) Chromosomes() ([]string, error) {
	rows, err := l.db.Query("SELECT DISTINCT chrom FROM transcripts ORDER BY chrom")
	if err != nil {
		return nil, err
	}
	defer rows.Close()

	var chroms []string
	for rows.Next() {
		var chrom string
		if err := rows.Scan(&chrom); err != nil {
			return nil, err
		}
		chroms = append(chroms, chrom)
	}
	return chroms, rows.Err()
}

// ExportToJSON exports all transcripts to JSON format.
func (l *DuckDBLoader) ExportToJSON() ([]byte, error) {
	c := New()
	if err := l.LoadAll(c); err != nil {
		return nil, err
	}

	var all []*Transcript
	for _, chrom := range c.Chromosomes() {
		for _, t := range c.transcripts[chrom] {
			all = append(all, t)
		}
	}
	return json.MarshalIndent(all, "", "  ")
}

// nullString returns nil if s is empty, otherwise s.
func nullString(s string) interface{} {
	if s == "" {
		return nil
	}
	return s
}

// nullInt64 returns nil if n is 0, otherwise n.
func nullInt64(n int64) interface{} {
	if n == 0 {
		return nil
	}
	return n
}

// IsDuckDB checks if a path is a DuckDB database file.
func IsDuckDB(path string) bool {
	return strings.HasSuffix(path, ".duckdb") ||
		strings.HasSuffix(path, ".db") ||
		strings.HasPrefix(path, "s3://")
}
