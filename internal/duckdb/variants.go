package duckdb

import (
	"context"
	"database/sql/driver"
	"fmt"

	goduckdb "github.com/marcboeker/go-duckdb"

	"github.com/inodb/vibe-vep/internal/annotate"
)

// variantKey is the composite key for variant cache lookups.
type variantKey struct {
	Chrom string
	Pos   int64
	Ref   string
	Alt   string
}

// VariantCache provides in-memory lookup of previously cached variant annotations.
type VariantCache struct {
	entries map[variantKey][]*annotate.Annotation
}

// Get returns cached annotations for a variant, or false if not cached.
func (vc *VariantCache) Get(chrom string, pos int64, ref, alt string) ([]*annotate.Annotation, bool) {
	key := variantKey{Chrom: chrom, Pos: pos, Ref: ref, Alt: alt}
	anns, ok := vc.entries[key]
	return anns, ok
}

// Len returns the number of cached variants.
func (vc *VariantCache) Len() int {
	return len(vc.entries)
}

// LoadVariantCache reads all variant results from DuckDB into an in-memory cache.
func (s *Store) LoadVariantCache() (*VariantCache, error) {
	rows, err := s.db.Query(`SELECT
		chrom, pos, ref, alt, transcript_id,
		gene_name, gene_id, consequence, impact,
		cds_position, protein_position, amino_acid_change, codon_change,
		is_canonical, allele, biotype, exon_number, intron_number,
		cdna_position, hgvsp, hgvsc, gene_type
		FROM variant_results`)
	if err != nil {
		return nil, fmt.Errorf("query variant_results: %w", err)
	}
	defer rows.Close()

	vc := &VariantCache{entries: make(map[variantKey][]*annotate.Annotation)}

	for rows.Next() {
		var chrom, ref, alt string
		var pos int64
		var ann annotate.Annotation

		if err := rows.Scan(
			&chrom, &pos, &ref, &alt, &ann.TranscriptID,
			&ann.GeneName, &ann.GeneID, &ann.Consequence, &ann.Impact,
			&ann.CDSPosition, &ann.ProteinPosition, &ann.AminoAcidChange, &ann.CodonChange,
			&ann.IsCanonical, &ann.Allele, &ann.Biotype, &ann.ExonNumber, &ann.IntronNumber,
			&ann.CDNAPosition, &ann.HGVSp, &ann.HGVSc, &ann.GeneType,
		); err != nil {
			return nil, fmt.Errorf("scan variant_result: %w", err)
		}

		ann.VariantID = annotate.FormatVariantID(chrom, pos, ref, alt)
		key := variantKey{Chrom: chrom, Pos: pos, Ref: ref, Alt: alt}
		a := ann // copy for pointer stability
		vc.entries[key] = append(vc.entries[key], &a)
	}

	if err := rows.Err(); err != nil {
		return nil, fmt.Errorf("iterate variant_results: %w", err)
	}

	return vc, nil
}

// VariantResult holds the data needed to write a variant result to DuckDB.
type VariantResult struct {
	Chrom string
	Pos   int64
	Ref   string
	Alt   string
	Ann   *annotate.Annotation
}

// WriteVariantResults batch-inserts variant results into DuckDB using the Appender API.
func (s *Store) WriteVariantResults(results []VariantResult) error {
	if len(results) == 0 {
		return nil
	}

	conn, err := s.db.Conn(context.Background())
	if err != nil {
		return fmt.Errorf("get connection: %w", err)
	}
	defer conn.Close()

	var appender *goduckdb.Appender
	if err := conn.Raw(func(driverConn any) error {
		var err error
		appender, err = goduckdb.NewAppenderFromConn(driverConn.(driver.Conn), "", "variant_results")
		return err
	}); err != nil {
		return fmt.Errorf("create appender: %w", err)
	}
	defer appender.Close()

	for _, r := range results {
		a := r.Ann
		if err := appender.AppendRow(
			r.Chrom, r.Pos, r.Ref, r.Alt, a.TranscriptID,
			a.GeneName, a.GeneID, a.Consequence, a.Impact,
			a.CDSPosition, a.ProteinPosition, a.AminoAcidChange, a.CodonChange,
			a.IsCanonical, a.Allele, a.Biotype, a.ExonNumber, a.IntronNumber,
			a.CDNAPosition, a.HGVSp, a.HGVSc, a.GeneType,
		); err != nil {
			return fmt.Errorf("append variant result: %w", err)
		}
	}

	return appender.Flush()
}

// ClearVariantResults removes all cached variant results.
func (s *Store) ClearVariantResults() error {
	_, err := s.db.Exec("DELETE FROM variant_results")
	return err
}
