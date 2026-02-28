package duckdb

import (
	"context"
	"database/sql/driver"
	"fmt"

	goduckdb "github.com/marcboeker/go-duckdb"

	"github.com/inodb/vibe-vep/internal/annotate"
)

// VariantResult holds the data needed to write a variant result to DuckDB.
type VariantResult struct {
	Chrom string
	Pos   int64
	Ref   string
	Alt   string
	Ann   *annotate.Annotation
}

// resultKey is the composite key for deduplicating variant results before writing.
type resultKey struct {
	chrom, ref, alt, transcriptID string
	pos                           int64
}

// WriteVariantResults batch-inserts variant results into DuckDB using the Appender API.
// Duplicate (chrom, pos, ref, alt, transcript_id) entries are deduplicated before writing.
func (s *Store) WriteVariantResults(results []VariantResult) error {
	if len(results) == 0 {
		return nil
	}

	// Deduplicate by primary key (same variant from multiple MAF rows)
	seen := make(map[resultKey]bool, len(results))
	deduped := make([]VariantResult, 0, len(results))
	for _, r := range results {
		k := resultKey{r.Chrom, r.Ref, r.Alt, r.Ann.TranscriptID, r.Pos}
		if !seen[k] {
			seen[k] = true
			deduped = append(deduped, r)
		}
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

	for _, r := range deduped {
		a := r.Ann
		if err := appender.AppendRow(
			r.Chrom, r.Pos, r.Ref, r.Alt, a.TranscriptID,
			a.GeneName, a.GeneID, a.Consequence, a.Impact,
			a.CDSPosition, a.ProteinPosition, a.AminoAcidChange, a.CodonChange,
			a.IsCanonical, a.Allele, a.Biotype, a.ExonNumber, a.IntronNumber,
			a.CDNAPosition, a.HGVSp, a.HGVSc, a.GeneType,
			float32(a.AlphaMissenseScore), a.AlphaMissenseClass,
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

// LookupVariant queries DuckDB for previously cached annotations of a specific variant.
func (s *Store) LookupVariant(chrom string, pos int64, ref, alt string) ([]*annotate.Annotation, error) {
	rows, err := s.db.Query(`SELECT
		transcript_id, gene_name, gene_id, consequence, impact,
		cds_position, protein_position, amino_acid_change, codon_change,
		is_canonical, allele, biotype, exon_number, intron_number,
		cdna_position, hgvsp, hgvsc, gene_type,
		am_score, am_class
		FROM variant_results
		WHERE chrom=? AND pos=? AND ref=? AND alt=?`,
		chrom, pos, ref, alt)
	if err != nil {
		return nil, fmt.Errorf("query variant: %w", err)
	}
	defer rows.Close()

	var anns []*annotate.Annotation
	for rows.Next() {
		var ann annotate.Annotation
		if err := rows.Scan(
			&ann.TranscriptID, &ann.GeneName, &ann.GeneID, &ann.Consequence, &ann.Impact,
			&ann.CDSPosition, &ann.ProteinPosition, &ann.AminoAcidChange, &ann.CodonChange,
			&ann.IsCanonical, &ann.Allele, &ann.Biotype, &ann.ExonNumber, &ann.IntronNumber,
			&ann.CDNAPosition, &ann.HGVSp, &ann.HGVSc, &ann.GeneType,
			&ann.AlphaMissenseScore, &ann.AlphaMissenseClass,
		); err != nil {
			return nil, fmt.Errorf("scan variant: %w", err)
		}
		ann.VariantID = annotate.FormatVariantID(chrom, pos, ref, alt)
		anns = append(anns, &ann)
	}
	if err := rows.Err(); err != nil {
		return nil, fmt.Errorf("iterate variants: %w", err)
	}
	return anns, nil
}

// SearchByGene queries DuckDB for all cached variant results for a gene.
func (s *Store) SearchByGene(geneName string) ([]VariantResult, error) {
	rows, err := s.db.Query(`SELECT
		chrom, pos, ref, alt, transcript_id,
		gene_name, gene_id, consequence, impact,
		cds_position, protein_position, amino_acid_change, codon_change,
		is_canonical, allele, biotype, exon_number, intron_number,
		cdna_position, hgvsp, hgvsc, gene_type,
		am_score, am_class
		FROM variant_results
		WHERE gene_name=?`, geneName)
	if err != nil {
		return nil, fmt.Errorf("query by gene: %w", err)
	}
	defer rows.Close()

	return scanVariantResults(rows)
}

// SearchByProteinChange queries DuckDB for cached results matching a specific
// gene and amino acid change (e.g. "G12C").
func (s *Store) SearchByProteinChange(geneName, aaChange string) ([]VariantResult, error) {
	rows, err := s.db.Query(`SELECT
		chrom, pos, ref, alt, transcript_id,
		gene_name, gene_id, consequence, impact,
		cds_position, protein_position, amino_acid_change, codon_change,
		is_canonical, allele, biotype, exon_number, intron_number,
		cdna_position, hgvsp, hgvsc, gene_type,
		am_score, am_class
		FROM variant_results
		WHERE gene_name=? AND amino_acid_change=?`, geneName, aaChange)
	if err != nil {
		return nil, fmt.Errorf("query by protein change: %w", err)
	}
	defer rows.Close()

	return scanVariantResults(rows)
}

// scanVariantResults scans rows into VariantResult slices.
func scanVariantResults(rows interface {
	Next() bool
	Scan(dest ...any) error
	Err() error
}) ([]VariantResult, error) {
	var results []VariantResult
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
			&ann.AlphaMissenseScore, &ann.AlphaMissenseClass,
		); err != nil {
			return nil, fmt.Errorf("scan variant result: %w", err)
		}

		ann.VariantID = annotate.FormatVariantID(chrom, pos, ref, alt)
		a := ann
		results = append(results, VariantResult{
			Chrom: chrom, Pos: pos, Ref: ref, Alt: alt, Ann: &a,
		})
	}
	if err := rows.Err(); err != nil {
		return nil, fmt.Errorf("iterate variant results: %w", err)
	}
	return results, nil
}
