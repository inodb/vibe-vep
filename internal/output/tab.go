// Package output provides annotation output formatters.
package output

import (
	"bufio"
	"fmt"
	"io"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// TabWriter writes annotations in tab-delimited format.
type TabWriter struct {
	w       *bufio.Writer
	columns []string
}

// NewTabWriter creates a new tab-delimited writer.
func NewTabWriter(w io.Writer) *TabWriter {
	return &TabWriter{
		w: bufio.NewWriter(w),
		columns: []string{
			"#Uploaded_variation",
			"Location",
			"Allele",
			"Gene",
			"Feature",
			"Feature_type",
			"Consequence",
			"cDNA_position",
			"CDS_position",
			"Protein_position",
			"Amino_acids",
			"Codons",
			"Existing_variation",
			"IMPACT",
			"BIOTYPE",
			"CANONICAL",
			"EXON",
			"INTRON",
			"HGVSp",
		},
	}
}

// WriteHeader writes the header line.
func (tw *TabWriter) WriteHeader() error {
	_, err := tw.w.WriteString(strings.Join(tw.columns, "\t") + "\n")
	return err
}

// Write writes a single annotation.
func (tw *TabWriter) Write(v *vcf.Variant, ann *annotate.Annotation) error {
	// Format location
	location := fmt.Sprintf("%s:%d", v.Chrom, v.Pos)

	// Format amino acids (if present)
	aminoAcids := "-"
	if ann.AminoAcidChange != "" {
		// Extract ref and alt AA from change like "G12C"
		aminoAcids = ann.AminoAcidChange
	}

	// Format codons (if present)
	codons := "-"
	if ann.CodonChange != "" {
		codons = ann.CodonChange
	}

	// Format canonical
	canonical := "-"
	if ann.IsCanonical {
		canonical = "YES"
	}

	// Format positions
	cdsPos := "-"
	if ann.CDSPosition > 0 {
		cdsPos = fmt.Sprintf("%d", ann.CDSPosition)
	}

	proteinPos := "-"
	if ann.ProteinPosition > 0 {
		proteinPos = fmt.Sprintf("%d", ann.ProteinPosition)
	}

	cdnaPos := "-"
	if ann.CDNAPosition > 0 {
		cdnaPos = fmt.Sprintf("%d", ann.CDNAPosition)
	}

	// Feature type
	featureType := "Transcript"
	if ann.TranscriptID == "" {
		featureType = "-"
	}

	// Gene and feature
	gene := ann.GeneName
	if gene == "" {
		gene = "-"
	}
	feature := ann.TranscriptID
	if feature == "" {
		feature = "-"
	}

	// Biotype
	biotype := ann.Biotype
	if biotype == "" {
		biotype = "-"
	}

	// Exon/Intron
	exon := ann.ExonNumber
	if exon == "" {
		exon = "-"
	}
	intron := ann.IntronNumber
	if intron == "" {
		intron = "-"
	}
	hgvsp := ann.HGVSp
	if hgvsp == "" {
		hgvsp = "-"
	}

	// Build the row
	values := []string{
		v.ID,
		location,
		ann.Allele,
		gene,
		feature,
		featureType,
		ann.Consequence,
		cdnaPos,
		cdsPos,
		proteinPos,
		aminoAcids,
		codons,
		"-", // Existing_variation (not implemented)
		ann.Impact,
		biotype,
		canonical,
		exon,
		intron,
		hgvsp,
	}

	_, err := tw.w.WriteString(strings.Join(values, "\t") + "\n")
	return err
}

// Flush flushes any buffered data to the underlying writer.
func (tw *TabWriter) Flush() error {
	return tw.w.Flush()
}
