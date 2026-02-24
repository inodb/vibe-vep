package output

import (
	"bytes"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

func TestTabWriter_WriteHeader(t *testing.T) {
	var buf bytes.Buffer
	w := NewTabWriter(&buf)

	require.NoError(t, w.WriteHeader())
	require.NoError(t, w.Flush())

	header := buf.String()

	// Check for expected columns
	expectedCols := []string{
		"#Uploaded_variation",
		"Location",
		"Allele",
		"Gene",
		"Feature",
		"Consequence",
		"IMPACT",
		"CANONICAL",
	}

	for _, col := range expectedCols {
		assert.Contains(t, header, col)
	}
}

func TestTabWriter_Write_KRASG12C(t *testing.T) {
	var buf bytes.Buffer
	w := NewTabWriter(&buf)

	v := &vcf.Variant{
		Chrom:  "12",
		Pos:    25245351,
		ID:     ".",
		Ref:    "C",
		Alt:    "A",
		Filter: "PASS",
	}

	ann := &annotate.Annotation{
		VariantID:       "12_25245351_C/A",
		TranscriptID:    "ENST00000311936",
		GeneName:        "KRAS",
		GeneID:          "ENSG00000133703",
		Consequence:     "missense_variant",
		Impact:          "MODERATE",
		CDSPosition:     34,
		ProteinPosition: 12,
		AminoAcidChange: "G12C",
		CodonChange:     "gGt/gCt",
		IsCanonical:     true,
		Allele:          "C",
		Biotype:         "protein_coding",
		ExonNumber:      "2/5",
	}

	require.NoError(t, w.WriteHeader())
	require.NoError(t, w.Write(v, ann))
	require.NoError(t, w.Flush())

	output := buf.String()
	lines := strings.Split(output, "\n")
	require.GreaterOrEqual(t, len(lines), 2)

	dataLine := lines[1]

	// Verify key fields
	checks := []struct {
		name  string
		value string
	}{
		{"location", "12:25245351"},
		{"gene", "KRAS"},
		{"transcript", "ENST00000311936"},
		{"consequence", "missense_variant"},
		{"impact", "MODERATE"},
		{"canonical", "YES"},
		{"amino acid change", "G12C"},
	}

	for _, check := range checks {
		assert.Contains(t, dataLine, check.value, check.name)
	}
}

func TestTabWriter_Write_IntergenicVariant(t *testing.T) {
	var buf bytes.Buffer
	w := NewTabWriter(&buf)

	v := &vcf.Variant{
		Chrom:  "1",
		Pos:    1000000,
		ID:     "rs123",
		Ref:    "A",
		Alt:    "G",
		Filter: "PASS",
	}

	ann := &annotate.Annotation{
		VariantID:   "1_1000000_A/G",
		Consequence: "intergenic_variant",
		Impact:      "MODIFIER",
		Allele:      "G",
	}

	require.NoError(t, w.Write(v, ann))
	require.NoError(t, w.Flush())

	output := buf.String()

	assert.Contains(t, output, "intergenic_variant")
	assert.Contains(t, output, "MODIFIER")
}
