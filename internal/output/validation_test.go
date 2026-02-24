package output

import (
	"bytes"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/vcf"
)

func TestValidationWriter_Match(t *testing.T) {
	var buf bytes.Buffer
	w := NewValidationWriter(&buf, true)

	require.NoError(t, w.WriteHeader())

	variant := &vcf.Variant{
		Chrom: "12",
		Pos:   25398284,
		Ref:   "C",
		Alt:   "A",
	}

	mafAnn := &maf.MAFAnnotation{
		HugoSymbol:   "KRAS",
		Consequence:  "missense_variant",
		HGVSpShort:   "p.G12C",
		TranscriptID: "ENST00000256078",
	}

	vepAnns := []*annotate.Annotation{
		{
			TranscriptID:    "ENST00000256078",
			GeneName:        "KRAS",
			Consequence:     "missense_variant",
			AminoAcidChange: "G/C",
			IsCanonical:     true,
		},
	}

	require.NoError(t, w.WriteComparison(variant, mafAnn, vepAnns))
	require.NoError(t, w.Flush())

	output := buf.String()
	assert.Contains(t, output, "KRAS")
	assert.Contains(t, output, "Y")

	total, matches, mismatches := w.Summary()
	assert.Equal(t, 1, total)
	assert.Equal(t, 1, matches)
	assert.Equal(t, 0, mismatches)
}

func TestValidationWriter_Mismatch(t *testing.T) {
	var buf bytes.Buffer
	w := NewValidationWriter(&buf, true)

	w.WriteHeader()

	variant := &vcf.Variant{
		Chrom: "12",
		Pos:   25398284,
		Ref:   "C",
		Alt:   "A",
	}

	mafAnn := &maf.MAFAnnotation{
		HugoSymbol:  "KRAS",
		Consequence: "missense_variant",
		HGVSpShort:  "p.G12C",
	}

	vepAnns := []*annotate.Annotation{
		{
			Consequence: "synonymous_variant", // Wrong prediction
		},
	}

	w.WriteComparison(variant, mafAnn, vepAnns)
	w.Flush()

	output := buf.String()
	assert.Contains(t, output, "N")

	_, matches, mismatches := w.Summary()
	assert.Equal(t, 0, matches)
	assert.Equal(t, 1, mismatches)
}

func TestNormalizeConsequence(t *testing.T) {
	tests := []struct {
		input    string
		expected string
	}{
		{"missense_variant", "missense_variant"},
		{"Missense_Mutation", "missense_variant"},
		{"MISSENSE_MUTATION", "missense_variant"},
		{"Nonsense_Mutation", "stop_gained"},
		{"stop_gained", "stop_gained"},
		{"Silent", "synonymous_variant"},
		{"Frame_Shift_Del", "frameshift_variant"},
		{"Frame_Shift_Ins", "frameshift_variant"},
		{"In_Frame_Del", "inframe_deletion"},
		{"unknown_type", "unknown_type"},
	}

	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			result := normalizeConsequence(tt.input)
			assert.Equal(t, tt.expected, result)
		})
	}
}

func TestValidationWriter_MismatchesOnly(t *testing.T) {
	var buf bytes.Buffer
	w := NewValidationWriter(&buf, false) // showAll = false

	w.WriteHeader()

	variant := &vcf.Variant{Chrom: "12", Pos: 100, Ref: "A", Alt: "T"}

	// Write a match - should not appear in output
	mafAnn := &maf.MAFAnnotation{Consequence: "missense_variant"}
	vepAnns := []*annotate.Annotation{{Consequence: "missense_variant"}}
	w.WriteComparison(variant, mafAnn, vepAnns)

	// Write a mismatch - should appear
	mafAnn2 := &maf.MAFAnnotation{Consequence: "stop_gained", HugoSymbol: "TP53"}
	vepAnns2 := []*annotate.Annotation{{Consequence: "missense_variant"}}
	w.WriteComparison(variant, mafAnn2, vepAnns2)

	w.Flush()

	output := buf.String()
	lines := strings.Split(strings.TrimSpace(output), "\n")

	// Should have header + 1 mismatch line only
	assert.Len(t, lines, 2)
	assert.Contains(t, output, "TP53")
}
