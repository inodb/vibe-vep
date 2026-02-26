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

func TestTranscriptBaseID(t *testing.T) {
	tests := []struct {
		input string
		want  string
	}{
		{"ENST00000333418.4", "ENST00000333418"},
		{"ENST00000333418.5", "ENST00000333418"},
		{"ENST00000333418", "ENST00000333418"},
		{"", ""},
	}
	for _, tt := range tests {
		assert.Equal(t, tt.want, transcriptBaseID(tt.input), "transcriptBaseID(%q)", tt.input)
	}
}

func TestSelectBestAnnotation_VersionMismatch(t *testing.T) {
	// MAF has transcript v4, VEP has v5 — should still match
	mafAnn := &maf.MAFAnnotation{
		HugoSymbol:   "KRAS",
		Consequence:  "missense_variant",
		TranscriptID: "ENST00000256078.4",
	}
	vepAnns := []*annotate.Annotation{
		{
			TranscriptID: "ENST00000311936.8",
			GeneName:     "KRAS",
			Consequence:  "missense_variant",
			Biotype:      "protein_coding",
			IsCanonical:  true,
		},
		{
			TranscriptID: "ENST00000256078.5",
			GeneName:     "KRAS",
			Consequence:  "missense_variant",
			Biotype:      "protein_coding",
		},
	}

	best := SelectBestAnnotation(mafAnn, vepAnns)
	require.NotNil(t, best)
	assert.Equal(t, "ENST00000256078.5", best.TranscriptID)
}

func TestPrimaryConsequence(t *testing.T) {
	tests := []struct {
		input string
		want  string
	}{
		{"missense_variant", "missense_variant"},
		{"frameshift_variant,stop_lost", "frameshift_variant"}, // both HIGH, first wins
		{"intron_variant,splice_region_variant", "splice_region_variant"},
		{"missense_variant,splice_region_variant", "missense_variant"},
	}
	for _, tt := range tests {
		assert.Equal(t, tt.want, primaryConsequence(tt.input), "primaryConsequence(%q)", tt.input)
	}
}

func TestConsequencesMatch_SubAnnotations(t *testing.T) {
	// Same primary term but different sub-annotations should match
	assert.True(t, consequencesMatch("missense_variant", "missense_variant,NMD_transcript_variant"))
	assert.True(t, consequencesMatch("frameshift_variant", "frameshift_variant,splice_donor_5th_base_variant"))
	// Different primary terms should not match
	assert.False(t, consequencesMatch("missense_variant", "synonymous_variant"))
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

func TestIsNonStandardIntronicHGVSp(t *testing.T) {
	tests := []struct {
		input string
		want  bool
	}{
		{"p.*130*", true},     // typical MAF intronic notation
		{"p.*1*", true},       // single digit
		{"p.*12345*", true},   // large number
		{"p.G12C", false},     // missense
		{"p.*130fs", false},   // frameshift starting at stop
		{"p.*130=", false},    // synonymous stop
		{"", false},           // empty
		{"p.*", false},        // just p.* (3 chars, prefix+suffix overlap)
	}
	for _, tt := range tests {
		assert.Equal(t, tt.want, isNonStandardIntronicHGVSp(tt.input), "isNonStandardIntronicHGVSp(%q)", tt.input)
	}
}

func TestIsFrameshiftHGVSp(t *testing.T) {
	assert.True(t, isFrameshiftHGVSp("p.G12fs"))
	assert.True(t, isFrameshiftHGVSp("p.G12Afs*33"))
	assert.False(t, isFrameshiftHGVSp("p.G12C"))
	assert.False(t, isFrameshiftHGVSp(""))
}

func TestIsSpliceConsequence(t *testing.T) {
	assert.True(t, isSpliceConsequence("splice_donor_variant"))
	assert.True(t, isSpliceConsequence("splice_acceptor_variant"))
	assert.True(t, isSpliceConsequence("splice_donor_variant,intron_variant"))
	assert.False(t, isSpliceConsequence("splice_region_variant"))
	assert.False(t, isSpliceConsequence("missense_variant"))
}

func TestHgvspValuesMatch_FrameshiftFuzzy(t *testing.T) {
	// Both frameshifts with different positions → match (GENCODE version shift)
	assert.True(t, hgvspValuesMatch("p.G12fs", "p.Gly12fs"))
	assert.True(t, hgvspValuesMatch("p.A10fs*33", "p.Ala15fs*40"))
	// One frameshift, one not → no match
	assert.False(t, hgvspValuesMatch("p.G12fs", "p.Gly12Cys"))
	assert.False(t, hgvspValuesMatch("p.G12C", "p.Gly12fs"))
	// Exact match still works
	assert.True(t, hgvspValuesMatch("p.G12C", "p.Gly12Cys"))
	// Both empty
	assert.True(t, hgvspValuesMatch("", ""))
}

func TestWriteComparison_SkipsNonStandardIntronic(t *testing.T) {
	var buf bytes.Buffer
	w := NewValidationWriter(&buf, false)
	w.WriteHeader()

	variant := &vcf.Variant{Chrom: "1", Pos: 100, Ref: "A", Alt: "T"}
	mafAnn := &maf.MAFAnnotation{
		HugoSymbol:  "TP53",
		Consequence: "intron",
		HGVSpShort:  "p.*130*",
	}
	vepAnns := []*annotate.Annotation{{
		Consequence: "intron_variant",
		HGVSp:      "",
	}}

	w.WriteComparison(variant, mafAnn, vepAnns)
	w.Flush()

	_, mismatches, skipped := w.HGVSpSummary()
	assert.Equal(t, 0, mismatches)
	assert.Equal(t, 1, skipped)
}

func TestWriteComparison_SkipsSpliceEmptyHGVSp(t *testing.T) {
	var buf bytes.Buffer
	w := NewValidationWriter(&buf, false)
	w.WriteHeader()

	variant := &vcf.Variant{Chrom: "1", Pos: 100, Ref: "ACGT", Alt: "A"}
	mafAnn := &maf.MAFAnnotation{
		HugoSymbol:  "TP53",
		Consequence: "Splice_Site",
		HGVSpShort:  "p.X125_splice",
	}
	vepAnns := []*annotate.Annotation{{
		Consequence: "splice_acceptor_variant,intron_variant",
		HGVSp:      "",
	}}

	w.WriteComparison(variant, mafAnn, vepAnns)
	w.Flush()

	_, mismatches, skipped := w.HGVSpSummary()
	assert.Equal(t, 0, mismatches)
	assert.Equal(t, 1, skipped)
}
