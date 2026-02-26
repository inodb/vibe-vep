package output

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/maf"
)

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
		{"frameshift_variant,stop_lost", "frameshift_variant"},
		{"intron_variant,splice_region_variant", "splice_region_variant"},
		{"missense_variant,splice_region_variant", "missense_variant"},
	}
	for _, tt := range tests {
		assert.Equal(t, tt.want, primaryConsequence(tt.input), "primaryConsequence(%q)", tt.input)
	}
}

func TestConsequencesMatch_SubAnnotations(t *testing.T) {
	assert.True(t, consequencesMatch("missense_variant", "missense_variant,NMD_transcript_variant"))
	assert.True(t, consequencesMatch("frameshift_variant", "frameshift_variant,splice_donor_5th_base_variant"))
	assert.False(t, consequencesMatch("missense_variant", "synonymous_variant"))
}

func TestIsNonStandardIntronicHGVSp(t *testing.T) {
	tests := []struct {
		input string
		want  bool
	}{
		{"p.*130*", true},
		{"p.*1*", true},
		{"p.*12345*", true},
		{"p.G12C", false},
		{"p.*130fs", false},
		{"p.*130=", false},
		{"", false},
		{"p.*", false},
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
	assert.True(t, hgvspValuesMatch("p.G12fs", "p.Gly12fs"))
	assert.True(t, hgvspValuesMatch("p.A10fs*33", "p.Ala15fs*40"))
	assert.False(t, hgvspValuesMatch("p.G12fs", "p.Gly12Cys"))
	assert.False(t, hgvspValuesMatch("p.G12C", "p.Gly12fs"))
	assert.True(t, hgvspValuesMatch("p.G12C", "p.Gly12Cys"))
	assert.True(t, hgvspValuesMatch("", ""))
}
