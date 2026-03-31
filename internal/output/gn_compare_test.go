package output

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestParseGNAnnotation(t *testing.T) {
	input := `{"variant":"7:g.140753336A>T","originalVariantQuery":"7,140753336,140753336,A,T","hgvsg":"7:g.140753336A>T","assembly_name":"GRCh38","seq_region_name":"7","start":140753336,"end":140753336,"allele_string":"A/T","strand":1,"most_severe_consequence":"missense_variant","transcript_consequences":[{"transcript_id":"ENST00000288602","gene_symbol":"BRAF","gene_id":"ENSG00000157764","consequence_terms":["missense_variant"],"amino_acids":"V/E","protein_start":640,"protein_end":640,"hgvsp":"ENSP00000288602.7:p.Val640Glu","hgvsc":"ENST00000288602.11:c.1919T>A","canonical":"1","sift_score":0.0,"sift_prediction":"deleterious_low_confidence","polyphen_score":0.43,"polyphen_prediction":"benign"}],"successfully_annotated":true}`

	gn, err := ParseGNAnnotation([]byte(input))
	require.NoError(t, err)

	assert.Equal(t, "7:g.140753336A>T", gn.Variant)
	assert.Equal(t, "GRCh38", gn.AssemblyName)
	assert.Equal(t, "missense_variant", gn.MostSevereConsequence)
	assert.True(t, gn.SuccessfullyAnnotated)
	require.Len(t, gn.TranscriptConsequences, 1)

	tc := gn.TranscriptConsequences[0]
	assert.Equal(t, "ENST00000288602", tc.TranscriptID)
	assert.Equal(t, "BRAF", tc.GeneSymbol)
	assert.Equal(t, []string{"missense_variant"}, tc.ConsequenceTerms)
	assert.Equal(t, "V/E", tc.AminoAcids)
	assert.Equal(t, int64(640), tc.ProteinStart)
	assert.Equal(t, "1", tc.Canonical)
	require.NotNil(t, tc.SIFTScore)
	assert.InDelta(t, 0.0, *tc.SIFTScore, 0.001)
	assert.Equal(t, "deleterious_low_confidence", tc.SIFTPrediction)
}

func TestPickGNCanonical(t *testing.T) {
	tests := []struct {
		name string
		tcs  []GNTranscriptConsequence
		want string
	}{
		{
			name: "single canonical",
			tcs: []GNTranscriptConsequence{
				{TranscriptID: "ENST00000001", Canonical: ""},
				{TranscriptID: "ENST00000002", Canonical: "1", Biotype: "protein_coding"},
			},
			want: "ENST00000002",
		},
		{
			name: "no canonical, prefer protein_coding",
			tcs: []GNTranscriptConsequence{
				{TranscriptID: "ENST00000001", Biotype: "lncRNA"},
				{TranscriptID: "ENST00000002", Biotype: "protein_coding"},
			},
			want: "ENST00000002",
		},
		{
			name: "first if all equal",
			tcs: []GNTranscriptConsequence{
				{TranscriptID: "ENST00000001", Biotype: "lncRNA"},
				{TranscriptID: "ENST00000002", Biotype: "lncRNA"},
			},
			want: "ENST00000001",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := pickGNCanonical(tt.tcs)
			require.NotNil(t, got)
			assert.Equal(t, tt.want, got.TranscriptID)
		})
	}
}

func TestCompareGNToVEP(t *testing.T) {
	siftScore := 0.0
	ppScore := 0.43

	gn := &GNAnnotation{
		Variant:              "7:g.140753336A>T",
		OriginalVariantQuery: "7,140753336,140753336,A,T",
		TranscriptConsequences: []GNTranscriptConsequence{
			{
				TranscriptID:       "ENST00000288602",
				GeneSymbol:         "BRAF",
				ConsequenceTerms:   []string{"missense_variant"},
				AminoAcids:         "V/E",
				ProteinStart:       640,
				HGVSc:              "ENST00000288602.11:c.1919T>A",
				HGVSp:              "ENSP00000288602.7:p.Val640Glu",
				Canonical:          "1",
				Biotype:            "protein_coding",
				SIFTScore:          &siftScore,
				SIFTPrediction:     "deleterious_low_confidence",
				PolyPhenScore:      &ppScore,
				PolyPhenPrediction: "benign",
			},
		},
	}

	vepAnns := []*annotate.Annotation{
		{
			TranscriptID:    "ENST00000288602.11",
			GeneName:        "BRAF",
			Consequence:     "missense_variant",
			AminoAcidChange: "V640E",
			ProteinPosition: 640,
			HGVSc:           "c.1919T>A",
			HGVSp:           "p.Val640Glu",
			Extra: map[string]string{
				"sift.score":          "0",
				"sift.prediction":     "deleterious_low_confidence",
				"polyphen.score":      "0.43",
				"polyphen.prediction": "benign",
			},
		},
	}

	cmp := CompareGNToVEP(gn, vepAnns)

	assert.Equal(t, "7:g.140753336A>T", cmp.VariantID)
	assert.Equal(t, "ENST00000288602", cmp.TranscriptID)
	assert.Empty(t, cmp.Error)
	assert.True(t, cmp.AllMatch, "expected all fields to match, mismatches: %+v", cmp.Fields)

	for _, f := range cmp.Fields {
		if !f.Match && f.Category != CatBothEmpty {
			t.Errorf("field %s: GN=%q vibe-vep=%q category=%s", f.Field, f.GNValue, f.VEPValue, f.Category)
		}
	}
}

func TestCompareGNToVEP_NoTranscripts(t *testing.T) {
	gn := &GNAnnotation{Variant: "1:g.100A>T"}
	cmp := CompareGNToVEP(gn, nil)
	assert.Contains(t, cmp.Error, "no transcript consequences")

	gn.TranscriptConsequences = []GNTranscriptConsequence{{TranscriptID: "ENST00000001"}}
	cmp = CompareGNToVEP(gn, nil)
	assert.Contains(t, cmp.Error, "no annotations from vibe-vep")
}

func TestFloatScoresMatch(t *testing.T) {
	assert.True(t, floatScoresMatch("0", "0"))
	assert.True(t, floatScoresMatch("0.43", "0.43"))
	assert.True(t, floatScoresMatch("0.431", "0.433"))  // within 0.005
	assert.False(t, floatScoresMatch("0.43", "0.50"))
	assert.True(t, floatScoresMatch("", ""))
	assert.False(t, floatScoresMatch("0.5", ""))
}

func TestExtractFieldValues_HGVSc(t *testing.T) {
	gn := &GNTranscriptConsequence{
		HGVSc: "ENST00000288602.11:c.1919T>A",
	}
	vep := &annotate.Annotation{
		HGVSc: "c.1919T>A",
	}

	gnVal, vepVal := extractFieldValues("hgvsc", gn, vep)
	assert.Equal(t, "c.1919T>A", gnVal)
	assert.Equal(t, "c.1919T>A", vepVal)
}

func TestExtractFieldValues_HGVSp(t *testing.T) {
	gn := &GNTranscriptConsequence{
		HGVSp: "ENSP00000288602.7:p.Val640Glu",
	}
	vep := &annotate.Annotation{
		HGVSp: "p.Val640Glu",
	}

	gnVal, vepVal := extractFieldValues("hgvsp", gn, vep)
	assert.Equal(t, "p.V640E", gnVal)
	assert.Equal(t, "p.V640E", vepVal)
}

func TestGNComparisonReport(t *testing.T) {
	report := NewGNComparisonReport()

	// Add a full match.
	report.AddComparison(GNVariantComparison{
		VariantID: "1:g.100A>T",
		AllMatch:  true,
		Fields: []GNFieldComparison{
			{Field: "gene_symbol", GNValue: "TP53", VEPValue: "TP53", Match: true, Category: CatMatch},
		},
	})

	// Add a partial match.
	report.AddComparison(GNVariantComparison{
		VariantID: "2:g.200C>G",
		AllMatch:  false,
		Fields: []GNFieldComparison{
			{Field: "gene_symbol", GNValue: "KRAS", VEPValue: "KRAS", Match: true, Category: CatMatch},
			{Field: "consequence_terms", GNValue: "missense_variant", VEPValue: "synonymous_variant", Match: false, Category: CatMismatch},
		},
	})

	// Add an error.
	report.AddComparison(GNVariantComparison{
		VariantID: "3:g.300T>A",
		Error:     "annotation failed",
	})

	assert.Equal(t, 3, report.TotalVariants)
	assert.Equal(t, 1, report.FullMatches)
	assert.Equal(t, 1, report.PartialMatches)
	assert.Equal(t, 1, report.Errors)

	geneStats := report.FieldStats["gene_symbol"]
	require.NotNil(t, geneStats)
	assert.Equal(t, 2, geneStats.Total)
	assert.Equal(t, 2, geneStats.Matches)
	assert.InDelta(t, 100.0, geneStats.MatchRate(), 0.1)

	md := FormatGNReport(report)
	assert.Contains(t, md, "Total variants:** 3")
	assert.Contains(t, md, "gene_symbol")
}

func TestMatchVEPTranscript(t *testing.T) {
	anns := []*annotate.Annotation{
		{TranscriptID: "ENST00000001.5", GeneName: "A"},
		{TranscriptID: "ENST00000002.3", GeneName: "B"},
	}

	got := matchVEPTranscript("ENST00000002", anns)
	require.NotNil(t, got)
	assert.Equal(t, "B", got.GeneName)

	got = matchVEPTranscript("ENST00000999", anns)
	assert.Nil(t, got)
}
