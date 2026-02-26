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

func TestCategorizeConsequence(t *testing.T) {
	tests := []struct {
		name     string
		maf, vep string
		want     Category
	}{
		{"exact match", "missense_variant", "missense_variant", CatMatch},
		{"normalized match", "Missense_Mutation", "missense_variant", CatMatch},
		{"primary term match", "missense_variant", "missense_variant,NMD_transcript_variant", CatMatch},
		{"upstream reclassified", "5'Flank", "5_prime_utr_variant", CatUpstreamReclass},
		{"downstream reclassified", "3'Flank", "intron_variant", CatUpstreamReclass},
		{"mismatch", "missense_variant", "synonymous_variant", CatMismatch},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			assert.Equal(t, tt.want, categorizeConsequence(tt.maf, tt.vep))
		})
	}
}

func TestCategorizeHGVSp(t *testing.T) {
	tests := []struct {
		name                             string
		mafHGVSp, vepHGVSp, vepConseq   string
		mafHGVSc, vepHGVSc              string
		want                             Category
	}{
		{
			name: "exact match",
			mafHGVSp: "p.G12C", vepHGVSp: "p.G12C",
			want: CatMatch,
		},
		{
			name: "both empty",
			mafHGVSp: "", vepHGVSp: "",
			want: CatBothEmpty,
		},
		{
			name: "fuzzy frameshift",
			mafHGVSp: "p.G12fs", vepHGVSp: "p.G15fs*33",
			want: CatFuzzyFS,
		},
		{
			name: "splice vs syn",
			mafHGVSp: "p.X125_splice", vepHGVSp: "p.G125=",
			want: CatSpliceVsSyn,
		},
		{
			name: "maf nonstandard intronic",
			mafHGVSp: "p.*130*", vepHGVSp: "",
			want: CatMafNonstandard,
		},
		{
			name: "splice no protein",
			mafHGVSp: "p.X125_splice", vepHGVSp: "",
			vepConseq: "splice_acceptor_variant",
			want: CatSpliceNoProtein,
		},
		{
			name: "maf empty",
			mafHGVSp: "", vepHGVSp: "p.G12C",
			want: CatMafEmpty,
		},
		{
			name: "vep empty",
			mafHGVSp: "p.G12C", vepHGVSp: "",
			want: CatVepEmpty,
		},
		{
			name: "position shift - same change type",
			mafHGVSp: "p.Y1145C", vepHGVSp: "p.Y1215C",
			want: CatPositionShift,
		},
		{
			name: "position shift - synonymous",
			mafHGVSp: "p.A735=", vepHGVSp: "p.A744=",
			want: CatPositionShift,
		},
		{
			name: "position shift - hgvsc matches",
			mafHGVSp: "p.G12C", vepHGVSp: "p.G15C",
			mafHGVSc: "ENST00000256078.4:c.34G>T", vepHGVSc: "c.34G>T",
			want: CatPositionShift,
		},
		{
			name: "mismatch",
			mafHGVSp: "p.G12C", vepHGVSp: "p.A15T",
			mafHGVSc: "ENST00000256078.4:c.34G>T", vepHGVSc: "c.99A>C",
			want: CatMismatch,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := categorizeHGVSp(tt.mafHGVSp, tt.vepHGVSp, tt.vepConseq, tt.mafHGVSc, tt.vepHGVSc)
			assert.Equal(t, tt.want, got)
		})
	}
}

func TestCategorizeHGVSc(t *testing.T) {
	tests := []struct {
		name     string
		maf, vep string
		want     Category
	}{
		{"match with prefix", "ENST00000361923.2:c.1428C>G", "c.1428C>G", CatMatch},
		{"both empty", "", "", CatBothEmpty},
		{"maf empty", "", "c.1428C>G", CatMafEmpty},
		{"vep empty", "ENST00000361923.2:c.1428C>G", "", CatVepEmpty},
		{"mismatch", "ENST00000361923.2:c.1428C>G", "c.999A>T", CatMismatch},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			assert.Equal(t, tt.want, categorizeHGVSc(tt.maf, tt.vep))
		})
	}
}

func TestCompareWriter_AllColumns(t *testing.T) {
	var buf bytes.Buffer
	cols := map[string]bool{"consequence": true, "hgvsp": true, "hgvsc": true}
	w := NewCompareWriter(&buf, cols, true)

	require.NoError(t, w.WriteHeader())

	variant := &vcf.Variant{Chrom: "12", Pos: 25398284, Ref: "C", Alt: "A"}
	mafAnn := &maf.MAFAnnotation{
		HugoSymbol:   "KRAS",
		Consequence:  "missense_variant",
		HGVSpShort:   "p.G12C",
		TranscriptID: "ENST00000256078",
	}
	vepAnns := []*annotate.Annotation{{
		TranscriptID: "ENST00000256078",
		GeneName:     "KRAS",
		Consequence:  "missense_variant",
		HGVSp:        "p.Gly12Cys",
		IsCanonical:  true,
	}}

	require.NoError(t, w.WriteComparison(variant, mafAnn, vepAnns))
	require.NoError(t, w.Flush())

	out := buf.String()
	assert.Contains(t, out, "KRAS")
	assert.Contains(t, out, "match")
	assert.Equal(t, 1, w.Total())
	assert.Equal(t, 1, w.Counts()["consequence"][CatMatch])
}

func TestCompareWriter_DefaultHideMatches(t *testing.T) {
	var buf bytes.Buffer
	cols := map[string]bool{"consequence": true, "hgvsp": true, "hgvsc": true}
	w := NewCompareWriter(&buf, cols, false) // showAll = false

	w.WriteHeader()

	variant := &vcf.Variant{Chrom: "12", Pos: 100, Ref: "A", Alt: "T"}

	// Match row — should be hidden
	mafAnn := &maf.MAFAnnotation{Consequence: "missense_variant"}
	vepAnns := []*annotate.Annotation{{Consequence: "missense_variant"}}
	w.WriteComparison(variant, mafAnn, vepAnns)

	// Mismatch row — should be shown
	mafAnn2 := &maf.MAFAnnotation{Consequence: "stop_gained", HugoSymbol: "TP53"}
	vepAnns2 := []*annotate.Annotation{{Consequence: "missense_variant"}}
	w.WriteComparison(variant, mafAnn2, vepAnns2)

	w.Flush()

	out := buf.String()
	lines := strings.Split(strings.TrimSpace(out), "\n")
	// Header + 1 mismatch line
	assert.Len(t, lines, 2)
	assert.Contains(t, out, "TP53")
}

func TestCompareWriter_Summary(t *testing.T) {
	var buf, summary bytes.Buffer
	cols := map[string]bool{"consequence": true}
	w := NewCompareWriter(&buf, cols, true)

	variant := &vcf.Variant{Chrom: "1", Pos: 1, Ref: "A", Alt: "T"}
	w.WriteComparison(variant, &maf.MAFAnnotation{Consequence: "missense_variant"},
		[]*annotate.Annotation{{Consequence: "missense_variant"}})
	w.WriteComparison(variant, &maf.MAFAnnotation{Consequence: "stop_gained"},
		[]*annotate.Annotation{{Consequence: "missense_variant"}})

	w.WriteSummary(&summary)

	out := summary.String()
	assert.Contains(t, out, "Comparison Summary (2 variants)")
	assert.Contains(t, out, "Consequence:")
	assert.Contains(t, out, "match")
	assert.Contains(t, out, "mismatch")
}
