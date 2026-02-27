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
		{"no cds data", "synonymous_variant", "coding_sequence_variant", CatNoCDS},
		{"frameshift+splice_region → splice_acceptor", "frameshift_variant,splice_region_variant", "splice_acceptor_variant", CatMatch},
		{"frameshift+splice_region → splice_donor", "frameshift_variant,splice_region_variant", "splice_donor_variant", CatMatch},
		{"splice_region+intron → splice_acceptor", "splice_region_variant,intron_variant", "splice_acceptor_variant", CatMatch},
		{"splice_region+intron → splice_donor", "splice_region_variant,intron_variant", "splice_donor_variant", CatMatch},
		{"splice_region+polypyrimidine → splice_acceptor", "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant", "splice_acceptor_variant", CatMatch},
		{"inframe_deletion ↔ stop_gained", "stop_gained,inframe_deletion", "inframe_deletion", CatMatch},
		{"inframe_insertion ↔ stop_gained", "inframe_insertion", "stop_gained", CatMatch},
		{"synonymous ↔ stop_retained", "synonymous_variant", "stop_retained_variant", CatMatch},
		{"stop_retained ↔ synonymous", "stop_retained_variant", "synonymous_variant", CatMatch},
		{"compound match - vep primary in maf", "frameshift_variant,start_lost", "start_lost", CatMatch},
		{"coding → non_coding_exon biotype change", "missense_variant", "non_coding_transcript_exon_variant", CatTranscriptModelChange},
		{"synonymous → non_coding_exon biotype change", "synonymous_variant", "non_coding_transcript_exon_variant", CatTranscriptModelChange},
		{"5'UTR → intron boundary shift", "5_prime_UTR_variant", "intron_variant", CatUpstreamReclass},
		// Phase 1a: normalize modifiers
		{"drop start_retained with frameshift", "frameshift_variant,start_retained_variant", "start_lost", CatMatch},
		{"drop stop_retained with stop_gained", "stop_gained,stop_retained_variant", "stop_gained", CatMatch},
		{"drop coding_sequence_variant with splice HIGH", "splice_acceptor_variant,coding_sequence_variant,intron_variant", "splice_acceptor_variant", CatMatch},
		{"drop coding_sequence_variant with splice donor", "splice_donor_variant,coding_sequence_variant,intron_variant", "splice_donor_variant", CatMatch},
		// Phase 1b: consequence equivalences
		{"inframe ↔ stop_lost", "inframe_deletion", "stop_lost", CatMatch},
		{"stop_lost ↔ inframe", "stop_lost", "inframe_insertion", CatMatch},
		{"stop_lost ↔ stop_retained", "stop_lost", "stop_retained_variant", CatMatch},
		{"stop_retained ↔ stop_lost", "stop_retained_variant", "stop_lost", CatMatch},
		{"start_lost ↔ synonymous", "start_lost", "synonymous_variant", CatMatch},
		{"synonymous ↔ start_lost", "synonymous_variant", "start_lost", CatMatch},
		{"start_lost ↔ missense", "start_lost", "missense_variant", CatMatch},
		{"start_lost ↔ inframe", "start_lost", "inframe_deletion", CatMatch},
		{"inframe ↔ start_lost", "In_Frame_Ins", "start_lost", CatMatch},
		// Phase 2: transcript model change
		{"coding → non_coding_exon transcript model", "stop_gained", "non_coding_transcript_exon_variant", CatTranscriptModelChange},
		{"coding → intron non-coding", "missense_variant", "intron_variant,non_coding_transcript_exon_variant", CatTranscriptModelChange},
		{"coding → intron transcript model", "frameshift_variant", "intron_variant", CatTranscriptModelChange},
		// Phase 2: gene model change (coding ↔ intergenic)
		{"missense → intergenic gene model", "missense_variant", "intergenic_variant", CatGeneModelChange},
		{"synonymous → intergenic gene model", "synonymous_variant", "intergenic_variant", CatGeneModelChange},
		{"intergenic → missense gene model", "intergenic_variant", "missense_variant", CatGeneModelChange},
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
			name: "frameshift vs stop_gained",
			mafHGVSp: "p.E454Sfs*117", vepHGVSp: "p.K453*",
			want: CatFuzzyFS,
		},
		{
			name: "stop_gained vs frameshift",
			mafHGVSp: "p.K453*", vepHGVSp: "p.E454Sfs*117",
			want: CatFuzzyFS,
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
		{"position shift SNV", "ENST00000333418.4:c.3434A>G", "c.1145A>G", CatPositionShift},
		{"position shift del", "ENST00000295754.10:c.94+16291_94+16307del", "c.94+16290_94+16306del", CatPositionShift},
		{"position shift dup", "ENST00000357077.9:c.4426dup", "c.4427dup", CatPositionShift},
		{"delins normalized", "ENST00000324385.9:c.1382_1407delinsG", "c.1383_1407del", CatDelinsNorm},
		{"dup vs ins", "ENST00000397901.8:c.*46_*49dup", "c.*49_*50insAAGG", CatDupVsIns},
		{"mismatch different SNV", "ENST00000361923.2:c.1428C>G", "c.999A>T", CatMismatch},
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

func TestCompareWriter_CrossColumnPositionShift(t *testing.T) {
	// When HGVSc shows a position shift, a consequence mismatch should be
	// reclassified as position_shift since different CDS positions explain
	// the consequence difference (GENCODE version change).
	var buf bytes.Buffer
	cols := map[string]bool{"consequence": true, "hgvsc": true}
	w := NewCompareWriter(&buf, cols, true)

	variant := &vcf.Variant{Chrom: "22", Pos: 38026136, Ref: "T", Alt: "G"}

	// MAF says synonymous at c.390, VEP says missense at c.388 (position shift)
	mafAnn := &maf.MAFAnnotation{
		Consequence:  "synonymous_variant",
		HGVSc:        "ENST00000333418.4:c.390T>G",
		TranscriptID: "ENST00000333418",
	}
	vepAnns := []*annotate.Annotation{{
		TranscriptID: "ENST00000333418",
		Consequence:  "missense_variant",
		HGVSc:        "c.388T>G",
	}}

	w.WriteComparison(variant, mafAnn, vepAnns)

	counts := w.Counts()
	assert.Equal(t, 1, counts["consequence"][CatPositionShift],
		"consequence mismatch with HGVSc position_shift should be reclassified")
	assert.Equal(t, 0, counts["consequence"][CatMismatch],
		"should have no consequence mismatches")
}

func TestCompareWriter_CrossColumnHGVSpPositionShift(t *testing.T) {
	// When HGVSc shows a position shift, an HGVSp mismatch should be
	// reclassified as position_shift since different CDS positions explain
	// the amino acid difference.
	var buf bytes.Buffer
	cols := map[string]bool{"hgvsp": true, "hgvsc": true}
	w := NewCompareWriter(&buf, cols, true)

	variant := &vcf.Variant{Chrom: "22", Pos: 38026136, Ref: "T", Alt: "G"}

	mafAnn := &maf.MAFAnnotation{
		HGVSpShort:   "p.P130=",
		HGVSc:        "ENST00000333418.4:c.390T>G",
		TranscriptID: "ENST00000333418",
	}
	vepAnns := []*annotate.Annotation{{
		TranscriptID: "ENST00000333418",
		HGVSp:        "p.F130V",
		HGVSc:        "c.388T>G",
	}}

	w.WriteComparison(variant, mafAnn, vepAnns)

	counts := w.Counts()
	assert.Equal(t, 1, counts["hgvsp"][CatPositionShift],
		"hgvsp mismatch with HGVSc position_shift should be reclassified")
	assert.Equal(t, 0, counts["hgvsp"][CatMismatch],
		"should have no hgvsp mismatches")
}

func TestCompareWriter_CrossColumnConseqDelinsNorm(t *testing.T) {
	// When HGVSc is delins_normalized, a consequence mismatch should be
	// reclassified as delins_normalized (MNV missense→synonymous).
	var buf bytes.Buffer
	cols := map[string]bool{"consequence": true, "hgvsc": true}
	w := NewCompareWriter(&buf, cols, true)

	variant := &vcf.Variant{Chrom: "1", Pos: 100, Ref: "ACG", Alt: "GTA"}

	mafAnn := &maf.MAFAnnotation{
		Consequence:  "Missense_Mutation",
		HGVSc:        "ENST00000123456.1:c.100_102delinsGTA",
		TranscriptID: "ENST00000123456",
	}
	vepAnns := []*annotate.Annotation{{
		TranscriptID: "ENST00000123456",
		Consequence:  "synonymous_variant",
		HGVSc:        "c.102G>A",
	}}

	w.WriteComparison(variant, mafAnn, vepAnns)

	counts := w.Counts()
	assert.Equal(t, 1, counts["consequence"][CatDelinsNorm],
		"consequence mismatch with HGVSc delins_normalized should be reclassified")
	assert.Equal(t, 0, counts["consequence"][CatMismatch],
		"should have no consequence mismatches")
}

func TestCompareWriter_CrossColumnHGVSpDelinsNorm(t *testing.T) {
	// When HGVSc is delins_normalized, an HGVSp mismatch should be
	// reclassified as delins_normalized.
	var buf bytes.Buffer
	cols := map[string]bool{"hgvsp": true, "hgvsc": true}
	w := NewCompareWriter(&buf, cols, true)

	variant := &vcf.Variant{Chrom: "5", Pos: 140558745, Ref: "CGCAGGCGG", Alt: "GCC"}

	mafAnn := &maf.MAFAnnotation{
		HGVSpShort:   "p.R432_R434delinsA",
		HGVSc:        "ENST00000357560.9:c.1293_1301delinsGGC",
		TranscriptID: "ENST00000357560",
	}
	vepAnns := []*annotate.Annotation{{
		TranscriptID: "ENST00000357560",
		HGVSp:        "p.L435_K436del",
		HGVSc:        "c.1293_1300delinsGG",
	}}

	w.WriteComparison(variant, mafAnn, vepAnns)

	counts := w.Counts()
	assert.Equal(t, 1, counts["hgvsp"][CatDelinsNorm],
		"hgvsp mismatch with HGVSc delins_normalized should be reclassified")
	assert.Equal(t, 0, counts["hgvsp"][CatMismatch],
		"should have no hgvsp mismatches")
}
