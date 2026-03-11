package output

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/inodb/vibe-vep/internal/annotate"
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
		// Fix 2: non-coding boundary shifts
		{"non_coding_exon ↔ intron boundary shift", "non_coding_transcript_exon_variant", "intron_variant", CatTranscriptModelChange},
		{"intron ↔ non_coding_exon boundary shift", "intron_variant", "non_coding_transcript_exon_variant", CatTranscriptModelChange},
		{"intron ↔ 3'UTR exon boundary shift", "intron_variant", "3_prime_UTR_variant", CatUpstreamReclass},
		{"non_coding_exon ↔ 5'UTR", "non_coding_transcript_exon_variant", "5_prime_UTR_variant", CatTranscriptModelChange},
		{"intron ↔ intergenic non-coding gene model", "intron_variant", "intergenic_variant", CatGeneModelChange},
		// Fix 4: splice reclassification with inframe/start_lost
		{"inframe+splice → splice_acceptor", "inframe_insertion,splice_region_variant", "splice_acceptor_variant", CatMatch},
		{"protein_altering → splice_acceptor", "protein_altering_variant,splice_region_variant", "splice_acceptor_variant", CatMatch},
		{"start_lost → splice_acceptor", "splice_acceptor_variant,coding_sequence_variant,5_prime_UTR_variant,intron_variant", "start_lost", CatMatch},
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
		name                          string
		mafHGVSp, vepHGVSp, vepConseq string
		mafHGVSc, vepHGVSc            string
		want                          Category
	}{
		{
			name:     "exact match",
			mafHGVSp: "p.G12C", vepHGVSp: "p.G12C",
			want: CatMatch,
		},
		{
			name:     "both empty",
			mafHGVSp: "", vepHGVSp: "",
			want: CatBothEmpty,
		},
		{
			name:     "fuzzy frameshift",
			mafHGVSp: "p.G12fs", vepHGVSp: "p.G15fs*33",
			want: CatFuzzyFS,
		},
		{
			name:     "splice vs syn",
			mafHGVSp: "p.X125_splice", vepHGVSp: "p.G125=",
			want: CatSpliceVsSyn,
		},
		{
			name:     "maf nonstandard intronic",
			mafHGVSp: "p.*130*", vepHGVSp: "",
			want: CatMafNonstandard,
		},
		{
			name:     "splice no protein",
			mafHGVSp: "p.X125_splice", vepHGVSp: "",
			vepConseq: "splice_acceptor_variant",
			want:      CatSpliceNoProtein,
		},
		{
			name:     "maf empty",
			mafHGVSp: "", vepHGVSp: "p.G12C",
			want: CatMafEmpty,
		},
		{
			name:     "vep empty",
			mafHGVSp: "p.G12C", vepHGVSp: "",
			want: CatVepEmpty,
		},
		{
			name:     "position shift - same change type",
			mafHGVSp: "p.Y1145C", vepHGVSp: "p.Y1215C",
			want: CatPositionShift,
		},
		{
			name:     "position shift - synonymous",
			mafHGVSp: "p.A735=", vepHGVSp: "p.A744=",
			want: CatPositionShift,
		},
		{
			name:     "position shift - hgvsc matches",
			mafHGVSp: "p.G12C", vepHGVSp: "p.G15C",
			mafHGVSc: "ENST00000256078.4:c.34G>T", vepHGVSc: "c.34G>T",
			want: CatPositionShift,
		},
		{
			name:     "frameshift vs stop_gained",
			mafHGVSp: "p.E454Sfs*117", vepHGVSp: "p.K453*",
			want: CatFuzzyFS,
		},
		{
			name:     "stop_gained vs frameshift",
			mafHGVSp: "p.K453*", vepHGVSp: "p.E454Sfs*117",
			want: CatFuzzyFS,
		},
		{
			name:     "maf nonstandard intronic with vep prediction - transcript model change",
			mafHGVSp: "p.*130*", vepHGVSp: "p.G12C",
			want: CatTranscriptModelChange,
		},
		{
			name:     "mismatch",
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
		// Fix 1: n./c. transcript model change
		{"n. vs c. transcript model change", "ENST00000123456.1:n.100C>G", "c.50C>G", CatTranscriptModelChange},
		{"c. vs n. transcript model change", "ENST00000123456.1:c.100C>G", "n.50C>G", CatTranscriptModelChange},
		{"both n. different coords", "ENST00000123456.1:n.100C>G", "n.200G>C", CatTranscriptModelChange},
		// Fix 6: insertion cyclic rotation
		{"cyclic rotation ins", "ENST00000123456.1:c.100_101insTA", "c.101_102insAT", CatPositionShift},
		{"non-rotation ins stays mismatch", "ENST00000123456.1:c.100_101insAA", "c.101_102insGG", CatMismatch},
		{"mismatch different SNV", "ENST00000361923.2:c.1428C>G", "c.999A>T", CatMismatch},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			assert.Equal(t, tt.want, categorizeHGVSc(tt.maf, tt.vep))
		})
	}
}

func TestPickBestAnnotation(t *testing.T) {
	tests := []struct {
		name string
		anns []*annotate.Annotation
		want string // expected TranscriptID
	}{
		{
			name: "empty",
			anns: nil,
			want: "",
		},
		{
			name: "single",
			anns: []*annotate.Annotation{{TranscriptID: "T1"}},
			want: "T1",
		},
		{
			name: "prefer canonical",
			anns: []*annotate.Annotation{
				{TranscriptID: "T1", Impact: "HIGH", Biotype: "protein_coding"},
				{TranscriptID: "T2", Impact: "MODERATE", Biotype: "protein_coding", IsCanonicalMSK: true},
			},
			want: "T2",
		},
		{
			name: "prefer protein-coding over non-coding",
			anns: []*annotate.Annotation{
				{TranscriptID: "T1", Impact: "HIGH", Biotype: "retained_intron"},
				{TranscriptID: "T2", Impact: "MODERATE", Biotype: "protein_coding"},
			},
			want: "T2",
		},
		{
			name: "prefer higher impact when both coding",
			anns: []*annotate.Annotation{
				{TranscriptID: "T1", Impact: "MODERATE", Biotype: "protein_coding"},
				{TranscriptID: "T2", Impact: "HIGH", Biotype: "protein_coding"},
			},
			want: "T2",
		},
		{
			name: "prefer HGVSp when all else equal",
			anns: []*annotate.Annotation{
				{TranscriptID: "T1", Impact: "MODERATE", Biotype: "protein_coding"},
				{TranscriptID: "T2", Impact: "MODERATE", Biotype: "protein_coding", HGVSp: "p.G12C"},
			},
			want: "T2",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := PickBestAnnotation(tt.anns)
			if tt.want == "" {
				assert.Nil(t, result)
			} else {
				require.NotNil(t, result)
				assert.Equal(t, tt.want, result.TranscriptID)
			}
		})
	}
}

func TestPickMostSevere(t *testing.T) {
	tests := []struct {
		name string
		anns []*annotate.Annotation
		want string
	}{
		{
			name: "empty",
			anns: nil,
			want: "",
		},
		{
			name: "prefer highest impact",
			anns: []*annotate.Annotation{
				{TranscriptID: "T1", Impact: "MODERATE", Biotype: "protein_coding", IsCanonicalMSK: true},
				{TranscriptID: "T2", Impact: "HIGH", Biotype: "retained_intron"},
			},
			want: "T2",
		},
		{
			name: "tie-break canonical",
			anns: []*annotate.Annotation{
				{TranscriptID: "T1", Impact: "HIGH", Biotype: "protein_coding"},
				{TranscriptID: "T2", Impact: "HIGH", Biotype: "protein_coding", IsCanonicalMSK: true},
			},
			want: "T2",
		},
		{
			name: "tie-break protein-coding",
			anns: []*annotate.Annotation{
				{TranscriptID: "T1", Impact: "MODERATE", Biotype: "retained_intron"},
				{TranscriptID: "T2", Impact: "MODERATE", Biotype: "protein_coding"},
			},
			want: "T2",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := PickMostSevere(tt.anns)
			if tt.want == "" {
				assert.Nil(t, result)
			} else {
				require.NotNil(t, result)
				assert.Equal(t, tt.want, result.TranscriptID)
			}
		})
	}
}

func TestFormatAllEffects(t *testing.T) {
	tests := []struct {
		name string
		anns []*annotate.Annotation
		want string
	}{
		{
			name: "nil",
			anns: nil,
			want: "",
		},
		{
			name: "empty",
			anns: []*annotate.Annotation{},
			want: "",
		},
		{
			name: "single annotation",
			anns: []*annotate.Annotation{{
				GeneName:       "KRAS",
				Consequence:    "missense_variant",
				HGVSp:          "p.Gly12Cys",
				TranscriptID:   "ENST00000311936",
				HGVSc:          "c.34G>T",
				Impact:         "MODERATE",
				IsCanonicalMSK: true,
			}},
			want: "KRAS,missense_variant,p.G12C,ENST00000311936,c.34G>T,MODERATE,YES",
		},
		{
			name: "multiple annotations",
			anns: []*annotate.Annotation{
				{
					GeneName:       "KRAS",
					Consequence:    "missense_variant",
					HGVSp:          "p.Gly12Cys",
					TranscriptID:   "ENST00000311936",
					HGVSc:          "c.34G>T",
					Impact:         "MODERATE",
					IsCanonicalMSK: true,
				},
				{
					GeneName:     "KRAS",
					Consequence:  "intron_variant",
					TranscriptID: "ENST00000556131",
					Impact:       "MODIFIER",
				},
			},
			want: "KRAS,missense_variant,p.G12C,ENST00000311936,c.34G>T,MODERATE,YES;KRAS,intron_variant,,ENST00000556131,,MODIFIER,",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			assert.Equal(t, tt.want, FormatAllEffects(tt.anns))
		})
	}
}

func TestCategorizer_CategorizeRow(t *testing.T) {
	cat := &Categorizer{}

	t.Run("cross-column position_shift reclassifies consequence", func(t *testing.T) {
		columns := []string{"Consequence", "HGVSc"}
		left := map[string]string{
			"Consequence": "synonymous_variant",
			"HGVSc":       "ENST00000333418.4:c.390T>G",
		}
		right := map[string]string{
			"Consequence": "missense_variant",
			"HGVSc":       "c.388T>G",
		}
		cats := cat.CategorizeRow(columns, left, right,
			[]string{"Consequence", "HGVSc"}, []string{"Consequence", "HGVSc"})

		assert.Equal(t, CatPositionShift, cats["Consequence"])
		assert.Equal(t, CatPositionShift, cats["HGVSc"])
	})

	t.Run("cross-column delins_normalized reclassifies consequence", func(t *testing.T) {
		columns := []string{"Consequence", "HGVSc"}
		left := map[string]string{
			"Consequence": "Missense_Mutation",
			"HGVSc":       "ENST00000123456.1:c.100_102delinsGTA",
		}
		right := map[string]string{
			"Consequence": "synonymous_variant",
			"HGVSc":       "c.102G>A",
		}
		cats := cat.CategorizeRow(columns, left, right,
			[]string{"Consequence", "HGVSc"}, []string{"Consequence", "HGVSc"})

		assert.Equal(t, CatDelinsNorm, cats["Consequence"])
	})

	t.Run("cross-column dup_vs_ins reclassifies hgvsp", func(t *testing.T) {
		columns := []string{"HGVSp_Short", "HGVSc"}
		left := map[string]string{
			"HGVSp_Short": "p.V78dup",
			"HGVSc":       "ENST00000123456.1:c.232_233dup",
		}
		right := map[string]string{
			"HGVSp_Short": "p.V79_80insG",
			"HGVSc":       "c.233_234insGT",
		}
		cats := cat.CategorizeRow(columns, left, right,
			[]string{"HGVSp_Short", "HGVSc"}, []string{"HGVSp_Short", "HGVSc"})

		assert.Equal(t, CatDupVsIns, cats["HGVSp_Short"])
	})

	t.Run("simple string equality for unknown columns", func(t *testing.T) {
		columns := []string{"Hugo_Symbol"}
		left := map[string]string{"Hugo_Symbol": "KRAS"}
		right := map[string]string{"Hugo_Symbol": "BRAF"}
		cats := cat.CategorizeRow(columns, left, right,
			[]string{"Hugo_Symbol"}, []string{"Hugo_Symbol"})

		assert.Equal(t, CatMismatch, cats["Hugo_Symbol"])
	})

	t.Run("simple match for unknown columns", func(t *testing.T) {
		columns := []string{"Hugo_Symbol"}
		left := map[string]string{"Hugo_Symbol": "KRAS"}
		right := map[string]string{"Hugo_Symbol": "KRAS"}
		cats := cat.CategorizeRow(columns, left, right,
			[]string{"Hugo_Symbol"}, []string{"Hugo_Symbol"})

		assert.Equal(t, CatMatch, cats["Hugo_Symbol"])
	})
}

func TestValidateExcludeColumns(t *testing.T) {
	t.Run("valid columns", func(t *testing.T) {
		err := ValidateExcludeColumns([]string{"hugo_symbol", "all_effects", "canonical_ensembl"})
		assert.NoError(t, err)
	})

	t.Run("all valid core columns", func(t *testing.T) {
		err := ValidateExcludeColumns([]string{
			"hugo_symbol", "consequence", "variant_classification", "transcript_id",
			"hgvsc", "hgvsp", "hgvsp_short", "canonical_mskcc", "canonical_ensembl",
			"all_effects",
		})
		assert.NoError(t, err)
	})

	t.Run("invalid column", func(t *testing.T) {
		err := ValidateExcludeColumns([]string{"hugo_symbol", "not_a_column"})
		require.Error(t, err)
		assert.Contains(t, err.Error(), "not_a_column")
		assert.Contains(t, err.Error(), "valid columns")
	})

	t.Run("empty list", func(t *testing.T) {
		err := ValidateExcludeColumns(nil)
		assert.NoError(t, err)
	})
}

func TestValidOutputColumns(t *testing.T) {
	cols := ValidOutputColumns()
	// 9 core + all_effects = 10
	assert.Len(t, cols, 10)
	assert.True(t, cols["hugo_symbol"])
	assert.True(t, cols["all_effects"])
	assert.False(t, cols["not_a_column"])
}
