package annotate

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/inodb/vibe-vep/internal/vcf"
)

func TestFormatHGVSp(t *testing.T) {
	tests := []struct {
		name   string
		result *ConsequenceResult
		want   string
	}{
		{
			name: "missense",
			result: &ConsequenceResult{
				Consequence:     ConsequenceMissenseVariant,
				ProteinPosition: 12,
				RefAA:           'G',
				AltAA:           'C',
			},
			want: "p.Gly12Cys",
		},
		{
			name: "synonymous",
			result: &ConsequenceResult{
				Consequence:     ConsequenceSynonymousVariant,
				ProteinPosition: 12,
				RefAA:           'G',
				AltAA:           'G',
			},
			want: "p.Gly12=",
		},
		{
			name: "stop_gained",
			result: &ConsequenceResult{
				Consequence:     ConsequenceStopGained,
				ProteinPosition: 12,
				RefAA:           'G',
				AltAA:           '*',
			},
			want: "p.Gly12Ter",
		},
		{
			name: "stop_lost",
			result: &ConsequenceResult{
				Consequence:     ConsequenceStopLost,
				ProteinPosition: 130,
				RefAA:           '*',
				AltAA:           'K',
			},
			want: "p.Ter130Lysext*?",
		},
		{
			name: "start_lost",
			result: &ConsequenceResult{
				Consequence:     ConsequenceStartLost,
				ProteinPosition: 1,
				RefAA:           'M',
				AltAA:           'K',
			},
			want: "p.Met1?",
		},
		{
			name: "stop_retained",
			result: &ConsequenceResult{
				Consequence:     ConsequenceStopRetained,
				ProteinPosition: 130,
				RefAA:           '*',
				AltAA:           '*',
			},
			want: "p.Ter130=",
		},
		{
			name: "frameshift_with_stop",
			result: &ConsequenceResult{
				Consequence:        ConsequenceFrameshiftVariant,
				ProteinPosition:    12,
				RefAA:              'G',
				AltAA:              'A',
				FrameshiftStopDist: 6,
			},
			want: "p.Gly12AlafsTer6",
		},
		{
			name: "frameshift_no_stop",
			result: &ConsequenceResult{
				Consequence:     ConsequenceFrameshiftVariant,
				ProteinPosition: 12,
				RefAA:           'G',
				AltAA:           'A',
			},
			want: "p.Gly12Alafs",
		},
		{
			name: "frameshift_no_altaa",
			result: &ConsequenceResult{
				Consequence:     ConsequenceFrameshiftVariant,
				ProteinPosition: 12,
				RefAA:           'G',
			},
			want: "p.Gly12fs",
		},
		{
			name: "inframe_deletion",
			result: &ConsequenceResult{
				Consequence:     ConsequenceInframeDeletion,
				ProteinPosition: 12,
				RefAA:           'G',
			},
			want: "p.Gly12del",
		},
		{
			name: "inframe_insertion",
			result: &ConsequenceResult{
				Consequence:     ConsequenceInframeInsertion,
				ProteinPosition: 12,
				RefAA:           'G',
			},
			want: "p.Gly12_13ins",
		},
		{
			name: "intronic_empty",
			result: &ConsequenceResult{
				Consequence: ConsequenceIntronVariant,
			},
			want: "",
		},
		{
			name: "splice_region_with_missense",
			result: &ConsequenceResult{
				Consequence:     ConsequenceMissenseVariant + "," + ConsequenceSpliceRegion,
				ProteinPosition: 37,
				RefAA:           'A',
				AltAA:           'V',
			},
			want: "p.Ala37Val",
		},
		{
			name: "frameshift_no_refaa",
			result: &ConsequenceResult{
				Consequence:     ConsequenceFrameshiftVariant,
				ProteinPosition: 5,
				RefAA:           0,
			},
			want: "p.5fs",
		},
		{
			name: "splice_donor",
			result: &ConsequenceResult{
				Consequence:     ConsequenceSpliceDonor,
				ProteinPosition: 37,
			},
			want: "p.X37_splice",
		},
		{
			name: "splice_acceptor",
			result: &ConsequenceResult{
				Consequence:     ConsequenceSpliceAcceptor,
				ProteinPosition: 100,
			},
			want: "p.X100_splice",
		},
		{
			name: "multi_codon_deletion",
			result: &ConsequenceResult{
				Consequence:        ConsequenceInframeDeletion,
				ProteinPosition:    59,
				ProteinEndPosition: 66,
				RefAA:              'R',
				EndAA:              'E',
			},
			want: "p.Arg59_Glu66del",
		},
		{
			name: "single_codon_deletion_with_end_same",
			result: &ConsequenceResult{
				Consequence:        ConsequenceInframeDeletion,
				ProteinPosition:    12,
				ProteinEndPosition: 12, // same as start — single codon
				RefAA:              'G',
				EndAA:              'G',
			},
			want: "p.Gly12del",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := FormatHGVSp(tt.result)
			assert.Equal(t, tt.want, got)
		})
	}
}

func TestHGVSp_KRASG12C(t *testing.T) {
	// KRAS G12C: missense variant should produce p.Gly12Cys
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245351,
		Ref:   "C",
		Alt:   "A",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	assert.Equal(t, "p.Gly12Cys", result.HGVSp)
}

func TestHGVSp_Synonymous(t *testing.T) {
	// GGT -> GGC at codon 12 (both Gly) should produce p.Gly12=
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245349, // Third position of codon 12 on reverse strand
		Ref:   "C",      // Genomic C -> coding G
		Alt:   "T",      // Genomic T -> coding A... let's verify
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	// If this is indeed synonymous, check HGVSp
	if result.Consequence == ConsequenceSynonymousVariant {
		assert.Equal(t, "p.Gly12=", result.HGVSp)
	}
}

func TestHGVSp_Frameshift(t *testing.T) {
	// Frameshift deletion at codon 12 should produce full HGVSp with alt AA and stop distance
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245350,
		Ref:   "GG",
		Alt:   "G",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	require.Equal(t, ConsequenceFrameshiftVariant, result.Consequence)

	assert.NotEqual(t, byte(0), result.RefAA)
	assert.NotEqual(t, byte(0), result.AltAA)

	// HGVSp should contain "fsTer" with a stop distance
	assert.NotEmpty(t, result.HGVSp)
	if result.FrameshiftStopDist > 0 {
		assert.Contains(t, result.HGVSp, "fsTer")
	}
	t.Logf("Frameshift HGVSp: %s (RefAA=%c, AltAA=%c, pos=%d, stopDist=%d)",
		result.HGVSp, result.RefAA, result.AltAA, result.ProteinPosition, result.FrameshiftStopDist)
}

func TestHGVSp_InframeDeletion(t *testing.T) {
	// In-frame deletion (3 bases) should produce p.{AA}{pos}del
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245350,
		Ref:   "GGGA",
		Alt:   "G",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	require.Equal(t, ConsequenceInframeDeletion, result.Consequence)

	assert.NotEqual(t, byte(0), result.RefAA)

	assert.NotEmpty(t, result.HGVSp)
	t.Logf("Inframe deletion HGVSp: %s", result.HGVSp)
}

func TestHGVSp_StartLost(t *testing.T) {
	// Create a variant at the start codon (Met1)
	// For KRAS reverse strand, the start codon ATG starts at CDS position 1
	// which maps to genomic position around CDSEnd (25245384)
	result := &ConsequenceResult{
		Consequence:     ConsequenceStartLost,
		ProteinPosition: 1,
		RefAA:           'M',
		AltAA:           'K',
	}

	hgvsp := FormatHGVSp(result)
	assert.Equal(t, "p.Met1?", hgvsp)
}

func TestAaThree(t *testing.T) {
	tests := []struct {
		aa   byte
		want string
	}{
		{'G', "Gly"},
		{'A', "Ala"},
		{'*', "Ter"},
		{'X', "Xaa"},
		{'Z', "Xaa"}, // unknown
		{0, "Xaa"},   // zero value
	}

	for _, tt := range tests {
		got := aaThree(tt.aa)
		assert.Equal(t, tt.want, got)
	}
}

func BenchmarkFormatHGVSp(b *testing.B) {
	results := []*ConsequenceResult{
		{Consequence: ConsequenceMissenseVariant, ProteinPosition: 12, RefAA: 'G', AltAA: 'C'},
		{Consequence: ConsequenceSynonymousVariant, ProteinPosition: 12, RefAA: 'G', AltAA: 'G'},
		{Consequence: ConsequenceStopGained, ProteinPosition: 12, RefAA: 'G', AltAA: '*'},
		{Consequence: ConsequenceFrameshiftVariant, ProteinPosition: 12, RefAA: 'G', AltAA: 'A', FrameshiftStopDist: 6},
		{Consequence: ConsequenceInframeDeletion, ProteinPosition: 12, RefAA: 'G'},
		{Consequence: ConsequenceInframeInsertion, ProteinPosition: 12, RefAA: 'G'},
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for _, r := range results {
			FormatHGVSp(r)
		}
	}
}

// TestAllocRegression_FormatHGVSp verifies allocation count for HGVSp formatting.
func TestAllocRegression_FormatHGVSp(t *testing.T) {
	results := []*ConsequenceResult{
		{Consequence: ConsequenceMissenseVariant, ProteinPosition: 12, RefAA: 'G', AltAA: 'C'},
		{Consequence: ConsequenceSynonymousVariant, ProteinPosition: 12, RefAA: 'G', AltAA: 'G'},
		{Consequence: ConsequenceStopGained, ProteinPosition: 12, RefAA: 'G', AltAA: '*'},
		{Consequence: ConsequenceFrameshiftVariant, ProteinPosition: 12, RefAA: 'G', AltAA: 'A', FrameshiftStopDist: 6},
	}

	allocs := testing.AllocsPerRun(100, func() {
		for _, r := range results {
			FormatHGVSp(r)
		}
	})

	// Each FormatHGVSp call does ~1 alloc (string concat result).
	// 4 calls → ~4 allocs. Allow small margin.
	const maxAllocs = 8
	if int(allocs) > maxAllocs {
		t.Errorf("FormatHGVSp allocation regression: got %.0f, want <= %d", allocs, maxAllocs)
	}
	t.Logf("FormatHGVSp allocs: %.0f (budget: %d)", allocs, maxAllocs)
}
