package annotate

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Real variants from an MSK-IMPACT MAF, all typed INS but carrying a
// reference allele that the alternate replaces rather than extends.
func TestLongerAltThatReplacesTheReference(t *testing.T) {
	tr := &cache.Transcript{
		ID: "ENST_DELINS", GeneName: "DELINS", Chrom: "1", Strand: 1,
		Start: 1000, End: 1023, CDSStart: 1000, CDSEnd: 1023,
		Biotype:     "protein_coding",
		CDSSequence: "ATGCACCATGAATAAATGCATGAA",
		Exons: []cache.Exon{
			{Number: 1, Start: 1000, End: 1023, CDSStart: 1000, CDSEnd: 1023, Frame: 0},
		},
	}
	tr.BuildCDSIndex()
	tests := []struct {
		name     string
		pos      int64
		ref, alt string
		want     string
	}{
		// Nothing shared: the whole span is replaced. Shaped like NF1 ref C
		// alt GG, which used to come out as an insertion of a single G.
		{"one base replaced by two", 1003, "C", "GG", "c.4delinsGG"},
		{"four bases replaced by five", 1003, "CACC", "TTTCG", "c.4_7delinsTTTCG"},
		// The reference survives as a prefix, so this really is an insertion.
		{"anchored insertion stays an insertion", 1003, "C", "CGG", "c.4_5insGG"},
		// And so is the stripped form the pipeline sends.
		{"bare insertion stays an insertion", 1003, "", "GG", "c.4_5insGG"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := &vcf.Variant{Chrom: "1", Pos: tt.pos, Ref: tt.ref, Alt: tt.alt}
			result := PredictConsequence(v, tr)
			if got := FormatHGVSc(v, tr, result); got != tt.want {
				t.Errorf("ref %q alt %q -> %q, want %q", tt.ref, tt.alt, got, tt.want)
			}
		})
	}
}
