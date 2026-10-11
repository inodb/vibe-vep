package annotate

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// A reverse-strand transcript, where exon 1 is the rightmost genomic exon.
// This is the orientation that made the old start-position lookup fail: the
// leftmost genomic base of a deletion written c.N_N+Kdel is the intronic one.
func revStrandTranscript() *cache.Transcript {
	return &cache.Transcript{
		ID: "ENST_EXONNUM", GeneName: "EXONNUM", Chrom: "1", Strand: -1,
		Start: 100, End: 600, Biotype: "protein_coding",
		Exons: []cache.Exon{
			{Number: 3, Start: 100, End: 200},
			{Number: 2, Start: 300, End: 400},
			{Number: 1, Start: 500, End: 600},
		},
	}
}

func TestSpannedExonNumber(t *testing.T) {
	tr := revStrandTranscript()
	tests := []struct {
		name string
		pos  int64
		ref  string
		want string
	}{
		{"inside an exon", 350, "A", "2/3"},
		{"deletion running off the exon into the intron", 395, "AAAAAAAAAA", "2/3"},
		{"deletion starting in the intron and reaching the exon", 295, "AAAAAAAAAA", "2/3"},
		{"span covering two exons and the intron between", 150, string(make([]byte, 201)), "2-3/3"},
		{"span covering every exon", 100, string(make([]byte, 501)), "1-3/3"},
		{"wholly intronic", 220, "AAAA", ""},
		{"outside the transcript", 50, "A", ""},
		// A pure insertion covers no genomic base; it is in an exon only when
		// both of the bases it sits between are, so neither side of a
		// junction counts.
		{"junction insertion, intronic side", 201, "", ""},
		{"junction insertion, exonic side", 200, "", ""},
		{"junction insertion before an exon", 299, "", ""},
		{"insertion inside an exon", 350, "", "2/3"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := &vcf.Variant{Chrom: "1", Pos: tt.pos, Ref: tt.ref, Alt: ""}
			if got := spannedExonNumber(v, tr); got != tt.want {
				t.Errorf("pos %d ref len %d -> %q, want %q", tt.pos, len(tt.ref), got, tt.want)
			}
		})
	}
}

func TestExonNumberWhenTheVariantStartsOutsideTheTranscript(t *testing.T) {
	tr := revStrandTranscript() // exons at 100-200, 300-400, 500-600
	// A deletion beginning before the transcript and reaching exon 3.
	v := &vcf.Variant{Chrom: "1", Pos: 50, Ref: string(make([]byte, 101)), Alt: ""}
	if got := PredictConsequence(v, tr).ExonNumber; got != "3/3" {
		t.Errorf("ExonNumber = %q, want %q", got, "3/3")
	}
	// One that genuinely stays outside names no exon.
	v2 := &vcf.Variant{Chrom: "1", Pos: 50, Ref: "AAAA", Alt: ""}
	if got := PredictConsequence(v2, tr).ExonNumber; got != "" {
		t.Errorf("ExonNumber = %q, want empty", got)
	}
}
