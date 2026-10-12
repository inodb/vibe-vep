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

// A deletion that begins before the transcript but deletes part of its first
// exon is not upstream of it.
func TestDeletionReachingIntoTheTranscriptIsNotUpstream(t *testing.T) {
	tr := &cache.Transcript{
		ID: "ENST_REACH", GeneName: "REACH", Chrom: "1", Strand: 1,
		Start: 1000, End: 1023, CDSStart: 1006, CDSEnd: 1017,
		Biotype:     "protein_coding",
		CDSSequence: "ATGCATGAATAA",
		Exons: []cache.Exon{
			{Number: 1, Start: 1000, End: 1023, CDSStart: 1006, CDSEnd: 1017, Frame: 0},
		},
	}
	tr.BuildCDSIndex()
	// Starts 10 bases before the transcript and runs to CDS base 3, so it
	// covers upstream sequence, the 5'UTR and coding sequence.
	v := &vcf.Variant{Chrom: "1", Pos: 990, Ref: string(make([]byte, 19)), Alt: ""}
	result := PredictConsequence(v, tr)
	if result.Consequence == ConsequenceUpstreamGene {
		t.Fatalf("called upstream_gene_variant while naming exon %q", result.ExonNumber)
	}
	want := ConsequenceCodingSequenceVariant + "," + Consequence5PrimeUTR
	if result.Consequence != want {
		t.Errorf("Consequence = %q, want %q", result.Consequence, want)
	}
	if result.ExonNumber != "1/1" {
		t.Errorf("ExonNumber = %q, want %q", result.ExonNumber, "1/1")
	}
	// Its first base is outside the transcript, so there is no c. coordinate
	// to name it from and VEP writes no HGVSc either.
	if got := FormatHGVSc(v, tr, result); got != "" {
		t.Errorf("HGVSc = %q, want empty", got)
	}
	// One that really stays outside is still upstream.
	v2 := &vcf.Variant{Chrom: "1", Pos: 980, Ref: "AAAA", Alt: ""}
	if got := PredictConsequence(v2, tr).Consequence; got != ConsequenceUpstreamGene {
		t.Errorf("Consequence = %q, want %q", got, ConsequenceUpstreamGene)
	}
}
