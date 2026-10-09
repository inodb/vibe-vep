package annotate

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

func threeExonTranscript() *cache.Transcript {
	t := &cache.Transcript{
		ID: "ENST_SR", GeneName: "SR", Chrom: "1", Strand: 1,
		Start: 1000, End: 1311, CDSStart: 1000, CDSEnd: 1311,
		Biotype:     "protein_coding",
		CDSSequence: "ATGCATGAATAAATGCATGAATAAATGCATGAATAA",
		Exons: []cache.Exon{
			{Number: 1, Start: 1000, End: 1011, CDSStart: 1000, CDSEnd: 1011, Frame: 0},
			{Number: 2, Start: 1100, End: 1111, CDSStart: 1100, CDSEnd: 1111, Frame: 0},
			{Number: 3, Start: 1300, End: 1311, CDSStart: 1300, CDSEnd: 1311, Frame: 0},
		},
	}
	t.BuildCDSIndex()
	return t
}

func TestSpliceRegionNeedsAnAdjacentIntron(t *testing.T) {
	tr := threeExonTranscript()
	cases := []struct {
		name string
		pos  int64
		want bool
	}{
		// The transcript's outer edges have no intron beyond them.
		{"start of the first exon", 1001, false},
		{"end of the last exon", 1310, false},
		// Every interior boundary does.
		{"end of the first exon", 1010, true},
		{"start of the second exon", 1101, true},
		{"end of the second exon", 1110, true},
		{"start of the third exon", 1301, true},
		// 3-8 bases into an intron, on either side.
		{"donor side of the first intron", 1015, true},
		{"acceptor side of the first intron", 1095, true},
		// Deep inside an exon.
		{"middle of the second exon", 1106, false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if got := isSpliceRegion(c.pos, tr); got != c.want {
				t.Errorf("isSpliceRegion(%d) = %v, want %v", c.pos, got, c.want)
			}
		})
	}
}

func TestSpliceRegionSingleExonTranscript(t *testing.T) {
	// A gene with one exon has no splice sites anywhere, like the histone
	// H3C4 whose c.409T>C vibe-vep used to call a splice region.
	tr := &cache.Transcript{
		ID: "ENST_ONE", GeneName: "ONE", Chrom: "1", Strand: 1,
		Start: 1000, End: 1011, CDSStart: 1000, CDSEnd: 1011,
		Biotype:     "protein_coding",
		CDSSequence: "ATGCATGAATAA",
		Exons: []cache.Exon{
			{Number: 1, Start: 1000, End: 1011, CDSStart: 1000, CDSEnd: 1011, Frame: 0},
		},
	}
	tr.BuildCDSIndex()
	for _, pos := range []int64{1000, 1001, 1010, 1011} {
		if isSpliceRegion(pos, tr) {
			t.Errorf("isSpliceRegion(%d) = true for a single-exon transcript", pos)
		}
	}
}

func TestSpliceRegionUsesTheWholeSpan(t *testing.T) {
	tr := threeExonTranscript()
	// A deletion starting in the middle of exon 2 and running to its last
	// base: the start is nowhere near a boundary, the end is on one.
	v := &vcf.Variant{Chrom: "1", Pos: 1104, Ref: "AAAAAAA", Alt: ""}
	if !spanTouchesSpliceRegion(v, tr) {
		t.Error("a deletion reaching the exon boundary should touch the splice region")
	}
	// The same deletion kept clear of both boundaries.
	v2 := &vcf.Variant{Chrom: "1", Pos: 1104, Ref: "AAA", Alt: ""}
	if spanTouchesSpliceRegion(v2, tr) {
		t.Error("a deletion inside the exon should not touch the splice region")
	}
}
