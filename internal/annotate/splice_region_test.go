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

	// An insertion needs both of the bases it sits between to be in the
	// region. Exon 2 runs 1100..1111, so its first three bases are the
	// region: between 1101 and 1102 is inside it, between 1102 and 1103 is
	// not.
	ins := &vcf.Variant{Chrom: "1", Pos: 1101, Ref: "", Alt: "T"}
	if !spanTouchesSpliceRegion(ins, tr) {
		t.Error("an insertion between the second and third exon bases is in the splice region")
	}
	insOut := &vcf.Variant{Chrom: "1", Pos: 1102, Ref: "", Alt: "T"}
	if spanTouchesSpliceRegion(insOut, tr) {
		t.Error("an insertion between the third and fourth exon bases is past the splice region")
	}

	// An insertion on the junction itself: exon 1 ends at 1010, so this sits
	// between its last base and the intron's +1. Neither window holds both
	// bases, but VEP reports the splice region.
	junction := &vcf.Variant{Chrom: "1", Pos: 1011, Ref: "", Alt: "T"}
	if !spanTouchesSpliceRegion(junction, tr) {
		t.Error("an insertion on an exon/intron junction is in the splice region")
	}
	// The same shape at the transcript's outer edge borders no intron.
	outer := &vcf.Variant{Chrom: "1", Pos: 999, Ref: "", Alt: "T"}
	if spanTouchesSpliceRegion(outer, tr) {
		t.Error("the start of the first exon borders no intron")
	}
}

// Exon 1 ends at 1010 (CDS 12 is 1011), so the donor's +1 and +2 are 1012 and
// 1013; the acceptor's -1 and -2 before exon 2 are 1099 and 1098.
func TestInsertionSpliceTerms(t *testing.T) {
	tr := threeExonTranscript()
	cases := []struct {
		name string
		pos  int64
		want string
	}{
		// Between +1 and +2: the insertion splits the donor GT.
		{"splits the donor", 1012, ConsequenceSpliceDonor},
		// Between +2 and +3: the GT survives, so this is only a region.
		{"past the donor", 1013, ConsequenceSpliceRegion + "," + ConsequenceIntronVariant},
		// Between -2 and -1: splits the acceptor AG.
		{"splits the acceptor", 1098, ConsequenceSpliceAcceptor},
		// Between -3 and -2: the AG survives.
		{"before the acceptor", 1097, ConsequenceSpliceRegion + "," +
			ConsequencePolypyrimidineTract + "," + ConsequenceIntronVariant},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			v := &vcf.Variant{Chrom: "1", Pos: c.pos, Ref: "", Alt: "T"}
			terms, ok := insertionSpliceTerms(v, tr)
			if !ok {
				t.Fatalf("expected both flanks to be intronic at %d", c.pos)
			}
			if got := joinConsequenceTerms(terms); got != c.want {
				t.Errorf("terms = %q, want %q", got, c.want)
			}
		})
	}
	// An insertion with an exonic flank is not this function's business.
	v := &vcf.Variant{Chrom: "1", Pos: 1011, Ref: "", Alt: "T"}
	if _, ok := insertionSpliceTerms(v, tr); ok {
		t.Error("an insertion at the exon edge should be left to the exonic path")
	}
}

// Exon 2 runs 1100..1111 with the CDS starting at 1106, so 1101 is a 5'UTR
// base two inside the exon - in the splice region, and in the UTR.
func TestUTRVariantInSpliceRegion(t *testing.T) {
	tr := threeExonTranscript()
	tr.CDSStart = 1106
	tr.BuildCDSIndex()
	v := &vcf.Variant{Chrom: "1", Pos: 1101, Ref: "G", Alt: "C"}
	result := PredictConsequence(v, tr)
	want := ConsequenceSpliceRegion + "," + Consequence5PrimeUTR
	if result.Consequence != want {
		t.Errorf("Consequence = %q, want %q", result.Consequence, want)
	}
	// Deeper into the exon it is only a UTR variant.
	v2 := &vcf.Variant{Chrom: "1", Pos: 1104, Ref: "G", Alt: "C"}
	if got := PredictConsequence(v2, tr).Consequence; got != Consequence5PrimeUTR {
		t.Errorf("Consequence = %q, want %q", got, Consequence5PrimeUTR)
	}
}

// Ensembl orders consequence terms by severity, so splice_region_variant lands
// before the terms it outranks and after those that outrank it.
func TestAppendSpliceRegionOrdersBySeverity(t *testing.T) {
	tests := map[string]string{
		ConsequenceMissenseVariant:       ConsequenceMissenseVariant + "," + ConsequenceSpliceRegion,
		ConsequenceFrameshiftVariant:     ConsequenceFrameshiftVariant + "," + ConsequenceSpliceRegion,
		ConsequenceSynonymousVariant:     ConsequenceSpliceRegion + "," + ConsequenceSynonymousVariant,
		ConsequenceCodingSequenceVariant: ConsequenceSpliceRegion + "," + ConsequenceCodingSequenceVariant,
		// Not in the fast-path table: the term goes in the middle.
		ConsequenceStopLost + "," + Consequence3PrimeUTR: ConsequenceStopLost + "," +
			ConsequenceSpliceRegion + "," + Consequence3PrimeUTR,
	}
	for in, want := range tests {
		if got := appendSpliceRegion(in); got != want {
			t.Errorf("appendSpliceRegion(%q) = %q, want %q", in, got, want)
		}
	}
}
