package annotate

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Two coding exons with an intron between them. Exon 1 holds CDS bases 1..12
// at genomic 1000..1011, so the intron's donor side starts at 1012: +1 is the
// G and +2 the T of the GT dinucleotide.
func intronFlankTranscript(strand int8) *cache.Transcript {
	after, afterLen := cache.PackDNA2Bit("GTAAGCCCCC")
	before, beforeLen := cache.PackDNA2Bit("TTTTTTTCAG")
	e1 := cache.Exon{Number: 1, Start: 1000, End: 1011, CDSStart: 1000, CDSEnd: 1011, Frame: 0,
		IntronAfterPacked: after, IntronAfterLen: uint8(afterLen)}
	e2 := cache.Exon{Number: 2, Start: 1100, End: 1111, CDSStart: 1100, CDSEnd: 1111, Frame: 0,
		IntronBeforePacked: before, IntronBeforeLen: uint8(beforeLen)}
	if strand < 0 {
		// On the reverse strand the first exon is the rightmost one, and each
		// exon's flanks keep their coding-orientation meaning.
		e1 = cache.Exon{Number: 1, Start: 1100, End: 1111, CDSStart: 1100, CDSEnd: 1111, Frame: 0,
			IntronAfterPacked: after, IntronAfterLen: uint8(afterLen)}
		e2 = cache.Exon{Number: 2, Start: 1000, End: 1011, CDSStart: 1000, CDSEnd: 1011, Frame: 0,
			IntronBeforePacked: before, IntronBeforeLen: uint8(beforeLen)}
	}
	t := &cache.Transcript{
		ID: "ENST_FLANK", GeneName: "FLANK", Chrom: "1", Strand: strand,
		Start: 1000, End: 1111, CDSStart: 1000, CDSEnd: 1111,
		Biotype:     "protein_coding",
		CDSSequence: "ATGCATGAATAAATGCATGAATAA",
		Exons:       []cache.Exon{e1, e2},
	}
	if strand < 0 {
		t.Exons = []cache.Exon{e2, e1} // exons stay in genomic order
	}
	t.BuildCDSIndex()
	return t
}

func TestIntronicInsertionDup(t *testing.T) {
	tr := intronFlankTranscript(1)
	// Inserting a T between +1 and +2 shifts one base 3', onto the T already
	// there, and is then a duplication of it.
	v := &vcf.Variant{Chrom: "1", Pos: 1012, Ref: "", Alt: "T"}
	result := PredictConsequence(v, tr)
	got := FormatHGVSc(v, tr, result)
	if want := "c.12+2dup"; got != want {
		t.Errorf("HGVSc = %q, want %q", got, want)
	}
	if result.HGVSOffset != 1 {
		t.Errorf("HGVS_Offset = %d, want 1", result.HGVSOffset)
	}
}

func TestIntronicInsertionNotADup(t *testing.T) {
	tr := intronFlankTranscript(1)
	// A C duplicates neither neighbour, so it stays a plain insertion.
	v := &vcf.Variant{Chrom: "1", Pos: 1012, Ref: "", Alt: "C"}
	result := PredictConsequence(v, tr)
	if got, want := FormatHGVSc(v, tr, result), "c.12+1_12+2insC"; got != want {
		t.Errorf("HGVSc = %q, want %q", got, want)
	}
	if result.HGVSOffset != 0 {
		t.Errorf("HGVS_Offset = %d, want 0", result.HGVSOffset)
	}
}

func TestIntronFlankBaseAtOrientation(t *testing.T) {
	for _, strand := range []int8{1, -1} {
		tr := intronFlankTranscript(strand)
		// +1 and +2 of the donor are G and T whichever way the gene runs.
		donor1, donor2 := int64(1012), int64(1013)
		if strand < 0 {
			donor1, donor2 = 1099, 1098
		}
		for pos, want := range map[int64]byte{donor1: 'G', donor2: 'T'} {
			got, ok := intronFlankBaseAt(pos, tr)
			if !ok || got != want {
				t.Errorf("strand %d: base at %d = %q (ok=%v), want %q", strand, pos, got, ok, want)
			}
		}
		// An exonic position is not in any flank.
		if _, ok := intronFlankBaseAt(1005, tr); ok {
			t.Errorf("strand %d: exonic position reported as intron flank", strand)
		}
	}
}
