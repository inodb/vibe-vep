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

// A deletion starting inside the acceptor flank and ending in the exon shifts
// 3' when its first deleted base repeats just past its end. The flank before
// exon 2 is TTTTTTTCAG and the exon opens ATGCATGAATAA, so deleting the two
// bases at -4 and -3 (T, T) can move along the T run.
func TestGenomicDeletionShiftAcrossAcceptor(t *testing.T) {
	tr := intronFlankTranscript(1)
	// Flank positions: 1090..1099 hold TTTTTTTCAG, so 1093 and 1094 are both T
	// and the base past the end (1095) is T as well.
	if got := genomicDeletionShift(1093, 1094, tr); got != 2 {
		t.Errorf("shift = %d, want 2 (the T run ends at 1096)", got)
	}
}

func TestGenomicDeletionShiftStopsWithoutSequence(t *testing.T) {
	tr := intronFlankTranscript(1)
	// 1050 is in the middle of the intron, past both flanks, so there is no
	// sequence to compare and the deletion cannot be shifted.
	if got := genomicDeletionShift(1050, 1052, tr); got != 0 {
		t.Errorf("shift = %d, want 0 beyond the cached flank", got)
	}
}

func TestGenomicDeletionShiftReverseStrand(t *testing.T) {
	// On the reverse strand 3' runs toward lower genomic positions. Exon 2
	// sits at 1000..1011 and its acceptor flank is above it, so genomic 1013
	// and 1012 are the A and G at -2 and -1, and the first base of the exon,
	// at genomic 1011, is the A that CDS position 13 holds.
	//
	// Deleting that AG therefore shifts one base: the deletion's first base
	// (A) repeats just past its end (A). It stops there, because the next
	// comparison is the G at 1012 against the T at CDS 14.
	tr := intronFlankTranscript(-1)
	if got := genomicDeletionShift(1012, 1013, tr); got != 1 {
		t.Errorf("shift = %d, want 1", got)
	}
}
