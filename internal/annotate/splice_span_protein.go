package annotate

import (
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// spliceSpanProteinTerms names the protein consequence VEP adds to a deletion
// that destroys a splice site.
//
// VEP decides whether to describe the protein at all by where the deletion
// ends: with both ends inside exons it has a CDS coordinate for each and names
// a length term, and with either end in an intron it names none. In a VEP111
// MSK-IMPACT MAF that split is exact -- 66 of 66 deletions hitting a splice
// site with both ends exonic carry a protein term, and 3467 of 3467 with an
// intronic end carry none.
//
// The terms themselves follow from the deleted CDS length and from what VEP
// takes as the reference codons. Those are the GENOMIC bases the deletion
// removes, introns and all: VEP's Codons column for CTNNB1 c.26_403del is the
// whole 578-base reference allele over "-", not the 378 spliced ones. A
// sequence that long almost always holds an in-frame stop, the empty
// alternative holds none, and VEP's stop_lost test -- reference peptide has a
// stop, alternative does not -- therefore fires: 61 of the 62 such rows carry
// stop_lost.
//
// Returns nil when the deletion is not of this shape.
func spliceSpanProteinTerms(v *vcf.Variant, t *cache.Transcript) []string {
	if len(v.Ref) < 2 || len(v.Alt) != 0 || !t.IsProteinCoding() {
		return nil
	}
	start, end := v.Pos, v.Pos+int64(len(v.Ref))-1
	if t.FindExon(start) == nil || t.FindExon(end) == nil {
		return nil
	}
	// A deletion reaching into a UTR is described by the boundary it crosses,
	// not by a length term, and those paths handle it already.
	if GenomicToCDS(start, t) < 1 || GenomicToCDS(end, t) < 1 {
		return nil
	}

	deleted := deletedCDSBases(start, end, t)
	if deleted == 0 {
		return nil
	}
	terms := []string{ConsequenceInframeDeletion}
	if deleted%3 != 0 {
		terms = []string{ConsequenceFrameshiftVariant}
	}
	if deletedSpanHasStop(v, t) {
		terms = append(terms, ConsequenceStopLost)
	}
	return terms
}

// deletedCDSBases counts the coding bases a genomic span removes, skipping the
// introns it crosses.
func deletedCDSBases(start, end int64, t *cache.Transcript) int64 {
	var n int64
	for i := range t.Exons {
		e := &t.Exons[i]
		if e.CDSStart == 0 || e.CDSEnd == 0 {
			continue
		}
		lo, hi := max64(start, e.CDSStart), min64(end, e.CDSEnd)
		if hi >= lo {
			n += hi - lo + 1
		}
	}
	return n
}

// deletedSpanHasStop reports whether the removed genomic sequence, read in the
// transcript's orientation and in the frame the first deleted base sits in,
// contains a stop codon. That sequence is what VEP translates as the reference
// peptide for these deletions (see spliceSpanProteinTerms).
func deletedSpanHasStop(v *vcf.Variant, t *cache.Transcript) bool {
	seq := v.Ref
	first := v.Pos
	if !t.IsForwardStrand() {
		seq = ReverseComplement(seq)
		first = v.Pos + int64(len(v.Ref)) - 1
	}
	offset := 0
	if cds := GenomicToCDS(first, t); cds >= 1 {
		offset = int((cds - 1) % 3)
	}
	// The partial codon the deletion starts inside is completed by bases that
	// survive, so the stop search begins at the first whole codon it removes.
	for i := (3 - offset) % 3; i+3 <= len(seq); i += 3 {
		if TranslateCodon(seq[i:i+3]) == '*' {
			return true
		}
	}
	return false
}

func max64(a, b int64) int64 {
	if a > b {
		return a
	}
	return b
}

func min64(a, b int64) int64 {
	if a < b {
		return a
	}
	return b
}
