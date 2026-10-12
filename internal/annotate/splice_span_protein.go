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
	if !t.IsProteinCoding() {
		return nil
	}
	seq, start := deletedBases(v)
	if len(seq) < 2 {
		return nil
	}
	end := start + int64(len(seq)) - 1
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
	if deletedSpanHasStop(seq, start, t) {
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
func deletedSpanHasStop(seq string, start int64, t *cache.Transcript) bool {
	first := start
	if !t.IsForwardStrand() {
		seq = ReverseComplement(seq)
		first = start + int64(len(seq)) - 1
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

// genomicSpanCodons renders the Codons string for a deletion that crosses an
// intron.
//
// VEP's Codons column reports the bases the deletion removes from the GENOME,
// introns included, padded out to codon boundaries with the surviving flanks:
// for CTNNB1 c.26_403del it is the whole 578-base reference allele over "-",
// not the 378 spliced ones. Building the window from the CDS alone, as every
// other deletion is built, left out the intronic bases and produced a shorter
// string.
//
// Returns "" for a deletion that stays inside one exon, where the two are the
// same and the ordinary CDS window is correct.
func genomicSpanCodons(v *vcf.Variant, t *cache.Transcript, delStartCDS, delEndCDS int64) string {
	n := int64(len(t.CDSSequence))
	if delStartCDS < 1 || delEndCDS < delStartCDS || delEndCDS > n {
		return ""
	}
	seq, start := deletedBases(v)
	if len(seq) == 0 || withinOneExon(start, start+int64(len(seq))-1, t) {
		return "" // no intron between the ends
	}
	lpad := (delStartCDS - 1) % 3
	rpad := (3 - delEndCDS%3) % 3
	if delEndCDS+rpad > n {
		rpad = n - delEndCDS
	}
	left := t.CDSSequence[delStartCDS-1-lpad : delStartCDS-1]
	right := t.CDSSequence[delEndCDS : delEndCDS+rpad]
	core := seq
	if t.IsReverseStrand() {
		core = ReverseComplement(core)
	}
	alt := lower(left + right)
	if alt == "" {
		alt = "-"
	}
	return lower(left) + upper(core) + lower(right) + "/" + alt
}

// deletedBases returns the bases a variant removes and the genomic position of
// the first of them. A VCF-style deletion carries an anchor base that both
// alleles share and that nothing happens to; a MAF-style one, which is what
// genome-nexus sends, carries only the deleted bases.
func deletedBases(v *vcf.Variant) (seq string, start int64) {
	shared := 0
	for shared < len(v.Ref) && shared < len(v.Alt) && v.Ref[shared] == v.Alt[shared] {
		shared++
	}
	if len(v.Alt)-shared > 0 {
		return "", 0 // delins, not a pure deletion
	}
	return v.Ref[shared:], v.Pos + int64(shared)
}

// withinOneExon reports whether a genomic range lies inside a single exon, in
// which case the deleted bases and the deleted CDS bases are the same sequence.
func withinOneExon(start, end int64, t *cache.Transcript) bool {
	for i := range t.Exons {
		if start >= t.Exons[i].Start && end <= t.Exons[i].End {
			return true
		}
	}
	return false
}
