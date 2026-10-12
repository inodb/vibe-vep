package annotate

import (
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// startWindowBack and startWindowFwd bound the sequence examined around the
// start codon. Three bases either side is enough to see an ATG that has slid
// by a base or two; the window is widened a little so an indel cannot shorten
// it below that.
const (
	startWindowBack = 6
	startWindowFwd  = 9
)

// addStartRetainedWhenTheStartReforms adds start_retained_variant where a
// frameshift has destroyed the start codon but the sequence still presents an
// ATG across the same place.
//
// The two terms read as a contradiction and VEP only ever emits
// start_retained_variant in this pairing -- 16 rows in a VEP111 MSK-IMPACT
// MAF, all frameshift_variant,start_lost,start_retained_variant. What decides
// it is the base in front of the start codon. ERBB4 ENST00000342788 and KDM5A
// ENST00000399788 both have an A at c.-1, so deleting the A of ATG leaves
// ...A|TG and an ATG reappears one base upstream; CUL3 ENST00000264414 has a
// C there, nothing reforms, and VEP gives the same c.1del no such term. An
// inserted A at c.1_2 does the same from the other side, making AATG.
func addStartRetainedWhenTheStartReforms(v *vcf.Variant, t *cache.Transcript, result *ConsequenceResult) {
	if !hasTerm(result.Consequence, ConsequenceFrameshiftVariant) ||
		!hasTerm(result.Consequence, ConsequenceStartLost) ||
		hasTerm(result.Consequence, ConsequenceStartRetained) {
		return
	}
	if !startCodonReforms(v, t) {
		return
	}
	result.Consequence = joinConsequenceTerms(
		append(splitConsequenceTerms(result.Consequence), ConsequenceStartRetained))
}

// startCodonReforms applies the variant to the sequence around the start codon
// and reports whether an ATG still sits across where the old one was.
func startCodonReforms(v *vcf.Variant, t *cache.Transcript) bool {
	if !t.IsProteinCoding() || len(t.CDSSequence) < startWindowFwd {
		return false
	}

	// The window runs on one axis: c.-startWindowBack .. c.startWindowFwd.
	var win []byte
	for k := startWindowBack; k >= 1; k-- {
		b, ok := t.UTR5BaseAt(k)
		if !ok {
			return false // no 5' sequence to reform a start out of
		}
		win = append(win, b)
	}
	anchor := len(win) // index of c.1 in win
	win = append(win, t.CDSSequence[:startWindowFwd]...)

	lo, hi, ok := variantAxisSpan(v, t)
	if !ok {
		return false
	}
	// The A of the ATG has to be the base disturbed. Everything downstream of
	// it leaves the start codon standing where it was, and VEP gives none of
	// those the term however the frame ends up -- CBFB c.3dup and PIK3C3
	// c.7del are both plain frameshift_variant,start_lost.
	if lo > 1 || hi < 1 {
		return false
	}
	alt := v.Alt
	if t.IsReverseStrand() {
		alt = ReverseComplement(alt)
	}

	// Rebuild the window with the variant applied. An axis position p sits at
	// win[anchor+p-1] for p >= 1 and win[anchor+p] for p <= -1.
	idx := func(p int64) int {
		if p >= 1 {
			return anchor + int(p) - 1
		}
		return anchor + int(p)
	}
	first, last := idx(lo), idx(hi)
	if len(v.Ref) == 0 {
		// A pure insertion falls between lo and hi rather than over them.
		last = first
		first++
	}
	if first < 0 || last >= len(win) || first > last+1 {
		return false
	}
	altered := make([]byte, 0, len(win)+len(alt))
	altered = append(altered, win[:first]...)
	altered = append(altered, alt...)
	altered = append(altered, win[last+1:]...)

	// An ATG counts only where it still covers the old start codon, which ran
	// from anchor to anchor+2 before the edit, and only where it has actually
	// moved: one left where c.-1 supplies the A, one right where an A was
	// inserted behind it.
	for i := anchor - 2; i <= anchor+2; i++ {
		if i == anchor || i < 0 || i+3 > len(altered) {
			continue
		}
		if altered[i] == 'A' && altered[i+1] == 'T' && altered[i+2] == 'G' {
			return true
		}
	}
	return false
}

// variantAxisSpan returns the variant's reference bases as positions on the
// start-codon axis, in transcript order: 1,2,3... for coding bases and
// -1,-2... for the 5'UTR bases in front of them. A pure insertion reports the
// two positions it falls between.
func variantAxisSpan(v *vcf.Variant, t *cache.Transcript) (lo, hi int64, ok bool) {
	first, last := v.Pos, v.Pos+int64(len(v.Ref))-1
	if len(v.Ref) == 0 {
		last = v.Pos + 1
	}
	if t.IsReverseStrand() {
		first, last = last, first
	}
	a, ok1 := startAxisPos(first, t)
	b, ok2 := startAxisPos(last, t)
	if !ok1 || !ok2 {
		return 0, 0, false
	}
	return a, b, true
}

// startAxisPos places a genomic position on the start-codon axis.
func startAxisPos(pos int64, t *cache.Transcript) (int64, bool) {
	if cds := GenomicToCDS(pos, t); cds >= 1 {
		return cds, true
	}
	if t.FindExon(pos) == nil {
		return 0, false
	}
	// Exonic but not coding: a 5'UTR base is c.-k, k bases from the CDS start.
	k := exonicDistance(pos, t.CDSStart, t)
	if t.IsReverseStrand() {
		k = exonicDistance(t.CDSEnd, pos, t)
	}
	if k < 1 {
		return 0, false
	}
	return -k, true
}
