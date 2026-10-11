package annotate

import (
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Additional SO terms VEP reports for positions inside an intron near a splice
// site. vibe-vep previously collapsed all of these into splice_region_variant
// or intron_variant.
const (
	ConsequenceSpliceDonor5thBase  = "splice_donor_5th_base_variant"
	ConsequenceSpliceDonorRegion   = "splice_donor_region_variant"
	ConsequencePolypyrimidineTract = "splice_polypyrimidine_tract_variant"
	ConsequenceProteinAltering     = "protein_altering_variant"
)

// consequenceRank orders terms the way Ensembl ranks consequence severity,
// which is the order VEP emits them in. Only the terms vibe-vep can produce
// need an entry; anything unknown sorts last but keeps a stable relative order.
var consequenceRank = map[string]int{
	ConsequenceSpliceAcceptor:        2,
	ConsequenceSpliceDonor:           3,
	ConsequenceStopGained:            4,
	ConsequenceFrameshiftVariant:     5,
	ConsequenceStopLost:              6,
	ConsequenceStartLost:             7,
	ConsequenceInframeInsertion:      9,
	ConsequenceInframeDeletion:       10,
	ConsequenceMissenseVariant:       11,
	ConsequenceProteinAltering:       12,
	ConsequenceSpliceDonor5thBase:    13,
	ConsequenceSpliceRegion:          14,
	ConsequenceSpliceDonorRegion:     15,
	ConsequencePolypyrimidineTract:   16,
	ConsequenceStartRetained:         18,
	ConsequenceStopRetained:          19,
	ConsequenceSynonymousVariant:     20,
	ConsequenceCodingSequenceVariant: 21,
	ConsequenceMatureMiRNA:           22,
	Consequence5PrimeUTR:             23,
	Consequence3PrimeUTR:             24,
	ConsequenceNonCodingExon:         25,
	ConsequenceIntronVariant:         26,
	ConsequenceUpstreamGene:          29,
	ConsequenceDownstreamGene:        30,
	ConsequenceIntergenicVariant:     31,
}

func rankOf(term string) int {
	if r, ok := consequenceRank[term]; ok {
		return r
	}
	return 100
}

// intronOffset returns the signed distance from pos to the nearest exon
// boundary, in transcript orientation: positive N means the Nth base into the
// intron on the donor side, negative N the Nth base before the next exon on the
// acceptor side. ok is false when pos is not intronic for this transcript.
func intronOffset(pos int64, t *cache.Transcript) (offset int64, ok bool) {
	idx := t.FindNearestExonIdx(pos)
	if idx < 0 {
		return 0, false
	}
	best := int64(0)
	found := false
	for _, i := range [3]int{idx - 1, idx, idx + 1} {
		if i < 0 || i >= len(t.Exons) {
			continue
		}
		exon := &t.Exons[i]
		if pos >= exon.Start && pos <= exon.End {
			return 0, false // exonic
		}
		var d int64
		if pos > exon.End {
			d = pos - exon.End // downstream of this exon
			if !t.IsForwardStrand() {
				d = -d
			}
		} else {
			d = exon.Start - pos // upstream of this exon
			if t.IsForwardStrand() {
				d = -d
			}
		}
		if !found || abs64(d) < abs64(best) {
			best, found = d, true
		}
	}
	return best, found
}

// spliceTermsForOffset returns the SO terms VEP reports for a single intronic
// position at the given signed offset.
//
// Derived from single-base substitutions in a VEP111 MSK-IMPACT MAF:
//
//	+1,+2   splice_donor_variant                                  (no intron_variant)
//	+3,+4,+6 splice_donor_region_variant, intron_variant
//	+5      splice_donor_5th_base_variant, intron_variant         (not donor_region)
//	+7,+8   splice_region_variant, intron_variant
//	-1,-2   splice_acceptor_variant                               (no intron_variant)
//	-3..-8  splice_region_variant, splice_polypyrimidine_tract_variant, intron_variant
//	-9..-17 splice_polypyrimidine_tract_variant, intron_variant
//	beyond  intron_variant
func spliceTermsForOffset(offset int64) []string {
	switch {
	case offset == 1 || offset == 2:
		return []string{ConsequenceSpliceDonor}
	case offset == 5:
		return []string{ConsequenceSpliceDonor5thBase, ConsequenceIntronVariant}
	case offset >= 3 && offset <= 6:
		return []string{ConsequenceSpliceDonorRegion, ConsequenceIntronVariant}
	case offset >= 7 && offset <= 8:
		return []string{ConsequenceSpliceRegion, ConsequenceIntronVariant}
	case offset == -1 || offset == -2:
		return []string{ConsequenceSpliceAcceptor}
	case offset <= -3 && offset >= -8:
		return []string{ConsequenceSpliceRegion, ConsequencePolypyrimidineTract, ConsequenceIntronVariant}
	case offset <= -9 && offset >= -17:
		return []string{ConsequencePolypyrimidineTract, ConsequenceIntronVariant}
	default:
		return []string{ConsequenceIntronVariant}
	}
}

// spanSpliceTerms enumerates the SO terms for every genomic position a variant
// covers: intronic positions through spliceTermsForOffset, exonic ones as
// coding_sequence_variant.
//
// A deletion reaching from an exon into an intron earns the terms of both
// sides. VEP reports SMARCA4 c.4860_4864+34del as splice_donor_variant,
// splice_donor_5th_base_variant,coding_sequence_variant,intron_variant - one
// term per region the deletion touches, not just the site it breaks.
func spanSpliceTerms(v *vcf.Variant, t *cache.Transcript) []string {
	spanEnd := v.Pos + int64(len(v.Ref)) - 1
	if spanEnd < v.Pos {
		spanEnd = v.Pos
	}
	var terms []string
	for pos := v.Pos; pos <= spanEnd; pos++ {
		if !t.Contains(pos) {
			// Outside the transcript there is nothing to name. Without this
			// the nearest-exon lookup below answers with a large offset and
			// the span picks up an intron_variant it has not earned.
			continue
		}
		if t.FindExon(pos) == nil {
			if off, ok := intronOffset(pos, t); ok {
				terms = append(terms, spliceTermsForOffset(off)...)
			}
			continue
		}
		// An exonic base outside the CDS is a UTR base, and VEP names the UTR
		// it falls in: a deletion running from the 5'UTR into the coding
		// sequence is coding_sequence_variant,5_prime_UTR_variant.
		if GenomicToCDS(pos, t) >= 1 {
			terms = append(terms, ConsequenceCodingSequenceVariant)
		} else if term := utrTermAt(pos, t); term != "" {
			terms = append(terms, term)
		}
	}
	return terms
}

// utrTermAt names the UTR an exonic non-coding position falls in.
func utrTermAt(pos int64, t *cache.Transcript) string {
	if !t.IsProteinCoding() {
		return ""
	}
	if t.IsForwardStrand() {
		if pos < t.CDSStart {
			return Consequence5PrimeUTR
		}
		if pos > t.CDSEnd {
			return Consequence3PrimeUTR
		}
		return ""
	}
	if pos > t.CDSEnd {
		return Consequence5PrimeUTR
	}
	if pos < t.CDSStart {
		return Consequence3PrimeUTR
	}
	return ""
}

// exonicTerms are the SO terms that say a variant alters transcribed exonic
// sequence. VEP suppresses splice_polypyrimidine_tract_variant whenever one of
// these is present: of the 504 variants carrying the polypyrimidine term in a
// VEP111 MSK-IMPACT MAF, every one is intronic-only. A deletion running from
// the tract into the exon is reported as
// splice_acceptor_variant,coding_sequence_variant,intron_variant, dropping the
// tract term that the same span would earn on its own.
var exonicTerms = map[string]bool{
	ConsequenceCodingSequenceVariant: true,
	ConsequenceFrameshiftVariant:     true,
	ConsequenceInframeInsertion:      true,
	ConsequenceInframeDeletion:       true,
	ConsequenceMissenseVariant:       true,
	ConsequenceSynonymousVariant:     true,
	ConsequenceProteinAltering:       true,
	ConsequenceStopGained:            true,
	ConsequenceStopLost:              true,
	ConsequenceStopRetained:          true,
	ConsequenceStartLost:             true,
	ConsequenceStartRetained:         true,
	Consequence5PrimeUTR:             true,
	Consequence3PrimeUTR:             true,
	ConsequenceNonCodingExon:         true,
}

// joinConsequenceTerms renders a term set in Ensembl severity order, dropping
// duplicates. splice_region_variant is suppressed when a splice donor or
// acceptor term is present, which is what VEP does.
func joinConsequenceTerms(terms []string) string {
	seen := make(map[string]bool, len(terms))
	uniq := make([]string, 0, len(terms))
	hasSite := false
	has5thBase := false
	hasExonic := false
	for _, t := range terms {
		if t == "" || seen[t] {
			continue
		}
		switch t {
		case ConsequenceSpliceDonor, ConsequenceSpliceAcceptor:
			hasSite = true
		case ConsequenceSpliceDonor5thBase:
			has5thBase = true
		}
		if exonicTerms[t] {
			hasExonic = true
		}
		seen[t] = true
		uniq = append(uniq, t)
	}
	out := make([]string, 0, len(uniq))
	for _, t := range uniq {
		// splice_region_variant is not reported alongside a donor/acceptor
		// site, and splice_donor_5th_base_variant is the more specific term
		// where it overlaps splice_donor_region_variant. Both suppressions are
		// what VEP emits for spans covering these positions.
		if hasSite && t == ConsequenceSpliceRegion {
			continue
		}
		if has5thBase && t == ConsequenceSpliceDonorRegion {
			continue
		}
		// The polypyrimidine tract is only named for variants that stay
		// inside the intron; reaching the exon replaces it with the exonic
		// term. splice_donor_region_variant and the donor 5th base carry no
		// such rule and do co-occur with coding_sequence_variant.
		if hasExonic && t == ConsequencePolypyrimidineTract {
			continue
		}
		out = append(out, t)
	}
	// insertion sort by rank; the lists are tiny
	for i := 1; i < len(out); i++ {
		for j := i; j > 0 && rankOf(out[j]) < rankOf(out[j-1]); j-- {
			out[j], out[j-1] = out[j-1], out[j]
		}
	}
	res := ""
	for i, t := range out {
		if i > 0 {
			res += ","
		}
		res += t
	}
	return res
}

// isDelIns reports whether the variant both deletes and inserts sequence, as
// opposed to a pure insertion or pure deletion. Common prefix and suffix are
// trimmed first so a VCF-style anchored representation ("AG" -> "A") is still
// recognised as a pure deletion.
//
// VEP reports an in-frame delins as protein_altering_variant rather than
// inframe_deletion or inframe_insertion: across 776 in-frame delins variants in
// a VEP111 MSK-IMPACT MAF it says protein_altering_variant (642) or that plus a
// splice or stop term, and never inframe_deletion/inframe_insertion.
func isDelIns(ref, alt string) bool {
	i := 0
	for i < len(ref) && i < len(alt) && ref[i] == alt[i] {
		i++
	}
	r, a := ref[i:], alt[i:]
	j := 0
	for j < len(r) && j < len(a) && r[len(r)-1-j] == a[len(a)-1-j] {
		j++
	}
	return len(r)-j > 0 && len(a)-j > 0
}

// relabelInframeAsProteinAltering swaps an inframe_insertion/inframe_deletion
// term for protein_altering_variant, preserving any co-terms and their order.
func relabelInframeAsProteinAltering(consequence string) string {
	var terms []string
	start := 0
	for i := 0; i <= len(consequence); i++ {
		if i == len(consequence) || consequence[i] == ',' {
			t := consequence[start:i]
			if t == ConsequenceInframeInsertion || t == ConsequenceInframeDeletion {
				t = ConsequenceProteinAltering
			}
			terms = append(terms, t)
			start = i + 1
		}
	}
	return joinConsequenceTerms(terms)
}

// insertionSpliceTerms returns the SO terms for a pure insertion sitting
// between two intronic positions, or false when either flank is exonic.
//
// An insertion alters no base, so it only breaks a splice site when it lands
// strictly between the two essential ones. Between +1 and +2 it splits the
// donor GT and VEP says splice_donor_variant; one base further along, between
// +2 and +3, the GT is intact and VEP says splice_region_variant,intron_variant
// - 26 times against 3 that still say donor. The same holds at the acceptor:
// between -2 and -1 is splice_acceptor_variant, between -3 and -2 is the
// region plus the polypyrimidine tract.
//
// The finer donor terms go with the site terms. VEP names no
// splice_donor_region_variant or donor 5th base for an insertion, which is
// what vibe-vep was reporting for the +2/+3 case.
//
// This is the principle already applied where a junction insertion is treated
// as exonic: inserting between two bases leaves both of them where they are.
func insertionSpliceTerms(v *vcf.Variant, t *cache.Transcript) ([]string, bool) {
	if len(v.Ref) != 0 {
		return nil, false
	}
	after, before := v.Pos, v.Pos+1
	if !t.IsForwardStrand() {
		after, before = v.Pos+1, v.Pos
	}
	offAfter, okA := intronOffset(after, t)
	offBefore, okB := intronOffset(before, t)
	if !okA || !okB {
		return nil, false
	}
	if (offAfter == 1 && offBefore == 2) || (offAfter == -2 && offBefore == -1) {
		return spliceTermsForOffset(offAfter), true
	}
	terms := []string{ConsequenceIntronVariant}
	for _, off := range [2]int64{offAfter, offBefore} {
		for _, term := range spliceTermsForOffset(off) {
			switch term {
			case ConsequenceSpliceDonor, ConsequenceSpliceAcceptor,
				ConsequenceSpliceDonorRegion, ConsequenceSpliceDonor5thBase:
				terms = append(terms, ConsequenceSpliceRegion)
			default:
				terms = append(terms, term)
			}
		}
	}
	return terms, true
}
