package annotate

import (
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// intronFlankBaseAt returns the base at a genomic position inside one of the
// cached intron flanks, in coding orientation. It reports false for an exonic
// position or one further into the intron than the flank reaches.
func intronFlankBaseAt(pos int64, t *cache.Transcript) (byte, bool) {
	idx := t.FindNearestExonIdx(pos)
	if idx < 0 {
		return 0, false
	}
	fwd := t.IsForwardStrand()
	for _, i := range [3]int{idx - 1, idx, idx + 1} {
		if i < 0 || i >= len(t.Exons) {
			continue
		}
		e := &t.Exons[i]
		if pos >= e.Start && pos <= e.End {
			return 0, false // exonic
		}
		// In coding orientation the flank 3' of an exon lies to the right of
		// it on the forward strand and to the left on the reverse strand.
		var after, before int64
		if fwd {
			after, before = pos-e.End, e.Start-pos
		} else {
			after, before = e.Start-pos, pos-e.End
		}
		if after > 0 {
			if b, ok := e.IntronAfterBase(int(after - 1)); ok {
				return b, true
			}
		}
		if before > 0 {
			// IntronBefore runs up to the exon, so the base adjacent to it is
			// the last one stored.
			if b, ok := e.IntronBeforeBase(int(e.IntronBeforeLen) - int(before)); ok {
				return b, true
			}
		}
	}
	return 0, false
}

// intronicInsertionDup names a pure insertion that sits inside an intron, with
// both flanking bases intronic.
//
// The CDS-based shift cannot reach these: it works on CDSSequence, which an
// intron is not part of, so an insertion between two intron bases was written
// out unshifted and never tested for duplication. VEP shifts it like any
// other, and the donor dinucleotide makes that common - inserting a T after
// the T of a GT donor gives USP8 c.541+2dup with hgvs_offset 1, where
// vibe-vep wrote c.541+1_541+2insT and no offset.
//
// The 100bp flank bounds how far this can shift, which is the cap it was
// sized for; intronFlankBaseAt simply stops when it runs out of sequence.
func intronicInsertionDup(v *vcf.Variant, t *cache.Transcript, insertedSeq string, result *ConsequenceResult) (string, bool) {
	n := len(insertedSeq)
	if n == 0 || n > 64 {
		return "", false
	}
	// The insertion sits between these two bases, named in coding order.
	after, before := v.Pos, v.Pos+1
	step := int64(1)
	if t.IsReverseStrand() {
		after, before = v.Pos+1, v.Pos
		step = -1
	}
	if _, ok := intronFlankBaseAt(after, t); !ok {
		return "", false
	}
	if _, ok := intronFlankBaseAt(before, t); !ok {
		return "", false
	}

	var seq [64]byte
	copy(seq[:n], insertedSeq)

	// 3' shift: while the inserted sequence starts with the base it would be
	// placed in front of, moving it one base along describes the same change.
	shift := 0
	for {
		b, ok := intronFlankBaseAt(before, t)
		if !ok || b != seq[0] {
			break
		}
		copy(seq[:n-1], seq[1:n])
		seq[n-1] = b
		after, before = after+step, before+step
		shift++
	}

	// A duplication is an insertion whose bases repeat the ones immediately
	// before it, read back from the insertion point in coding order.
	for i := 0; i < n; i++ {
		b, ok := intronFlankBaseAt(after-int64(i)*step, t)
		if !ok || b != seq[n-1-i] {
			return "", false
		}
	}

	dupStart := after - int64(n-1)*step
	startStr := genomicToHGVScPos(dupStart, t)
	endStr := genomicToHGVScPos(after, t)
	if startStr == "" || endStr == "" {
		return "", false
	}
	result.HGVSOffset = shift
	if startStr == endStr {
		return "c." + startStr + "dup", true
	}
	return "c." + startStr + "_" + endStr + "dup", true
}

// refAlleleBase reads a base out of a variant's own reference allele, in coding
// orientation. A deletion's bases are all spelled out there, so this reaches
// positions the transcript cache does not hold.
func refAlleleBase(pos int64, v *vcf.Variant, t *cache.Transcript) (byte, bool) {
	i := pos - v.Pos
	if i < 0 || int(i) >= len(v.Ref) {
		return 0, false
	}
	b := v.Ref[i]
	if t.IsReverseStrand() {
		b = Complement(b)
	}
	return b, true
}

// codingBaseAt returns the base at a genomic position in coding orientation,
// from the CDS where the position is exonic and from the cached intron flanks
// otherwise. It reports false wherever vibe-vep has no sequence - past the
// flank, or in a UTR, which the UTR-specific shifts cover.
func codingBaseAt(pos int64, t *cache.Transcript) (byte, bool) {
	if cds := GenomicToCDS(pos, t); cds >= 1 && int(cds) <= len(t.CDSSequence) {
		return t.CDSSequence[cds-1], true
	}
	if t.FindExon(pos) != nil {
		// Exonic but outside the CDS, so a UTR base. Both tails are stored in
		// coding orientation and addressed by distance from the start or stop
		// codon, which is the same number HGVS writes as -k or *k. A deletion
		// running off the end of the CDS needs these to shift: VEP moves
		// CALR c.1216_*32del three bases to c.1219_*35del.
		past3 := pos > t.CDSEnd
		if t.IsReverseStrand() {
			past3 = pos < t.CDSStart
		}
		if past3 {
			d := exonicDistance(t.CDSEnd, pos, t)
			if t.IsReverseStrand() {
				d = exonicDistance(pos, t.CDSStart, t)
			}
			if d >= 1 && int(d) <= len(t.UTR3Sequence) {
				return t.UTR3Sequence[d-1], true
			}
			return 0, false
		}
		k := exonicDistance(pos, t.CDSStart, t)
		if t.IsReverseStrand() {
			k = exonicDistance(t.CDSEnd, pos, t)
		}
		if k >= 1 {
			return t.UTR5BaseAt(int(k))
		}
		return 0, false
	}
	return intronFlankBaseAt(pos, t)
}

// genomicDeletionShift returns how far a deletion 3'-shifts, counted in bases
// along the transcript, working directly on genomic positions.
//
// A deletion may move one base 3' whenever its first deleted base equals the
// base just past its end, which is one rule whether those two bases sit in an
// exon, in an intron, or one of each. The shifts this backs up - the donor
// flank, the acceptor flank, the 5'UTR - each describe one such arrangement,
// and a deletion in none of those arrangements was simply not shifted: of the
// 874 rows where vibe-vep and VEP disagree on a deletion's HGVSc, 790 have
// vibe-vep reporting no shift at all against VEP's 1 to 7 bases.
//
// The 3'-most base it has to compare always lies just past the deletion, so it
// comes from the cache. The 5'-most one is a base the deletion removes, and
// deleted carries those: the variant's own reference allele spells out every
// one of them, however deep in an intron it sits. Without it a deletion
// starting past the 100bp flank could not shift at all, which is why STK11
// c.921-473_977del, whose first deleted base is 481 bases into the intron and
// whose eight-base shift lands entirely inside the next exon, came out
// unshifted. deleted may be nil, and the cache is still asked first.
func genomicDeletionShift(delStart, delEnd int64, t *cache.Transcript, deleted func(int64) (byte, bool)) int {
	// formatHGVScDeletion hands these over in transcript order, so on the
	// reverse strand the first is the larger. Put them back in genomic order
	// before picking the ends -- the early return this replaces meant no
	// reverse-strand deletion ever reached the loop below.
	if delStart > delEnd {
		delStart, delEnd = delEnd, delStart
	}
	// head is the 5'-most deleted base in transcript order, tail the 3'-most.
	head, tail, step := delStart, delEnd, int64(1)
	if t.IsReverseStrand() {
		head, tail, step = delEnd, delStart, -1
	}
	maxShift := len(t.CDSSequence) + 2*cache.IntronFlankBases
	shift := 0
	for shift < maxShift {
		hb, ok := codingBaseAt(head, t)
		if !ok && deleted != nil {
			hb, ok = deleted(head)
		}
		if !ok {
			break
		}
		nb, ok := codingBaseAt(tail+step, t)
		if !ok || hb != nb {
			break
		}
		head, tail = head+step, tail+step
		shift++
	}
	return shift
}
