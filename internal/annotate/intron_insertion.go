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
