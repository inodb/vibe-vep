package output

import (
	"strconv"

	"github.com/inodb/vibe-vep/internal/annotate"
)

// GenomeSequence supplies reference bases, 1-based and inclusive, upper case,
// with unknown bases as N. It is what hgvsg needs and the transcript cache
// cannot give: the 3' shift runs along the chromosome, through introns and
// past the ends of every transcript.
type GenomeSequence interface {
	Range(chrom string, start, end int64) ([]byte, error)
}

// normalizeHGVSg names a change the way HGVS requires, which is not the way
// it was reported. Two rules apply, and VEP applies both:
//
//   - a change in a repeat is written at the 3'-most position the reference
//     allows, so deleting one G from GGG names the last G;
//   - an insertion that repeats what already precedes it is a duplication,
//     so inserting G after a G is "dup", not "insG".
//
// It returns the empty string when the change cannot be normalized -- no
// reference to read, or a reference allele that does not match the genome --
// and the caller then keeps the name built from the query's own coordinates.
func normalizeHGVSg(chrom string, start, end int64, ref, alt string, g GenomeSequence) string {
	if g == nil || start < 1 {
		return ""
	}
	if isEmptyAllele(ref) {
		ref = ""
	}
	if isEmptyAllele(alt) {
		alt = ""
	}

	// The replaced interval, which is empty for an insertion: s = e+1 says
	// the change sits between e and s without covering either.
	s, e := start, end
	if ref == "" {
		s, e = start+1, start
	} else {
		e = start + int64(len(ref)) - 1
	}

	// A caller writing VCF-style alleles anchors an indel on the base before
	// it (GG>G), and one writing MAF-style does not (G>-). Trimming what the
	// two alleles share turns the first into the second and leaves the second
	// alone, so the rest of this works on one convention.
	for len(ref) > 0 && len(alt) > 0 && ref[0] == alt[0] {
		ref, alt = ref[1:], alt[1:]
		s++
	}
	for len(ref) > 0 && len(alt) > 0 && ref[len(ref)-1] == alt[len(alt)-1] {
		ref, alt = ref[:len(ref)-1], alt[:len(alt)-1]
		e--
	}
	if ref == "" && alt == "" {
		return ""
	}

	// Shifting a change the reference does not actually carry would move it
	// somewhere arbitrary. VEP refuses such a variant outright; we decline to
	// rename it and leave the query's own coordinates standing.
	if ref != "" {
		got, err := g.Range(chrom, s, e)
		if err != nil || string(got) != ref {
			return ""
		}
	}

	switch {
	case alt == "":
		s, e = shiftDeletion(chrom, s, e, ref, g)
		return formatHGVSg(chrom, s, e, ref, "-")

	case ref == "":
		// The insertion sits after anchor, between it and the next base.
		anchor, alt := shiftInsertion(chrom, s-1, alt, g)
		prefix := chrom + ":g."
		if prev, err := g.Range(chrom, anchor-int64(len(alt))+1, anchor); err == nil && string(prev) == alt {
			if len(alt) == 1 {
				return prefix + strconv.FormatInt(anchor, 10) + "dup"
			}
			return prefix + strconv.FormatInt(anchor-int64(len(alt))+1, 10) + "_" +
				strconv.FormatInt(anchor, 10) + "dup"
		}
		return prefix + strconv.FormatInt(anchor, 10) + "_" +
			strconv.FormatInt(anchor+1, 10) + "ins" + alt

	case len(ref) > 1 && len(ref) == len(alt) && annotate.ReverseComplement(ref) == alt:
		// A stretch replaced by its own reverse complement is an inversion,
		// which HGVS names as such rather than as a delins.
		return chrom + ":g." + strconv.FormatInt(s, 10) + "_" + strconv.FormatInt(e, 10) + "inv"

	default:
		return formatHGVSg(chrom, s, e, ref, alt)
	}
}

// shiftDeletion walks a deleted span 3' for as long as the base just past it
// repeats the first base it deletes.
func shiftDeletion(chrom string, s, e int64, ref string, g GenomeSequence) (int64, int64) {
	w := refWalker{g: g, chrom: chrom, pos: e + 1}
	for i := 0; ; i++ {
		b, ok := w.at(i)
		if !ok || b != ref[0] {
			return s, e
		}
		ref = ref[1:] + ref[:1]
		s, e = s+1, e+1
	}
}

// shiftInsertion walks an insertion 3' for as long as the base just past it
// repeats the first base it inserts, rotating the inserted sequence to match.
// anchor is the base the insertion follows.
func shiftInsertion(chrom string, anchor int64, alt string, g GenomeSequence) (int64, string) {
	w := refWalker{g: g, chrom: chrom, pos: anchor + 1}
	for i := 0; ; i++ {
		b, ok := w.at(i)
		if !ok || b != alt[0] {
			return anchor, alt
		}
		alt = alt[1:] + alt[:1]
		anchor++
	}
}

// refWalker reads the reference forward from a position, fetching more only
// when the shift actually gets that far. Most shifts stop within a base or
// two, and the few that run on are repeats, so the window grows rather than
// being guessed at up front.
type refWalker struct {
	g     GenomeSequence
	chrom string
	pos   int64 // 1-based position of buf[0]
	buf   []byte
	done  bool
}

// at returns the i-th base from pos, and whether the reference has one.
func (w *refWalker) at(i int) (byte, bool) {
	const window = 256
	for i >= len(w.buf) {
		if w.done {
			return 0, false
		}
		next, err := w.g.Range(w.chrom, w.pos+int64(len(w.buf)), w.pos+int64(len(w.buf))+window-1)
		if err != nil || len(next) == 0 {
			w.done = true
			return 0, false
		}
		if len(next) < window {
			w.done = true // the chromosome ended
		}
		w.buf = append(w.buf, next...)
	}
	// A shift must stop at a gap rather than walk into one.
	if w.buf[i] == 'N' {
		return 0, false
	}
	return w.buf[i], true
}
