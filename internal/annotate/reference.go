package annotate

// Reference supplies reference bases, 1-based and inclusive, upper case.
//
// The transcript cache answers for the coding sequence, both UTR tails and a
// hundred bases of each intron flank, and the shifts read through it. Where
// it has nothing -- deeper in an intron than the stored flank, or anywhere at
// all on a transcript with no coding sequence -- a shift used to stop,
// because there was no third place to ask. This is that third place.
//
// It is optional. Without one the shifts behave exactly as they did before,
// which is right as far as the cache reaches and short of VEP beyond it.
type Reference interface {
	Range(chrom string, start, end int64) ([]byte, error)
}

// referenceBase reads one base from the reference, or reports false where the
// reference has none: an assembly without that chromosome, a position past
// its end, or a gap, which no shift should walk into.
func referenceBase(ref Reference, chrom string, pos int64) (byte, bool) {
	if ref == nil || pos < 1 {
		return 0, false
	}
	b, err := ref.Range(chrom, pos, pos)
	if err != nil || len(b) != 1 || b[0] == 'N' {
		return 0, false
	}
	return b[0], true
}

// firstReference picks the reference out of an optional argument list, so a
// caller with none writes nothing and one holding a nil interface is the same
// as a caller with none.
func firstReference(refs []Reference) Reference {
	for _, r := range refs {
		if r != nil {
			return r
		}
	}
	return nil
}

// windowBases is how much reference a buffered read fetches at once. A shift
// walks a run of adjacent positions, so one read of a few hundred bases
// covers the whole walk; asking for them one at a time is a system call each,
// which on a box with no page cache left is a disk seek each.
const windowBases = 256

// bufferedReference reads the reference in windows and answers from them.
//
// It is created per shift and used by one goroutine, so it holds its window
// without a lock. Walks run forward on a forward-strand transcript and
// backward on a reverse-strand one, so the window is placed around the
// position asked for rather than ahead of it.
type bufferedReference struct {
	ref    Reference
	chrom  string
	lo, hi int64 // 1-based inclusive bounds of buf, lo > hi when empty
	buf    []byte
}

// buffered wraps a reference so a walk reads it in windows. A nil reference
// stays nil, so callers keep testing for one the same way.
func buffered(ref Reference) Reference {
	if ref == nil {
		return nil
	}
	return &bufferedReference{ref: ref, lo: 1, hi: 0}
}

func (b *bufferedReference) Range(chrom string, start, end int64) ([]byte, error) {
	if chrom != b.chrom || start < b.lo || end > b.hi {
		lo := start - windowBases/4
		if lo < 1 {
			lo = 1
		}
		hi := max(end, lo+windowBases-1)
		buf, err := b.ref.Range(chrom, lo, hi)
		if err != nil {
			return nil, err
		}
		b.chrom, b.lo, b.buf = chrom, lo, buf
		b.hi = lo + int64(len(buf)) - 1
		if start < b.lo || end > b.hi {
			// Clipped at a chromosome end, so answer from the source rather
			// than from a window that does not reach.
			return b.ref.Range(chrom, start, end)
		}
	}
	return b.buf[start-b.lo : end-b.lo+1], nil
}
