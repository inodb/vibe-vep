package output

import (
	"errors"
	"testing"
)

// fakeGenome serves one stretch of one chromosome, placed at a real
// coordinate so the expected names can be the ones VEP actually returns.
type fakeGenome struct {
	chrom string
	first int64 // 1-based position of seq[0]
	seq   string
}

func (f fakeGenome) Range(chrom string, start, end int64) ([]byte, error) {
	if chrom != f.chrom {
		return nil, errors.New("no such chromosome")
	}
	last := f.first + int64(len(f.seq)) - 1
	if start < f.first {
		start = f.first
	}
	if end > last {
		end = last
	}
	if start > end {
		return nil, nil
	}
	return []byte(f.seq[start-f.first : end-f.first+1]), nil
}

func TestNormalizeHGVSgShiftsAndRenames(t *testing.T) {
	// 13:41134320-41134345 of GRCh37, which carries the GGG run at 330-332
	// that VEP shifts 13:41134330 G>- into.
	brca2 := fakeGenome{chrom: "13", first: 41134320, seq: "CTGGGGCAAAGGGCTCATGCTGGATT"}

	for _, tc := range []struct {
		name       string
		g          GenomeSequence
		chrom      string
		start, end int64
		ref, alt   string
		want       string
	}{
		{
			name: "deletion walks to the 3' end of its run",
			g:    brca2, chrom: "13", start: 41134330, end: 41134330, ref: "G", alt: "-",
			want: "13:g.41134332del",
		},
		{
			name: "deletion with nowhere to go keeps its place",
			g:    brca2, chrom: "13", start: 41134333, end: 41134333, ref: "C", alt: "-",
			want: "13:g.41134333del",
		},
		{
			name: "multi-base deletion shifts as a block",
			// 41134322-41134325 is GGGG, so deleting GG from 322 names 324_325.
			g: brca2, chrom: "13", start: 41134322, end: 41134323, ref: "GG", alt: "-",
			want: "13:g.41134324_41134325del",
		},
		{
			name: "insertion repeating what precedes it is a duplication",
			// inserting G after 41134329 shifts to 332 and duplicates it.
			g: brca2, chrom: "13", start: 41134329, end: 41134330, ref: "-", alt: "G",
			want: "13:g.41134332dup",
		},
		{
			name: "insertion that repeats nothing stays an insertion",
			// 41134334 is T and 41134335 is C, so an inserted A neither
			// shifts nor duplicates.
			g: brca2, chrom: "13", start: 41134334, end: 41134335, ref: "-", alt: "A",
			want: "13:g.41134334_41134335insA",
		},
		{
			name: "insertion rotates as it shifts",
			// 327 is A, so inserting AG after 326 shifts one base and the
			// inserted pair rotates to GA. 328 is A again but the rotated
			// sequence now leads with G, so it goes no further.
			g: brca2, chrom: "13", start: 41134326, end: 41134327, ref: "-", alt: "AG",
			want: "13:g.41134327_41134328insGA",
		},
		{
			name: "an allele that changes nothing is not renamed",
			g:    brca2, chrom: "13", start: 41134336, end: 41134337, ref: "AT", alt: "AT",
			want: "",
		},
		{
			name: "substitution is left alone",
			g:    brca2, chrom: "13", start: 41134330, end: 41134330, ref: "G", alt: "A",
			want: "13:g.41134330G>A",
		},
		{
			name: "VCF-style anchor base is trimmed before shifting",
			// 41134329_41134330 GG>G is the same deletion as 41134330 G>-.
			g: brca2, chrom: "13", start: 41134329, end: 41134330, ref: "GG", alt: "G",
			want: "13:g.41134332del",
		},
		{
			name: "a reference allele the genome does not carry is not renamed",
			g:    brca2, chrom: "13", start: 41134330, end: 41134330, ref: "A", alt: "-",
			want: "",
		},
		{
			name: "a chromosome the reference does not hold is not renamed",
			g:    brca2, chrom: "99", start: 41134330, end: 41134330, ref: "G", alt: "-",
			want: "",
		},
		{
			name: "without a reference nothing is renamed",
			g:    nil, chrom: "13", start: 41134330, end: 41134330, ref: "G", alt: "-",
			want: "",
		},
	} {
		t.Run(tc.name, func(t *testing.T) {
			got := normalizeHGVSg(tc.chrom, tc.start, tc.end, tc.ref, tc.alt, tc.g)
			if got != tc.want {
				t.Errorf("normalizeHGVSg = %q, want %q", got, tc.want)
			}
		})
	}
}

func TestNormalizeHGVSgNamesAnInversion(t *testing.T) {
	g := fakeGenome{chrom: "X", first: 76938200, seq: "CCCCCCCAATTTTTTT"}
	// 76938207_76938208 is AA, replaced by TT, its own reverse complement.
	if got := normalizeHGVSg("X", 76938207, 76938208, "AA", "TT", g); got != "X:g.76938207_76938208inv" {
		t.Errorf("normalizeHGVSg = %q, want an inversion", got)
	}
}

func TestNormalizeHGVSgStopsAtAGap(t *testing.T) {
	// A run of unknown bases is not a repeat to shift through.
	g := fakeGenome{chrom: "1", first: 100, seq: "ATNNNNNNNN"}
	if got := normalizeHGVSg("1", 101, 101, "T", "-", g); got != "1:g.101del" {
		t.Errorf("normalizeHGVSg = %q, want the unshifted name", got)
	}
}

func TestNormalizeHGVSgShiftsPastTheFirstWindow(t *testing.T) {
	// 600 A's: the shift must keep fetching rather than stop at the window.
	seq := make([]byte, 0, 601)
	for i := 0; i < 600; i++ {
		seq = append(seq, 'A')
	}
	seq = append(seq, 'C')
	g := fakeGenome{chrom: "1", first: 1, seq: string(seq)}
	if got := normalizeHGVSg("1", 1, 1, "A", "-", g); got != "1:g.600del" {
		t.Errorf("normalizeHGVSg = %q, want 1:g.600del", got)
	}
}
