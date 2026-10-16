package annotate

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestParseCDSCoord(t *testing.T) {
	cases := []struct {
		in   string
		want cdsCoord
		ok   bool
	}{
		{"1799", cdsCoord{Pos: 1799}, true},
		{"-29", cdsCoord{Pos: -29}, true},
		{"*8", cdsCoord{Pos: 8, UTR3: true}, true},
		{"1899-1", cdsCoord{Pos: 1899, Offset: -1}, true},
		{"1946+1", cdsCoord{Pos: 1946, Offset: 1}, true},
		{"-29+2", cdsCoord{Pos: -29, Offset: 2}, true},
		{"*8+2", cdsCoord{Pos: 8, UTR3: true, Offset: 2}, true},
		{"-29-2", cdsCoord{Pos: -29, Offset: -2}, true},
		{"", cdsCoord{}, false},
		{"abc", cdsCoord{}, false},
		{"12+", cdsCoord{}, false},
		{"*", cdsCoord{}, false},
	}
	for _, c := range cases {
		got, ok := parseCDSCoord(c.in)
		if ok != c.ok {
			t.Errorf("parseCDSCoord(%q) ok = %v, want %v", c.in, ok, c.ok)
			continue
		}
		if ok && got != c.want {
			t.Errorf("parseCDSCoord(%q) = %+v, want %+v", c.in, got, c.want)
		}
	}
}

// twoExonTranscript builds a transcript with two exons and a coding sequence
// that starts inside the first and ends inside the second, so that the 5'UTR,
// the CDS, the intron and the 3'UTR are all reachable.
//
//	exon 1: 100-119   exon 2: 200-219
//	CDS: 110 .. 209
func twoExonTranscript(strand int8) *cache.Transcript {
	t := &cache.Transcript{
		ID:     "ENSTTEST",
		Chrom:  "1",
		Strand: strand,
		Exons: []cache.Exon{
			{Number: 1, Start: 100, End: 119},
			{Number: 2, Start: 200, End: 219},
		},
		CDSStart: 110,
		CDSEnd:   209,
	}
	for i := range t.Exons {
		e := &t.Exons[i]
		if e.End >= t.CDSStart && e.Start <= t.CDSEnd {
			e.CDSStart, e.CDSEnd = max64(e.Start, t.CDSStart), min64(e.End, t.CDSEnd)
		}
	}
	t.BuildCDSIndex()
	return t
}

func TestTranscriptPosToGenomicRoundTrips(t *testing.T) {
	for _, strand := range []int8{1, -1} {
		tx := twoExonTranscript(strand)
		for _, pos := range []int64{100, 110, 119, 200, 209, 219} {
			txPos := GenomicToTranscriptPos(pos, tx)
			require.NotZero(t, txPos, "strand %d: %d has no transcript position", strand, pos)
			back := transcriptPosToGenomic(txPos, tx)
			assert.Equal(t, pos, back, "strand %d: %d -> tx %d -> %d", strand, pos, txPos, back)
		}
		// Positions off either end of the transcript have no genomic answer.
		assert.Zero(t, transcriptPosToGenomic(0, tx))
		assert.Zero(t, transcriptPosToGenomic(41, tx))
	}
}

func TestCDSCoordToGenomicForward(t *testing.T) {
	tx := twoExonTranscript(1)

	cases := []struct {
		name string
		in   string
		want int64
	}{
		{"first coding base", "1", 110},
		{"last base of exon 1", "10", 119},
		{"first base of exon 2", "11", 200},
		{"one base into the intron, from the exon before", "10+1", 120},
		{"one base into the intron, from the exon after", "11-1", 199},
		{"five bases into the intron", "10+5", 124},
		{"one base before the coding sequence", "-1", 109},
		{"ten bases before, still in exon 1", "-10", 100},
		{"one base past the coding sequence", "*1", 210},
		{"ten bases past, still in exon 2", "*10", 219},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			coord, ok := parseCDSCoord(c.in)
			require.True(t, ok, "parse %q", c.in)
			got, err := cdsCoordToGenomic(coord, tx)
			require.NoError(t, err)
			assert.Equal(t, c.want, got)
		})
	}
}

func TestCDSCoordToGenomicReverse(t *testing.T) {
	// On the reverse strand the transcript runs from high coordinates to low,
	// so every answer is mirrored and the intron offsets change direction.
	tx := twoExonTranscript(-1)

	cases := []struct {
		name string
		in   string
		want int64
	}{
		{"first coding base is the highest coordinate", "1", 209},
		{"last base of the first exon", "10", 200},
		{"first base of the second exon", "11", 119},
		{"one base into the intron, from the exon before", "10+1", 199},
		{"one base into the intron, from the exon after", "11-1", 120},
		{"one base before the coding sequence", "-1", 210},
		{"one base past the coding sequence", "*1", 109},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			coord, ok := parseCDSCoord(c.in)
			require.True(t, ok, "parse %q", c.in)
			got, err := cdsCoordToGenomic(coord, tx)
			require.NoError(t, err)
			assert.Equal(t, c.want, got)
		})
	}
}

// Every exonic position should survive a trip out through the forward
// implementation and back, which is the property the whole inverse rests on.
func TestForwardAndInverseAgree(t *testing.T) {
	for _, strand := range []int8{1, -1} {
		tx := twoExonTranscript(strand)
		for pos := int64(100); pos <= 219; pos++ {
			if pos > 119 && pos < 200 {
				continue // intronic; covered by the explicit cases above
			}
			written := genomicToHGVScPos(pos, tx)
			if written == "" {
				continue
			}
			coord, ok := parseCDSCoord(written)
			require.True(t, ok, "strand %d: cannot parse %q written for %d", strand, written, pos)
			back, err := cdsCoordToGenomic(coord, tx)
			require.NoError(t, err, "strand %d: %d written as %q", strand, pos, written)
			assert.Equal(t, pos, back, "strand %d: %d written as %q came back as %d", strand, pos, written, back)
		}
	}
}

func TestCDSCoordToGenomicRejectsOutOfRange(t *testing.T) {
	tx := twoExonTranscript(1)
	for _, in := range []string{"999", "-999", "*999"} {
		coord, ok := parseCDSCoord(in)
		require.True(t, ok)
		_, err := cdsCoordToGenomic(coord, tx)
		assert.Error(t, err, "c.%s should not resolve", in)
	}
}

// Guard the guard: if the forward implementation stopped producing positions,
// TestForwardAndInverseAgree would pass by skipping everything.
func TestForwardWritesPositionsForTheTestTranscript(t *testing.T) {
	for _, strand := range []int8{1, -1} {
		tx := twoExonTranscript(strand)
		written := 0
		for pos := int64(100); pos <= 219; pos++ {
			if pos > 119 && pos < 200 {
				continue
			}
			if genomicToHGVScPos(pos, tx) != "" {
				written++
			}
		}
		assert.Equal(t, 40, written, "strand %d: every exonic base should get a position", strand)
	}
}
