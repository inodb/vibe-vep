package cache

import (
	"strings"
	"testing"
)

func TestPackDNA2BitRoundTrip(t *testing.T) {
	cases := []string{"A", "ACGT", "ACGTACGTA", strings.Repeat("ACGT", 25), "TTTT", "GATTACA"}
	for _, want := range cases {
		buf, n := PackDNA2Bit(want)
		if int(n) != len(want) {
			t.Fatalf("PackDNA2Bit(%q) length = %d, want %d", want, n, len(want))
		}
		e := &Exon{IntronAfterPacked: buf, IntronAfterLen: uint8(n)}
		var got strings.Builder
		for i := 0; i < int(n); i++ {
			b, ok := e.IntronAfterBase(i)
			if !ok {
				t.Fatalf("IntronAfterBase(%d) not ok for %q", i, want)
			}
			got.WriteByte(b)
		}
		if got.String() != want {
			t.Errorf("round trip of %q = %q", want, got.String())
		}
		// Four bases to a byte.
		if wantBytes := (len(want) + 3) / 4; len(buf) != wantBytes {
			t.Errorf("packed %q into %d bytes, want %d", want, len(buf), wantBytes)
		}
	}
}

func TestPackDNA2BitTruncatesAtUnknownBase(t *testing.T) {
	// An N cannot be encoded, and a shift cannot be validated across a base we
	// do not know, so the flank stops there rather than guessing.
	buf, n := PackDNA2Bit("ACGNTT")
	if n != 3 {
		t.Fatalf("length = %d, want 3 (truncated at N)", n)
	}
	e := &Exon{IntronAfterPacked: buf, IntronAfterLen: uint8(n)}
	if _, ok := e.IntronAfterBase(3); ok {
		t.Error("base past the N should not be readable")
	}
	for i, want := range []byte("ACG") {
		if b, _ := e.IntronAfterBase(i); b != want {
			t.Errorf("base %d = %c, want %c", i, b, want)
		}
	}
}

func TestIntronFlankBaseOutOfRange(t *testing.T) {
	e := &Exon{}
	if _, ok := e.IntronAfterBase(0); ok {
		t.Error("empty flank should report no base")
	}
	if _, ok := e.IntronBeforeBase(-1); ok {
		t.Error("negative index should report no base")
	}
}
