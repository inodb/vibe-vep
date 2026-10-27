package output

import "testing"

// A change crossing a splice site has one end in an intron and one in an
// exon. VEP names a residue only for the end that has one, and which end that
// is cannot be read off the codons, because such a change has none.
func TestCodingEndsOfHGVSc(t *testing.T) {
	for _, c := range []struct {
		hgvsc              string
		wantStart, wantEnd bool
	}{
		{"ENST00000252674.5:c.421-1_421delinsAA", false, true}, // reaches in from the intron
		{"ENST00000220592.6:c.1216_1269+211del", true, false},  // reaches out into the intron
		{"ENST00000269305.4:c.1799T>A", true, true},            // wholly coding
		{"ENST00000269305.4:c.215_217del", true, true},         // a coding range
		{"ENST00000579081.1:c.*163dup", false, false},          // 3'UTR
		{"ENST00000540690.1:c.-148dup", false, false},          // 5'UTR
		{"ENST00000535955.1:c.43-3169del", false, false},       // deep intronic
		{"ENST00000502371.1:n.1552_1553del", false, false},     // non-coding transcript has no protein
		{"", false, false},
	} {
		start, end := codingEndsOfHGVSc(c.hgvsc)
		if start != c.wantStart || end != c.wantEnd {
			t.Errorf("codingEndsOfHGVSc(%q) = %v, %v; want %v, %v",
				c.hgvsc, start, end, c.wantStart, c.wantEnd)
		}
	}
}
