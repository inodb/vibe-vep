package output

import "testing"

// VEP writes amino_acids as the translation of the codons on either side of
// the change. These are real values from vep111.genomenexus.org.
func TestAminoAcidsFromCodons(t *testing.T) {
	for _, c := range []struct{ codons, want string }{
		{"aaA/aa", "K/X"},                      // frameshift deletion
		{"-/A", "-/X"},                         // frameshift insertion into nothing
		{"cct/ccCt", "P/PX"},                   // frameshift insertion keeping one codon
		{"ACC/-", "T/-"},                       // in-frame deletion of one codon
		{"gTg/gAg", "V/E"},                     // plain substitution
		{"cGa/ca", "R/X"},                      // frameshift on the reverse strand
		{"aAGGATCGGGCCCACAtc/atc", "KDRAHI/I"}, // in-frame deletion of five
		{"", ""},                               // nothing to translate
		{"-/-", ""},                            // neither side has codons
	} {
		if got := aminoAcidsFromCodons(c.codons); got != c.want {
			t.Errorf("aminoAcidsFromCodons(%q) = %q, want %q", c.codons, got, c.want)
		}
	}
}

// VEP names a synonymous change once: the codon changed, the amino acid did
// not. genome-nexus reads that single letter back as both sides.
func TestSynonymousAminoAcidIsNamedOnce(t *testing.T) {
	if got := collapseUnchangedAminoAcid("T/T"); got != "T" {
		t.Errorf("collapseUnchangedAminoAcid(\"T/T\") = %q, want \"T\"", got)
	}
	if got := collapseUnchangedAminoAcid("V/E"); got != "V/E" {
		t.Errorf("collapseUnchangedAminoAcid(\"V/E\") = %q, want \"V/E\"", got)
	}
	if ref, alt := splitAminoAcids("T"); ref != "T" || alt != "T" {
		t.Errorf("splitAminoAcids(\"T\") = %q, %q, want \"T\", \"T\"", ref, alt)
	}
	if ref, alt := splitAminoAcids("V/E"); ref != "V" || alt != "E" {
		t.Errorf("splitAminoAcids(\"V/E\") = %q, %q, want \"V\", \"E\"", ref, alt)
	}
}
