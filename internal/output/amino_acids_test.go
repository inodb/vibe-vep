package output

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
)

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
		// A change straddling a codon boundary reports both codons, even
		// where only one of the residues moves (HNF1A 12:121426810 CC>TT).
		{"gtCCgc/gtTTgc", "VR/VC"},
		{"gACAAa/ga", "DK/X"},            // two codons in, a part-codon out
		{"gaGGACCACAGg/gag", "EDHR/E"},   // APC c.1951_1958+1del
		{"GATGATGATGAAGAG/-", "DDDEE/-"}, // NPM1 c.511_524+1del
		{"tac/tAac", "Y/*"},              // a stop and nothing else drops the X
		{"tac/tGACTac", "Y/*LX"},         // but a stop with a residue after it keeps it
		{"aac/aTATGac", "N/I*X"},         // and so does one with a residue before
		{"", ""},                         // nothing to translate
		{"-/-", ""},                      // neither side has codons
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

// The codons cover every codon the change touches; the named protein change
// knows only the residue it renames. Where both are present the codons win,
// which is what makes a two-base change across a boundary read VR/VC.
func TestAminoAcidsPreferTheCodons(t *testing.T) {
	ann := &annotate.Annotation{
		CodonChange:     "gtCCgc/gtTTgc",
		AminoAcidChange: "R168C",
	}
	if got := vepAminoAcids(ann); got != "VR/VC" {
		t.Errorf("vepAminoAcids = %q, want \"VR/VC\"", got)
	}

	// With no codons to read -- an intronic or non-coding change -- the
	// named change is all there is.
	ann = &annotate.Annotation{AminoAcidChange: "R168C"}
	if got := vepAminoAcids(ann); got != "R/C" {
		t.Errorf("vepAminoAcids = %q, want \"R/C\"", got)
	}
}
