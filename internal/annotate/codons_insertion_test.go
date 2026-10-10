package annotate

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Inserting CAG before CDS base 4 of ATG CAT GAA TAA 3'-shifts twice, because
// the next two bases are C and A, and the inserted bases rotate CAG -> AGC ->
// GCA as they go. VEP reports the Codons at the position the change actually
// sits at, with the bases in the order the variant gives them.
func TestInsertionCodonsUseUnshiftedAllele(t *testing.T) {
	tr := &cache.Transcript{
		ID: "ENST_INSCODON", GeneName: "INSCODON", Chrom: "1", Strand: 1,
		Start: 1000, End: 1011, CDSStart: 1000, CDSEnd: 1011,
		Biotype:     "protein_coding",
		CDSSequence: "ATGCATGAATAA",
		Exons: []cache.Exon{
			{Number: 1, Start: 1000, End: 1011, CDSStart: 1000, CDSEnd: 1011, Frame: 0},
		},
	}
	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "", Alt: "CAG"}
	result := PredictConsequence(v, tr)
	result.HGVSc = FormatHGVSc(v, tr, result)

	if result.HGVSOffset != 2 {
		t.Fatalf("expected the insertion to 3'-shift 2 bases, got offset %d (HGVSc %q)",
			result.HGVSOffset, result.HGVSc)
	}
	// The insertion sits on a codon boundary, so there are no reference codons
	// to show and VEP writes "-" for that side.
	if want := "-/CAG"; result.CodonChange != want {
		t.Errorf("Codons = %q, want %q (GCA would be the bases after the shift)",
			result.CodonChange, want)
	}
}

// Codons strings taken from a VEP111 MSK-IMPACT MAF, with the consequence VEP
// reports alongside them.
func TestRelabelInsertionSplittingACodon(t *testing.T) {
	tests := []struct {
		name   string
		codons string
		want   string
	}{
		// Three bases one into the codon: the aspartate is replaced by a
		// glycine and a valine, so the protein change is a replacement.
		{"EGFR c.2308_2309insGTT", "gac/gGTTac", ConsequenceProteinAltering},
		{"FOXA1 c.767_768insGTT", "aac/aaGTTc", ConsequenceProteinAltering},
		{"VHL c.295_296insGTA", "cca/cGTAca", ConsequenceProteinAltering},
		// Whole codons added at a boundary.
		{"insertion at a codon boundary", "-/CAG", ConsequenceInframeInsertion},
		// Splits a codon, but the original residue survives at the end:
		// CCT is P, and cCAGct spells Q then P.
		{"original residue survives", "cct/cCAGct", ConsequenceInframeInsertion},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			r := &ConsequenceResult{Consequence: ConsequenceInframeInsertion, CodonChange: tt.codons}
			relabelInsertionSplittingACodon(r)
			if r.Consequence != tt.want {
				t.Errorf("Codons %q -> %q, want %q", tt.codons, r.Consequence, tt.want)
			}
		})
	}
}
