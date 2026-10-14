package annotate

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

func TestHGVScInversion(t *testing.T) {
	// Forward-strand transcript, CDS bases 1..12 at genomic 1000..1011.
	fwd := &cache.Transcript{
		ID: "ENST_INV", GeneName: "INV", Chrom: "1", Strand: 1,
		Start: 1000, End: 1011, CDSStart: 1000, CDSEnd: 1011,
		Biotype:     "protein_coding",
		CDSSequence: "ATGCATGAATAA",
		Exons: []cache.Exon{
			{Number: 1, Start: 1000, End: 1011, CDSStart: 1000, CDSEnd: 1011, Frame: 0},
		},
	}
	tests := []struct {
		name     string
		pos      int64
		ref, alt string
		want     string
	}{
		// CA -> TG is the reverse complement, so HGVS writes it as inv.
		{"two-base inversion", 1003, "CA", "TG", "c.4_5inv"},
		{"six-base inversion", 1003, "CATGAA", "TTCATG", "c.4_9inv"},
		// Same length but not the reverse complement.
		{"delins stays delins", 1003, "CA", "GG", "c.4_5delinsGG"},
		// The reverse complement of a single base is not an inversion; HGVS
		// has the substitution operator for that.
		{"single base is a substitution", 1003, "C", "G", "c.4C>G"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := &vcf.Variant{Chrom: "1", Pos: tt.pos, Ref: tt.ref, Alt: tt.alt}
			result := PredictConsequence(v, fwd)
			if got := FormatHGVSc(v, fwd, result); got != tt.want {
				t.Errorf("got %q, want %q", got, tt.want)
			}
		})
	}
}

func TestHGVScInversionReverseStrand(t *testing.T) {
	// Reverse-strand transcript: reverse-complementing both alleles preserves
	// the inversion relation, so the operator is the same.
	rev := &cache.Transcript{
		ID: "ENST_INVREV", GeneName: "INVREV", Chrom: "1", Strand: -1,
		Start: 1000, End: 1011, CDSStart: 1000, CDSEnd: 1011,
		Biotype:     "protein_coding",
		CDSSequence: "ATGCATGAATAA",
		Exons: []cache.Exon{
			{Number: 1, Start: 1000, End: 1011, CDSStart: 1000, CDSEnd: 1011, Frame: 0},
		},
	}
	v := &vcf.Variant{Chrom: "1", Pos: 1003, Ref: "CA", Alt: "TG"}
	result := PredictConsequence(v, rev)
	if got, want := FormatHGVSc(v, rev, result), "c.8_9inv"; got != want {
		t.Errorf("got %q, want %q", got, want)
	}
}
