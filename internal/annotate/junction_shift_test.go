package annotate

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// TestExonJunction_DeletionShift verifies 3' shift behavior at exon-exon junctions.
// HGVS spec (deletion.md): exception for exon/exon junctions when identical nucleotides
// flank the junction — the deletion should be assigned to the FIRST exon, not shifted across.
func TestExonJunction_DeletionShift(t *testing.T) {
	// 2-exon transcript, forward strand
	// Exon 1 CDS: pos 1-21 (genomic 1000-1020)
	// Exon 2 CDS: pos 22-102 (genomic 1100-1180)
	//
	// CDS constructed so pos 21 (index 20) = 'T' and pos 22 (index 21) = 'T'
	// Same base at junction → 3' shift could incorrectly cross it
	cds := "ATGATGATGATGATGATGATTTGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
	// Verify: CDS[20]='T', CDS[21]='T'
	if cds[20] != 'T' || cds[21] != 'T' {
		t.Fatalf("CDS setup error: CDS[20]=%c CDS[21]=%c, both should be T", cds[20], cds[21])
	}

	transcript := &cache.Transcript{
		ID: "ENST_JUNC", GeneName: "JUNC", Chrom: "1",
		Start: 990, End: 1300, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1180,
		Exons: []cache.Exon{
			{Number: 1, Start: 990, End: 1020, CDSStart: 1000, CDSEnd: 1020, Frame: 0},
			{Number: 2, Start: 1100, End: 1300, CDSStart: 1100, CDSEnd: 1180, Frame: 0},
		},
		CDSSequence: cds,
	}
	transcript.BuildCDSIndex()

	// Delete T at last base of exon 1 (genomic 1020, CDS pos 21)
	// VCF: pos=1019, ref=AT, alt=A → deletes genomic 1020
	v := &vcf.Variant{Chrom: "1", Pos: 1019, Ref: "AT", Alt: "A"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("Got: %s", hgvsc)
	// Per HGVS spec, deletion at exon boundary should NOT shift across the junction.
	// The deletion should stay at c.21del, not become c.22del.
	if hgvsc == "c.22del" {
		t.Errorf("BUG: 3' shift crossed exon-exon junction: got c.22del, want c.21del")
	}
}

// TestExonJunction_DupShift verifies dup 3' shift respects exon junctions.
func TestExonJunction_DupShift(t *testing.T) {
	// Same setup: T at exon 1 end (CDS 21) and T at exon 2 start (CDS 22)
	cds := "ATGATGATGATGATGATGATTTGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"

	transcript := &cache.Transcript{
		ID: "ENST_JUNC2", GeneName: "JUNC2", Chrom: "1",
		Start: 990, End: 1300, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1180,
		Exons: []cache.Exon{
			{Number: 1, Start: 990, End: 1020, CDSStart: 1000, CDSEnd: 1020, Frame: 0},
			{Number: 2, Start: 1100, End: 1300, CDSStart: 1100, CDSEnd: 1180, Frame: 0},
		},
		CDSSequence: cds,
	}
	transcript.BuildCDSIndex()

	// Insert T after last T of exon 1 → should be c.21dup, not c.22dup
	v := &vcf.Variant{Chrom: "1", Pos: 1020, Ref: "T", Alt: "TT"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("Got: %s", hgvsc)
	if hgvsc == "c.22dup" {
		t.Errorf("BUG: 3' shift crossed exon-exon junction for dup: got c.22dup, want c.21dup")
	}
}

// TestExonJunction_DeletionShift_PolyT tests with a longer run of T's at the junction.
// This is more likely to trigger the shift since there are multiple T's to shift through.
func TestExonJunction_DeletionShift_PolyT(t *testing.T) {
	// Exon 1 CDS ends with TTT (pos 19-21), exon 2 CDS starts with TT (pos 22-23)
	// Deleting one T from exon 1 should shift within exon 1 to pos 21, not cross to exon 2
	cds := "ATGATGATGATGATGATGTTTTTGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
	// pos 19-21 (0-based 18-20) = TTT, pos 22-23 (0-based 21-22) = TT
	if cds[18] != 'T' || cds[19] != 'T' || cds[20] != 'T' || cds[21] != 'T' || cds[22] != 'T' {
		t.Fatalf("CDS setup error")
	}

	transcript := &cache.Transcript{
		ID: "ENST_JUNC3", GeneName: "JUNC3", Chrom: "1",
		Start: 990, End: 1300, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1180,
		Exons: []cache.Exon{
			{Number: 1, Start: 990, End: 1020, CDSStart: 1000, CDSEnd: 1020, Frame: 0},
			{Number: 2, Start: 1100, End: 1300, CDSStart: 1100, CDSEnd: 1180, Frame: 0},
		},
		CDSSequence: cds,
	}
	transcript.BuildCDSIndex()

	// Delete T at CDS pos 19 (first T in the run within exon 1)
	// VCF: genomic 1018, ref=GT, alt=G → deletes genomic 1019 (CDS 20)
	// Without junction awareness: shifts from 20 to 22 (crossing junction to exon 2)
	// With junction awareness: shifts from 20 to 21 (last T in exon 1)
	v := &vcf.Variant{Chrom: "1", Pos: 1018, Ref: "GT", Alt: "G"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("Got: %s", hgvsc)
	// Should be c.21del (shifted within exon 1), not c.22del or c.23del (crossed to exon 2)
	if hgvsc != "c.21del" {
		t.Errorf("Expected c.21del (3' shift within exon 1), got %s", hgvsc)
	}
}
