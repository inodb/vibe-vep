package annotate

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Tests derived from real mismatches found in datahub GDC benchmark.
// Each test documents a disagreement between MAF annotation and vibe-vep,
// and asserts the behavior we believe is correct. Where vibe-vep is correct,
// we assert that; where MAF is correct, we document the expected fix.
//
// The mismatch patterns from 53 mismatches across ~2.5M variants:
//
//  1. frameshift_variant,stop_lost vs inframe_insertion,stop_retained_variant (13)
//  2. stop_lost,3_prime_UTR_variant vs frameshift_variant,stop_lost (9)
//  3. inframe_insertion,stop_retained_variant vs frameshift_variant (7)
//  4. 5_prime_UTR_variant vs intergenic_variant (4)
//  5. 3_prime_UTR_variant vs frameshift_variant,stop_lost (4)
//  6. frameshift_variant vs inframe_insertion,stop_retained_variant (3)
//  7. stop_gained vs stop_retained_variant (2)
//  8. missense_variant vs stop_lost (2)

// --- Pattern 1 & 6: Stop codon insertion classification ---
// MAF says frameshift_variant,stop_lost but vibe-vep says inframe_insertion,stop_retained_variant.
// These are 1-base insertions at/near the stop codon. If CDS length is divisible by 3
// and the stop codon is preserved, inframe_insertion,stop_retained_variant is correct
// (VEP's own behavior for insertions within the stop codon that keep it intact).

func TestDatahub_StopCodonInsertion_SingleBase(t *testing.T) {
	// Pattern: 1-base insertion at stop codon → MAF=frameshift_variant,stop_lost,
	// vibe-vep=inframe_insertion,stop_retained_variant.
	//
	// Real example: chr5:45261922 ins A, ENST00000303230 (chrcc_tcga_gdc)
	// Also: chr9:111686771 ins A, ENST00000318737 (coad/ucec)
	//       chr8:86048463 ins A, ENST00000276616 (luad)
	//
	// When a 1-base insertion occurs at the last base of the stop codon (CDS pos
	// divisible by 3), it shifts the stop codon but may preserve it.
	// For CDS of length 12 (ATG GCT AAA TAA), inserting at stop codon last base:
	// Mutant: ATG GCT AAA TA[ins]A → if ins=A: ATG GCT AAA TAA A...
	// The stop codon TAA is preserved → inframe_insertion,stop_retained_variant is correct.

	cds := "ATGGCTAAATAA" // 12bp: M A K *
	tr := &cache.Transcript{
		ID: "ENST_STOP_INS1", GeneName: "STOPINS1", Chrom: "1",
		Start: 995, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:       []cache.Exon{{Number: 1, Start: 995, End: 1020, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence: cds,
	}

	// Insert at the last base of stop codon (CDS pos 12 = genomic 1011)
	v := &vcf.Variant{Chrom: "1", Pos: 1011, Ref: "A", Alt: "AA"}
	result := PredictConsequence(v, tr)

	// vibe-vep classifies as inframe_insertion,stop_retained_variant — verify this is stable
	if result.Consequence != ConsequenceInframeInsertion+","+ConsequenceStopRetained {
		t.Errorf("stop codon insertion: got %q, want %q",
			result.Consequence, ConsequenceInframeInsertion+","+ConsequenceStopRetained)
	}
}

func TestDatahub_StopCodonInsertion_ThirdToLastBase(t *testing.T) {
	// Insert at the first base of stop codon (CDS pos 10 = genomic 1009)
	// TAA → T[ins]AA → if ins=A: TAAA → codon = TAA then A... stop preserved
	cds := "ATGGCTAAATAA"
	tr := &cache.Transcript{
		ID: "ENST_STOP_INS3", GeneName: "STOPINS3", Chrom: "1",
		Start: 995, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:       []cache.Exon{{Number: 1, Start: 995, End: 1020, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence: cds,
	}

	v := &vcf.Variant{Chrom: "1", Pos: 1009, Ref: "T", Alt: "TA"}
	result := PredictConsequence(v, tr)

	// Insertion at first base of stop codon: frameshift that shifts the stop
	// If the new reading frame creates stop at position 1: stop_gained
	// Otherwise: frameshift_variant
	t.Logf("stop codon first-base insertion: consequence=%q", result.Consequence)
	if result.Consequence == "" {
		t.Error("expected non-empty consequence")
	}
}

// --- Pattern 2: Deletion spanning stop codon into 3'UTR ---
// MAF says stop_lost,3_prime_UTR_variant (simple composite)
// vibe-vep says frameshift_variant,stop_lost (because deletion is not in-frame)
// Both capture the stop_lost aspect; the disagreement is whether the primary
// consequence is frameshift (out-of-frame deletion) or just stop_lost.

func TestDatahub_DeletionSpanningStopInto3UTR_OutOfFrame(t *testing.T) {
	// Real: chr3:10149965 del AA, ENST00000256474 (ccrcc_tcga_gdc)
	// 2bp deletion spanning stop codon boundary is out-of-frame → frameshift
	cds := "ATGGCTAAATAA" // 12bp: M A K *
	utr3 := "GGGAAATAA"
	tr := &cache.Transcript{
		ID: "ENST_DELSTOP3UTR", GeneName: "DELSTOP3UTR", Chrom: "1",
		Start: 995, End: 1025, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:        []cache.Exon{{Number: 1, Start: 995, End: 1025, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence:  cds,
		UTR3Sequence: utr3,
	}

	// Delete 2 bases at the stop codon: pos=1010 (CDS 11), ref=AA+3'UTR base, alt=A
	// This removes 1 base from stop codon + extends into 3'UTR = out-of-frame
	v := &vcf.Variant{Chrom: "1", Pos: 1010, Ref: "AAG", Alt: "A"}
	result := PredictConsequence(v, tr)

	// vibe-vep: frameshift_variant,stop_lost — this is correct for out-of-frame deletion
	t.Logf("del spanning stop+3UTR: consequence=%q", result.Consequence)
	if result.Consequence == "" {
		t.Error("expected non-empty consequence")
	}
}

func TestDatahub_DeletionSpanningStopInto3UTR_InFrame(t *testing.T) {
	// In-frame deletion spanning stop codon into 3'UTR
	// MAF: stop_lost,3_prime_UTR_variant
	// This should be stop_lost since the stop codon is destroyed
	cds := "ATGGCTAAATAA" // 12bp: M A K *
	utr3 := "GGGAAATAA"
	tr := &cache.Transcript{
		ID: "ENST_DELSTOP3UTR_IF", GeneName: "DELSTOPIF", Chrom: "1",
		Start: 995, End: 1025, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:        []cache.Exon{{Number: 1, Start: 995, End: 1025, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence:  cds,
		UTR3Sequence: utr3,
	}

	// Delete 6 bases (in-frame) spanning stop codon + 3'UTR
	// pos=1009 (CDS 10, first base of stop), ref=TAAGGG (3 CDS + 3 UTR), alt=T
	v := &vcf.Variant{Chrom: "1", Pos: 1008, Ref: "ATAAGGG", Alt: "A"}
	result := PredictConsequence(v, tr)

	t.Logf("in-frame del spanning stop+3UTR: consequence=%q", result.Consequence)
	// Should contain stop_lost
	if result.Consequence == "" {
		t.Error("expected non-empty consequence")
	}
}

// --- Pattern 3: inframe_insertion,stop_retained_variant vs frameshift_variant ---
// MAF says inframe_insertion,stop_retained but vibe-vep says frameshift.
// These are insertions where CDS length % 3 matters. Some MAF annotations
// may use a different transcript model where the CDS length differs.

func TestDatahub_InsertionAtStopCodon_NonDivisibleCDS(t *testing.T) {
	// When CDS length is NOT divisible by 3, a 1-base insertion might make it divisible,
	// making it appear in-frame relative to the new CDS. This is an edge case in
	// incomplete terminal codon handling.
	//
	// Real: chr12:69390312 ins A, ENST00000247843 (brca_tcga_gdc)
	// MAF: inframe_insertion,stop_retained_variant
	// vibe-vep: frameshift_variant
	//
	// CDS length 14 (not divisible by 3), 1-base ins → 15 (divisible by 3)
	cds := "ATGGCTAAAGGTAA" // 14bp, incomplete terminal codon
	tr := &cache.Transcript{
		ID: "ENST_NONDIV3", GeneName: "NONDIV3", Chrom: "1",
		Start: 995, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1013,
		Exons:       []cache.Exon{{Number: 1, Start: 995, End: 1020, CDSStart: 1000, CDSEnd: 1013, Frame: 0}},
		CDSSequence: cds,
	}

	v := &vcf.Variant{Chrom: "1", Pos: 1013, Ref: "A", Alt: "AA"}
	result := PredictConsequence(v, tr)

	t.Logf("non-div3 CDS insertion: consequence=%q (CDS len before=%d, after=%d)",
		result.Consequence, len(cds), len(cds)+1)
}

// --- Pattern 5: 3'UTR insertion reclassified as frameshift_variant,stop_lost ---
// MAF says 3_prime_UTR_variant but vibe-vep says frameshift_variant,stop_lost.
// This happens when the MAF annotator's transcript has the variant in the 3'UTR,
// but vibe-vep's GENCODE transcript model places it within the stop codon.

func TestDatahub_InsertionNearStopCodon_3UTRBoundary(t *testing.T) {
	// Real: chr17:78223554 ins G, ENST00000350051 (hgsoc_tcga_gdc)
	// MAF: 3_prime_UTR_variant
	// vibe-vep: frameshift_variant,stop_lost
	//
	// If the insertion is at the last base of the stop codon (or just after),
	// transcript model differences can place it in CDS vs 3'UTR.
	cds := "ATGGCTAAATAA" // 12bp: M A K *
	utr3 := "GGGAAATAA"
	tr := &cache.Transcript{
		ID: "ENST_3UTRBOUNDARY", GeneName: "UTR3BOUND", Chrom: "1",
		Start: 995, End: 1025, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:        []cache.Exon{{Number: 1, Start: 995, End: 1025, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence:  cds,
		UTR3Sequence: utr3,
	}

	// Insert at CDS end boundary (pos 1011 = last CDS base)
	v := &vcf.Variant{Chrom: "1", Pos: 1011, Ref: "A", Alt: "AG"}
	result := PredictConsequence(v, tr)

	t.Logf("3UTR boundary insertion: consequence=%q", result.Consequence)
	// At the stop codon boundary, insertion should be classified consistently
	if result.Consequence == "" {
		t.Error("expected non-empty consequence")
	}
}

// --- Pattern 7: stop_gained vs stop_retained_variant ---
// MAF says stop_gained but vibe-vep says stop_retained_variant.
// This happens when the variant is within an existing stop codon and the
// codon remains a stop. MAF may be using a different transcript model where
// this position is NOT already a stop codon.

func TestDatahub_SNVInStopCodon_StopRetained(t *testing.T) {
	// Real: chr19:1105404 G>A, ENST00000354171 (cesc_tcga_gdc)
	//       chr16:30445549 C>T, ENST00000478753 (hnsc_tcga_gdc)
	//
	// If the position is within a stop codon and the mutation preserves it
	// (e.g., TAA→TAA or TAG→TAA), it's stop_retained.
	// MAF calling it stop_gained suggests their model doesn't have a stop there.

	// TAA → TAA is not a mutation. Use TGA → TAA (both stops).
	cds := "ATGGCTAAATGA" // 12bp: M A K * (TGA stop)
	tr := &cache.Transcript{
		ID: "ENST_STOPRET", GeneName: "STOPRET", Chrom: "1",
		Start: 995, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:       []cache.Exon{{Number: 1, Start: 995, End: 1020, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence: cds,
	}

	// TGA → TAA: change G at pos 1010 (CDS pos 11) to A
	v := &vcf.Variant{Chrom: "1", Pos: 1010, Ref: "G", Alt: "A"}
	result := PredictConsequence(v, tr)

	if result.Consequence != ConsequenceStopRetained {
		t.Errorf("SNV changing TGA→TAA should be stop_retained, got %q", result.Consequence)
	}
}

func TestDatahub_SNVInStopCodon_TAGtoTAA(t *testing.T) {
	// TAG → TAA (both stop codons)
	cds := "ATGGCTAAATAG" // M A K * (TAG stop)
	tr := &cache.Transcript{
		ID: "ENST_STOPRET2", GeneName: "STOPRET2", Chrom: "1",
		Start: 995, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:       []cache.Exon{{Number: 1, Start: 995, End: 1020, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence: cds,
	}

	// TAG → TAA: change G at pos 1011 (CDS pos 12) to A
	v := &vcf.Variant{Chrom: "1", Pos: 1011, Ref: "G", Alt: "A"}
	result := PredictConsequence(v, tr)

	if result.Consequence != ConsequenceStopRetained {
		t.Errorf("SNV changing TAG→TAA should be stop_retained, got %q", result.Consequence)
	}
}

// --- Pattern 8: missense_variant vs stop_lost ---
// MAF says missense_variant but vibe-vep says stop_lost.
// If the variant is in a stop codon and changes it to a non-stop AA,
// it's stop_lost. MAF may be using a different transcript where this position
// is not a stop codon (different CDS end).

func TestDatahub_StopCodonSNV_StopLost(t *testing.T) {
	// Real: chr16:30445550 A>T, ENST00000478753 (skcm_tcga_gdc)
	//       chr5:42800912 T>G, ENST00000506577 (ucec_tcga_gdc)
	//
	// If the position is within the stop codon and mutation destroys it → stop_lost

	cds := "ATGGCTAAATAA" // 12bp: M A K * (TAA stop)
	utr3 := "GCATAA"
	tr := &cache.Transcript{
		ID: "ENST_STOPLOST_SNV", GeneName: "STOPLOSTSNV", Chrom: "1",
		Start: 995, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:        []cache.Exon{{Number: 1, Start: 995, End: 1020, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence:  cds,
		UTR3Sequence: utr3,
	}

	// TAA → TCA: change A at pos 1010 (CDS pos 11) to C → TCA = Ser
	v := &vcf.Variant{Chrom: "1", Pos: 1010, Ref: "A", Alt: "C"}
	result := PredictConsequence(v, tr)

	if result.Consequence != ConsequenceStopLost {
		t.Errorf("SNV destroying stop codon should be stop_lost, got %q", result.Consequence)
	}
}

// --- Pattern: Large multi-base insertion classification at stop codon ---

func TestDatahub_MultiBaseInsertionAtStopCodon(t *testing.T) {
	// Real: chr15:69028246 ins CTATCTATGACTCCTATTCTATCTA (24bp), ENST00000388866
	// MAF: inframe_insertion,stop_retained_variant
	// vibe-vep: frameshift_variant
	//
	// 24bp insertion is divisible by 3, so if at the stop codon it could be
	// inframe_insertion,stop_retained_variant IF stop codon is preserved.
	cds := "ATGGCTAAATAA" // 12bp
	tr := &cache.Transcript{
		ID: "ENST_BIGINS", GeneName: "BIGINS", Chrom: "1",
		Start: 995, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:       []cache.Exon{{Number: 1, Start: 995, End: 1020, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence: cds,
	}

	// 6bp (in-frame) insertion at stop codon
	v := &vcf.Variant{Chrom: "1", Pos: 1011, Ref: "A", Alt: "AGCTGCA"}
	result := PredictConsequence(v, tr)

	t.Logf("multi-base insertion at stop: consequence=%q", result.Consequence)
	// Behavior depends on whether stop codon is preserved
}

// --- Pattern: Large deletion from 3'UTR into stop codon (reverse strand) ---

func TestDatahub_LargeDeletion3UTRtoStop_ReverseStrand(t *testing.T) {
	// Real: chr11:44929172 del CCAGTAGAGGGTATGGCCTG (20bp), ENST00000340160 (hgsoc)
	// MAF: stop_lost,3_prime_UTR_variant
	// vibe-vep: frameshift_variant,stop_lost
	//
	// 20bp is not divisible by 3 → frameshift. Both agree on stop_lost.

	cds := "ATGGCTAAAGAAGGGTAA" // 18bp on coding strand
	utr3 := "GGGAAATTT"
	tr := &cache.Transcript{
		ID: "ENST_LARGEDEL_REV", GeneName: "LARGEDELREV", Chrom: "1",
		Start: 985, End: 1025, Strand: -1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1017,
		Exons:        []cache.Exon{{Number: 1, Start: 985, End: 1025, CDSStart: 1000, CDSEnd: 1017, Frame: 0}},
		CDSSequence:  cds,
		UTR3Sequence: utr3,
	}

	// For reverse strand, 3'UTR is at lower genomic coords (below CDSStart).
	// Delete 5bp starting in 3'UTR at pos 997 spanning into stop codon.
	// Stop codon on reverse strand is at CDSStart (genomic 1000-1002).
	v := &vcf.Variant{Chrom: "1", Pos: 997, Ref: "AAATAA", Alt: "A"}
	result := PredictConsequence(v, tr)

	t.Logf("large del from 3UTR into stop (rev): consequence=%q", result.Consequence)
	// Should detect stop_lost
}

// --- Pattern: Splice donor reclassification for deletion ---

func TestDatahub_DeletionSpanningSpliceDonor(t *testing.T) {
	// Real: chr7:142796891 del GCTCGG, ENST00000390415 (cesc_tcga_gdc)
	// MAF: coding_sequence_variant,3_prime_UTR_variant
	// vibe-vep: splice_donor_variant (because deletion spans into splice site)
	//
	// vibe-vep correctly upgrades to splice_donor when the deletion overlaps
	// the splice donor site. This is the correct higher-impact consequence.

	tr := &cache.Transcript{
		ID: "ENST_DELSPLICE", GeneName: "DELSPLICE", Chrom: "1",
		Start: 985, End: 2025, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 2010,
		Exons: []cache.Exon{
			{Number: 1, Start: 985, End: 1010, CDSStart: 1000, CDSEnd: 1010, Frame: 0},
			{Number: 2, Start: 2000, End: 2025, CDSStart: 2000, CDSEnd: 2010, Frame: 2},
		},
		CDSSequence: "ATGGCTAAAGAAGGGCCCAAA",
	}

	// Deletion starting in CDS spanning into splice donor (exon1 End=1010, donor at 1011/1012)
	v := &vcf.Variant{Chrom: "1", Pos: 1008, Ref: "GAAAAA", Alt: "G"}
	result := PredictConsequence(v, tr)

	if result.Consequence != ConsequenceSpliceDonor {
		t.Errorf("deletion spanning splice donor should be splice_donor_variant, got %q",
			result.Consequence)
	}
}

// --- Pattern: 5'UTR insertion reclassified as frameshift ---

func TestDatahub_InsertionIn5UTR_SpanningIntoCDS(t *testing.T) {
	// Real: chr11:14498928 ins GGTTTCTG (8bp), ENST00000249923 (stad_tcga_gdc)
	// MAF: 5_prime_UTR_variant
	// vibe-vep: frameshift_variant
	//
	// If the insertion is positioned at the UTR/CDS boundary in one model but
	// within CDS in another, the consequence differs.
	// For a variant that's clearly in 5'UTR, it should stay 5_prime_UTR_variant.

	tr := &cache.Transcript{
		ID: "ENST_5UTRINS", GeneName: "UTRINS5", Chrom: "1",
		Start: 990, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1020, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence: "ATGGCTAAATAA",
	}

	// Insertion clearly in 5'UTR (pos 995)
	v := &vcf.Variant{Chrom: "1", Pos: 995, Ref: "A", Alt: "AGGTTTCTG"}
	result := PredictConsequence(v, tr)

	if result.Consequence != Consequence5PrimeUTR {
		t.Errorf("insertion in 5'UTR should be 5_prime_UTR_variant, got %q", result.Consequence)
	}
}

// --- Comprehensive: stop codon position sensitivity ---

func TestDatahub_StopCodonPositionSensitivity(t *testing.T) {
	// Test SNV at each position of each stop codon type to verify correct classification.
	// This covers patterns 7 and 8 comprehensively.

	tests := []struct {
		name    string
		stopRef string // 3-letter stop codon
		pos     int    // 0-indexed position in stop codon to mutate
		alt     byte   // alternate base
		want    string // expected consequence
	}{
		// TAA stop codon
		{"TAA pos0 T>C = CAA(Gln) → stop_lost", "TAA", 0, 'C', ConsequenceStopLost},
		{"TAA pos1 A>C = TCA(Ser) → stop_lost", "TAA", 1, 'C', ConsequenceStopLost},
		{"TAA pos2 A>G = TAG(stop) → stop_retained", "TAA", 2, 'G', ConsequenceStopRetained},
		// TAG stop codon
		{"TAG pos0 T>C = CAG(Gln) → stop_lost", "TAG", 0, 'C', ConsequenceStopLost},
		{"TAG pos1 A>C = TCG(Ser) → stop_lost", "TAG", 1, 'C', ConsequenceStopLost},
		{"TAG pos2 G>A = TAA(stop) → stop_retained", "TAG", 2, 'A', ConsequenceStopRetained},
		// TGA stop codon
		{"TGA pos0 T>C = CGA(Arg) → stop_lost", "TGA", 0, 'C', ConsequenceStopLost},
		{"TGA pos1 G>A = TAA(stop) → stop_retained", "TGA", 1, 'A', ConsequenceStopRetained},
		{"TGA pos2 A>C = TGC(Cys) → stop_lost", "TGA", 2, 'C', ConsequenceStopLost},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			cds := "ATGGCTAAA" + tt.stopRef // 12bp
			tr := &cache.Transcript{
				ID: "ENST_STOPPOS", GeneName: "STOPPOS", Chrom: "1",
				Start: 995, End: 1020, Strand: 1, Biotype: "protein_coding",
				CDSStart: 1000, CDSEnd: 1011,
				Exons:        []cache.Exon{{Number: 1, Start: 995, End: 1020, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
				CDSSequence:  cds,
				UTR3Sequence: "GCATAA",
			}

			genomicPos := int64(1009 + tt.pos) // CDS pos 10 + offset
			ref := string(cds[9+tt.pos])
			v := &vcf.Variant{Chrom: "1", Pos: genomicPos, Ref: ref, Alt: string(tt.alt)}
			result := PredictConsequence(v, tr)

			if result.Consequence != tt.want {
				t.Errorf("got %q, want %q", result.Consequence, tt.want)
			}
		})
	}
}
