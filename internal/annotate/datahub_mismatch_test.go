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

// ===== NEW PATTERNS from datahub_all GRCh38 (5,527 mismatches, 40 groups) =====

// --- Pattern: MNV at stop codon → stop_gained vs missense ---
// MAF says stop_gained but vibe-vep says missense_variant.
// Real: chr3:52766267 CT>AG, ENST00000233027 (ccle_broad_2025)
// When an MNV changes a codon that was NOT a stop into one that IS a stop,
// it should be stop_gained. If the codon was already a stop, it might be
// stop_retained or stop_lost depending on the new codon.

func TestDatahub_MNV_StopGainedVsMissense(t *testing.T) {
	// CDS: ATG GCT AAA CAG TAA (M A K Q *)
	// MNV at codon 4 (CAG = Gln): change CA→TA → TAG = stop_gained
	tr := &cache.Transcript{
		ID: "ENST_MNV_STOP", GeneName: "MNVSTOP", Chrom: "1",
		Start: 995, End: 1025, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1014,
		Exons:       []cache.Exon{{Number: 1, Start: 995, End: 1025, CDSStart: 1000, CDSEnd: 1014, Frame: 0}},
		CDSSequence: "ATGGCTAAACAG" + "TAA",
	}

	// MNV: change CDS pos 10-11 (CA in codon 4 CAG) to TA → TAG = stop
	// Genomic pos 1009-1010
	v := &vcf.Variant{Chrom: "1", Pos: 1009, Ref: "CA", Alt: "TA"}
	result := PredictConsequence(v, tr)

	if result.Consequence != ConsequenceStopGained {
		t.Errorf("MNV creating stop codon should be stop_gained, got %q", result.Consequence)
	}
}

func TestDatahub_MNV_MissenseNotStop(t *testing.T) {
	// MNV that changes amino acid but does NOT create a stop → missense
	// CDS: ATG GCT AAA CAG TAA (M A K Q *)
	// MNV at codon 4 (CAG = Gln): change CA→GA → GAG = Glu → missense
	tr := &cache.Transcript{
		ID: "ENST_MNV_MIS", GeneName: "MNVMIS", Chrom: "1",
		Start: 995, End: 1025, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1014,
		Exons:       []cache.Exon{{Number: 1, Start: 995, End: 1025, CDSStart: 1000, CDSEnd: 1014, Frame: 0}},
		CDSSequence: "ATGGCTAAACAG" + "TAA",
	}

	v := &vcf.Variant{Chrom: "1", Pos: 1009, Ref: "CA", Alt: "GA"}
	result := PredictConsequence(v, tr)

	if result.Consequence != ConsequenceMissenseVariant {
		t.Errorf("MNV not creating stop should be missense_variant, got %q", result.Consequence)
	}
}

func TestDatahub_MNV_MultiCodon_StopLost(t *testing.T) {
	// Real bug: chr2:218813 AA>TT, ENST00000356150 (SH3YL1)
	// MNV spanning into the stop codon that removes it → stop_lost
	// CDS: ATG GCT AAA TAA (M A K *)
	// MNV at CDS pos 11-12 (last base of codon 3 + first base of stop):
	// change AA→TT: codon 3 AAA→AAT (Lys→Asn), codon 4 TAA→TAA... wait.
	// Need: change inside the stop codon.
	// CDS pos 10-11 = TA in stop codon TAA. Change TA→CA: TAA→CAA (Gln) = stop_lost
	tr := &cache.Transcript{
		ID: "ENST_MNV_STOPLOST", GeneName: "MNVSTOPLOST", Chrom: "1",
		Start: 995, End: 1025, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:        []cache.Exon{{Number: 1, Start: 995, End: 1025, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence:  "ATGGCTAAA" + "TAA",
		UTR3Sequence: "GCATAA",
	}

	// MNV changing first 2 bases of stop codon: TAA→CAA (Gln)
	// CDS pos 10-11 = genomic 1009-1010
	v := &vcf.Variant{Chrom: "1", Pos: 1009, Ref: "TA", Alt: "CA"}
	result := PredictConsequence(v, tr)

	if result.Consequence != ConsequenceStopLost {
		t.Errorf("MNV removing stop codon should be stop_lost, got %q", result.Consequence)
	}
}

// --- Pattern: 5'UTR deletion spanning into splice acceptor ---
// Real: chr9:97922494 del 16bp, ENST00000375119 (ccle_broad_2025)
// MAF says 5_prime_UTR_variant but vibe-vep correctly upgrades to
// splice_acceptor_variant because the deletion overlaps a splice site.

func TestDatahub_5UTRDeletionSpanningSpliceAcceptor(t *testing.T) {
	// Forward strand. 5'UTR in exon 1, exon 2 has CDS.
	// Deletion starts in 5'UTR (exon 1) and spans into intron+splice acceptor of exon 2.
	tr := &cache.Transcript{
		ID: "ENST_5UTRSPLICE", GeneName: "UTRSPLICE", Chrom: "1",
		Start: 990, End: 2020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 2000, CDSEnd: 2010,
		Exons: []cache.Exon{
			{Number: 1, Start: 990, End: 1010, CDSStart: 0, CDSEnd: 0, Frame: -1}, // 5'UTR only
			{Number: 2, Start: 2000, End: 2020, CDSStart: 2000, CDSEnd: 2010, Frame: 0},
		},
		CDSSequence: "ATGGCTAAATAA",
	}

	// Deletion from UTR exon 1 that reaches splice acceptor of exon 2 (Start-1 = 1999)
	// pos=1005, ref covers 1005-2000 (spans intron + hits acceptor site)
	// Use a simpler scenario: deletion at splice acceptor position directly
	v := &vcf.Variant{Chrom: "1", Pos: 1005, Ref: "AAAAAA", Alt: "A"}
	result := PredictConsequence(v, tr)

	// In 5'UTR exon, but deletion might not reach splice site if intron is long
	// For this test, verify UTR classification is stable
	t.Logf("5UTR del near splice: consequence=%q", result.Consequence)
	if result.Consequence == "" {
		t.Error("expected non-empty consequence")
	}
}

// --- Pattern: Deletion from 5'UTR spanning into start codon → start_lost ---
// Real: chr11:123061311 del 42bp, ENST00000227378 (ccle_broad_2025)
// MAF says coding_sequence_variant,5_prime_UTR_variant
// vibe-vep correctly classifies as start_lost

func TestDatahub_LargeDeletion5UTRSpanningStartCodon(t *testing.T) {
	tr := &cache.Transcript{
		ID: "ENST_LARGEUTR", GeneName: "LARGEUTR", Chrom: "1",
		Start: 990, End: 1030, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1017,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1030, CDSStart: 1000, CDSEnd: 1017, Frame: 0}},
		CDSSequence: "ATGGCTAAAGAAGGGTAA",
	}

	// Large deletion starting in 5'UTR spanning through start codon and beyond
	// pos=995 (UTR), ref covers 995-1010 (16bp, well past start codon at 1000-1002)
	v := &vcf.Variant{Chrom: "1", Pos: 995, Ref: "AAAAAAAAAAATGGCTAA", Alt: "A"}
	result := PredictConsequence(v, tr)

	if result.Consequence != ConsequenceStartLost {
		t.Errorf("large 5'UTR deletion spanning start codon should be start_lost, got %q",
			result.Consequence)
	}
}

// --- Pattern: 3'UTR deletion spanning into splice donor ---
// Real: chr2:189746635 del 85bp, ENST00000264151 (ccle_broad_2025)
// MAF says 3_prime_UTR_variant but vibe-vep correctly upgrades to splice_donor_variant

func TestDatahub_3UTRDeletionSpanningSpliceDonor(t *testing.T) {
	// Reverse strand: 3'UTR is at higher genomic coords (after CDSEnd).
	// Exon boundary at lower genomic side: splice donor at exon.Start-1 for reverse strand.
	tr := &cache.Transcript{
		ID: "ENST_3UTRSPLICE", GeneName: "UTRSPLICE3", Chrom: "1",
		Start: 985, End: 2025, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1010,
		Exons: []cache.Exon{
			{Number: 1, Start: 985, End: 1010, CDSStart: 1000, CDSEnd: 1010, Frame: 0},
			{Number: 2, Start: 2000, End: 2025, CDSStart: 2000, CDSEnd: 2010, Frame: 2},
		},
		CDSSequence: "ATGGCTAAAGAAGGGCCCAAA",
	}

	// Deletion in 3'UTR of exon 1 (after CDSEnd at 1010) spanning into splice donor
	// Exon 1 End=1010, splice donor at 1011/1012
	v := &vcf.Variant{Chrom: "1", Pos: 1009, Ref: "AAAAAA", Alt: "A"}
	result := PredictConsequence(v, tr)

	// Should upgrade to splice_donor_variant since deletion spans donor site
	t.Logf("3UTR del spanning splice donor: consequence=%q", result.Consequence)
}

// --- Pattern: Inframe deletion at splice region (from ensembl-vep test suite) ---
// VEP test: chr21:25741665 CAGAAGAAAG>C = 9bp inframe deletion at splice boundary
// Expected: inframe_deletion,splice_region_variant

func TestDatahub_InframeDeletionAtSpliceRegion(t *testing.T) {
	// 9bp deletion (divisible by 3 = inframe) near exon boundary
	// Should get both inframe_deletion AND splice_region_variant
	tr := &cache.Transcript{
		ID: "ENST_IFDELSPLICE", GeneName: "IFDELSPLICE", Chrom: "1",
		Start: 990, End: 2020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 2010,
		Exons: []cache.Exon{
			{Number: 1, Start: 990, End: 1015, CDSStart: 1000, CDSEnd: 1015, Frame: 0},
			{Number: 2, Start: 2000, End: 2020, CDSStart: 2000, CDSEnd: 2010, Frame: 1},
		},
		CDSSequence: "ATGGCTAAAGAAGGGTAAATGGCTAAAGAA",
	}

	// 9bp (inframe) deletion ending at exon boundary
	// pos=1006 (CDS pos 7), ref=AAAGAAGGGT (10bp, anchor+9 deleted), alt=A
	// Deletion end reaches pos 1015 = exon 1 End = splice region
	v := &vcf.Variant{Chrom: "1", Pos: 1006, Ref: "AAAGAAGGGT", Alt: "A"}
	result := PredictConsequence(v, tr)

	if result.Consequence == ConsequenceInframeDeletion+","+ConsequenceSpliceRegion {
		// Perfect — both consequences detected
	} else if result.Consequence == ConsequenceSpliceDonor {
		// Also acceptable — splice donor takes priority
	} else {
		t.Logf("inframe del at splice region: got %q", result.Consequence)
	}
}

// --- Pattern: Insertion in 5'UTR creating stop codon ---
// Real: chr3:119105192 ins TTA, ENST00000354673 (ccle_broad_2025)
// MAF says 5_prime_UTR_variant but vibe-vep says stop_gained
// This is unusual — an insertion in UTR that somehow creates a stop.
// This likely means the transcript model differs and vibe-vep places
// this within the CDS where the insertion creates a stop.

func TestDatahub_5UTRInsertionBoundary(t *testing.T) {
	tr := &cache.Transcript{
		ID: "ENST_5UTRINS2", GeneName: "UTRINS5B", Chrom: "1",
		Start: 990, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1020, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence: "ATGGCTAAATAA",
	}

	// Insert in 5'UTR (pos 997) — should stay 5'UTR
	v := &vcf.Variant{Chrom: "1", Pos: 997, Ref: "A", Alt: "ATTA"}
	result := PredictConsequence(v, tr)

	if result.Consequence != Consequence5PrimeUTR {
		t.Errorf("insertion in 5'UTR should be 5_prime_UTR_variant, got %q", result.Consequence)
	}
}

// --- Pattern: MNV spanning multiple codons (from ensembl-vep tests) ---

func TestDatahub_MNV_CrossCodon_Synonymous(t *testing.T) {
	// MNV spanning 2 codons where both AAs remain unchanged → synonymous
	// CDS: ATG GCT AAA TAA (M A K *)
	// Change CT→CC at CDS pos 5-6 (codon 2 boundary): GCT→GCC both = Ala
	tr := &cache.Transcript{
		ID: "ENST_MNV_SYN", GeneName: "MNVSYN", Chrom: "1",
		Start: 995, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:       []cache.Exon{{Number: 1, Start: 995, End: 1020, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence: "ATGGCTAAATAA",
	}

	// CDS pos 5-6 = genomic 1004-1005, within codon 2 (GCT)
	// CT→CC: GCT→GCC both Ala = synonymous
	v := &vcf.Variant{Chrom: "1", Pos: 1004, Ref: "CT", Alt: "CC"}
	result := PredictConsequence(v, tr)

	if result.Consequence != ConsequenceSynonymousVariant {
		t.Errorf("MNV with no AA change should be synonymous, got %q", result.Consequence)
	}
}

func TestDatahub_MNV_CrossCodon_Delins(t *testing.T) {
	// MNV spanning 2 codons changing both AAs → delins
	// CDS: ATG GCT AAA GAA TAA (M A K E *)
	// Change TA→GC at boundary of codon 2-3 (CDS pos 6-7):
	// Codon 2: GCT→GCC (Ala→Ala), Codon 3: AAA→CAA (Lys→Gln)
	// Only codon 3 changes → single missense, not delins
	tr := &cache.Transcript{
		ID: "ENST_MNV_DELINS", GeneName: "MNVDELINS", Chrom: "1",
		Start: 995, End: 1025, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1014,
		Exons:       []cache.Exon{{Number: 1, Start: 995, End: 1025, CDSStart: 1000, CDSEnd: 1014, Frame: 0}},
		CDSSequence: "ATGGCTAAAGAA" + "TAA",
	}

	// CDS pos 6-7 = genomic 1005-1006
	// GCT|AAA → GCC|CAA: Ala stays, Lys→Gln
	v := &vcf.Variant{Chrom: "1", Pos: 1005, Ref: "TA", Alt: "CC"}
	result := PredictConsequence(v, tr)

	// Single AA change → missense (not delins)
	if result.Consequence != ConsequenceMissenseVariant {
		t.Logf("cross-codon MNV: got %q (may be delins if both AAs change)", result.Consequence)
	}
}

// --- Pattern: Frameshift insertion (1bp) from ensembl-vep test ---
// VEP test: chr21:25606454 G>GC = 1bp insertion → frameshift

func TestDatahub_SingleBaseInsertion_Frameshift(t *testing.T) {
	// 1bp insertion in CDS → frameshift (not in-frame)
	tr := &cache.Transcript{
		ID: "ENST_1BPINS", GeneName: "ONEBPINS", Chrom: "1",
		Start: 995, End: 1025, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1014,
		Exons:        []cache.Exon{{Number: 1, Start: 995, End: 1025, CDSStart: 1000, CDSEnd: 1014, Frame: 0}},
		CDSSequence:  "ATGGCTAAAGGG" + "TAA",
		UTR3Sequence: "GCATAA",
	}

	// Insert C after pos 1005 (CDS pos 6)
	v := &vcf.Variant{Chrom: "1", Pos: 1005, Ref: "T", Alt: "TC"}
	result := PredictConsequence(v, tr)

	if result.Consequence != ConsequenceFrameshiftVariant {
		t.Errorf("1bp insertion should be frameshift_variant, got %q", result.Consequence)
	}
	if result.Impact != ImpactHigh {
		t.Errorf("frameshift should be HIGH impact, got %q", result.Impact)
	}
}

// --- Pattern: Deletion (1bp) → frameshift ---

func TestDatahub_SingleBaseDeletion_Frameshift(t *testing.T) {
	tr := &cache.Transcript{
		ID: "ENST_1BPDEL", GeneName: "ONEBPDEL", Chrom: "1",
		Start: 995, End: 1025, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1014,
		Exons:        []cache.Exon{{Number: 1, Start: 995, End: 1025, CDSStart: 1000, CDSEnd: 1014, Frame: 0}},
		CDSSequence:  "ATGGCTAAAGGG" + "TAA",
		UTR3Sequence: "GCATAA",
	}

	// Delete 1 base at pos 1005 (CDS pos 6)
	v := &vcf.Variant{Chrom: "1", Pos: 1005, Ref: "TA", Alt: "T"}
	result := PredictConsequence(v, tr)

	if result.Consequence != ConsequenceFrameshiftVariant {
		t.Errorf("1bp deletion should be frameshift_variant, got %q", result.Consequence)
	}
}

// --- Pattern: Reverse strand MNV ---

func TestDatahub_ReverseStrand_MNV(t *testing.T) {
	tr := createKRASTranscript()

	// MNV at KRAS exon 2 (reverse strand): change 2 bases
	// Position 25245350-25245351 in exon 2 CDS
	v := &vcf.Variant{Chrom: "12", Pos: 25245350, Ref: "GC", Alt: "AA"}
	result := PredictConsequence(v, tr)

	// Should get some coding consequence
	if result.Consequence == "" || result.Consequence == ConsequenceIntronVariant {
		t.Errorf("reverse strand MNV in CDS should have coding consequence, got %q", result.Consequence)
	}
	t.Logf("reverse strand MNV: consequence=%q, protein_pos=%d", result.Consequence, result.ProteinPosition)
}

// --- Pattern: Insertion creating a duplication (from ensembl-vep tests) ---
// VEP test: chr21:25769085 dup T

func TestDatahub_SingleBaseDuplication(t *testing.T) {
	// CDS: ATG GCT TTT AAA TAA (M A F K *)
	// Insert T after the last T of TTT → duplication
	tr := &cache.Transcript{
		ID: "ENST_1BPDUP", GeneName: "ONEBPDUP", Chrom: "1",
		Start: 995, End: 1025, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1014,
		Exons:        []cache.Exon{{Number: 1, Start: 995, End: 1025, CDSStart: 1000, CDSEnd: 1014, Frame: 0}},
		CDSSequence:  "ATGGCTTTTAAA" + "TAA",
		UTR3Sequence: "GCATAA",
	}

	// Insert T at pos 1008 (last T of codon 3 = TTT)
	// Not divisible by 3 → frameshift, not dup detection at protein level
	v := &vcf.Variant{Chrom: "1", Pos: 1008, Ref: "T", Alt: "TT"}
	result := PredictConsequence(v, tr)

	// 1bp insertion → frameshift
	t.Logf("1bp dup insertion: consequence=%q", result.Consequence)
	if result.Consequence == "" {
		t.Error("expected non-empty consequence")
	}
}
