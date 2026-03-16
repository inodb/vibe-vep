package annotate

// Tests derived from the official HGVS nomenclature specification
// (https://github.com/HGVSNomenclature/hgvs-nomenclature)
// covering HGVSc and HGVSp corner cases.

import (
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// === MNV / Delins tests ===
// HGVS rule: substitution affects exactly ONE nucleotide.
// Two or more consecutive nucleotide changes must use delins.

func TestHGVS_MNV_ForwardStrand_CDS(t *testing.T) {
	// Two-base MNV in CDS on forward strand: CC>TT at CDS positions 4-5
	// Expected: c.4_5delinsTT (NOT c.4CC>TT)
	cds := "ATGCCGAAAGGGGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
	transcript := &cache.Transcript{
		ID: "ENST_MNV_FWD", GeneName: "MNVFWD", Chrom: "1",
		Start: 990, End: 1210, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	v := &vcf.Variant{Chrom: "1", Pos: 1003, Ref: "CC", Alt: "TT"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.4_5delinsTT", hgvsc)
}

func TestHGVS_MNV_ForwardStrand_ThreeBase(t *testing.T) {
	// Three-base MNV spanning one codon: AGG>TCC at CDS positions 4-6
	// Expected: c.4_6delinsTCC
	cds := "ATGAGGAAAGGGGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
	transcript := &cache.Transcript{
		ID: "ENST_MNV3", GeneName: "MNV3", Chrom: "1",
		Start: 990, End: 1210, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	v := &vcf.Variant{Chrom: "1", Pos: 1003, Ref: "AGG", Alt: "TCC"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.4_6delinsTCC", hgvsc)
}

func TestHGVS_MNV_ReverseStrand_CDS(t *testing.T) {
	// Two-base MNV on KRAS (reverse strand)
	// Genomic: CC>TT at pos 25245350-25245351
	// Coding strand: RC(CC)=GG, RC(TT)=AA
	// CDS pos 34 maps to genomic 25245351, pos 35 maps to 25245350
	// Expected: c.34_35delinsAA
	transcript := createKRASTranscript()

	v := &vcf.Variant{Chrom: "12", Pos: 25245350, Ref: "CC", Alt: "TT"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.34_35delinsAA", hgvsc)
}

func TestHGVS_MNV_ReverseStrand_UTR(t *testing.T) {
	// TERT-like promoter variant: two-base MNV in 5'UTR on reverse strand
	// This is the reported TERT c.-139_-138delinsTT case
	//
	// Build a reverse strand transcript where:
	// - CDS starts at genomic position 1100 (CDSEnd for reverse)
	// - 5'UTR extends from 1101 upward
	// - Variant at genomic 1240-1241 is in 5'UTR
	cds := "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
	transcript := &cache.Transcript{
		ID: "ENST_TERT", GeneName: "TERT", Chrom: "5",
		Start: 990, End: 1300, Strand: -1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1300, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	// Two positions in 5'UTR: genomic 1240 and 1241
	// On reverse strand, these are upstream of CDSEnd (1101), so 5'UTR
	// Distance from CDSEnd: 1240-1101 = 139 exonic bases → c.-139
	//                       1241-1101 = 140 exonic bases → c.-140
	// Reverse strand: genomic 1241 (higher) → more 5' → c.-140
	//                 genomic 1240 (lower)  → less 5' → c.-139
	// After swap: c.-140_-139
	v := &vcf.Variant{Chrom: "5", Pos: 1240, Ref: "CC", Alt: "TT"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	// RC(CC)=GG, RC(TT)=AA
	// Reverse strand swap: positions are -140 and -139
	assert.Equal(t, "c.-140_-139delinsAA", hgvsc)
}

func TestHGVS_MNV_UTR_ForwardStrand(t *testing.T) {
	// MNV in 5'UTR on forward strand
	// Positions 998-999 are 5'UTR (CDSStart=1000)
	// c.-2 and c.-1
	// Expected: c.-2_-1delinsTT
	cds := "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
	transcript := &cache.Transcript{
		ID: "ENST_UTR_FWD", GeneName: "UTRFWD", Chrom: "1",
		Start: 990, End: 1210, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	v := &vcf.Variant{Chrom: "1", Pos: 998, Ref: "CC", Alt: "TT"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.-2_-1delinsTT", hgvsc)
}

// === SNV tests ===
// Single nucleotide substitution: c.123A>G

func TestHGVS_SNV_ForwardStrand(t *testing.T) {
	transcript := createDupTestTranscript()
	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "G", Alt: "C"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.3G>C", hgvsc)
}

func TestHGVS_SNV_ReverseStrand(t *testing.T) {
	// KRAS G12C: genomic C>A → coding G>T
	transcript := createKRASTranscript()
	v := &vcf.Variant{Chrom: "12", Pos: 25245351, Ref: "C", Alt: "A"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.34G>T", hgvsc)
}

func TestHGVS_SNV_FivePrimeUTR(t *testing.T) {
	// SNV in 5'UTR: c.-2A>G
	transcript := createForwardTranscript()
	v := &vcf.Variant{Chrom: "1", Pos: 998, Ref: "A", Alt: "G"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.-2A>G", hgvsc)
}

func TestHGVS_SNV_ThreePrimeUTR(t *testing.T) {
	// SNV in 3'UTR on forward strand
	// CDSEnd=1180, position 1182 → c.*2
	transcript := createForwardTranscript()
	v := &vcf.Variant{Chrom: "1", Pos: 1182, Ref: "A", Alt: "G"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Contains(t, hgvsc, "c.*")
	assert.Contains(t, hgvsc, ">")
}

// === Deletion tests ===
// HGVS: deletions always described at most 3' position

func TestHGVS_Deletion_SingleBase(t *testing.T) {
	// Single base deletion with 3' shift
	// CDS: "ATGAAAACCC..." delete first A → shift to last A in run
	cds := "ATGAAAACCCATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGA"
	transcript := &cache.Transcript{
		ID: "ENST_DEL1", GeneName: "DEL1", Chrom: "1",
		Start: 990, End: 1210, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	// Delete A at CDS pos 4 → should shift to CDS pos 7 (last A in AAAA run)
	v := &vcf.Variant{Chrom: "1", Pos: 1003, Ref: "AA", Alt: "A"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.7del", hgvsc)
}

func TestHGVS_Deletion_MultiBase_WithShift(t *testing.T) {
	// Multi-base deletion with 3' shift in repeat
	cds := "ATGATGATGCCCATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGA"
	transcript := &cache.Transcript{
		ID: "ENST_DEL_MULTI", GeneName: "DELMULTI", Chrom: "1",
		Start: 990, End: 1210, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	// Delete ATG at CDS 4-6 → shifts to 7-9 (last ATG in repeat)
	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "GATG", Alt: "G"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.7_9del", hgvsc)
}

func TestHGVS_Deletion_NoShift(t *testing.T) {
	// Deletion where next base differs → no shift
	cds := "ATGCCCAAAGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
	transcript := &cache.Transcript{
		ID: "ENST_DEL_NS", GeneName: "DELNS", Chrom: "1",
		Start: 990, End: 1210, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	// Delete CCC at CDS 4-6, next base A≠C → no shift
	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "GCCC", Alt: "G"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.4_6del", hgvsc)
}

// === Insertion tests ===
// HGVS: insertions use two flanking positions (c.123_124insXXX)

func TestHGVS_Insertion_Plain(t *testing.T) {
	// Non-dup insertion: inserted bases don't match adjacent
	transcript := createDupTestTranscript()
	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "G", Alt: "GCC"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Contains(t, hgvsc, "ins")
	assert.NotContains(t, hgvsc, "dup")
}

// === Duplication tests ===
// HGVS: tandem duplications must use dup, not ins

func TestHGVS_Dup_SingleBase(t *testing.T) {
	// Insert G that matches preceding base → dup
	transcript := createDupTestTranscript()
	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "G", Alt: "GG"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.3dup", hgvsc)
}

func TestHGVS_Dup_MultiBase(t *testing.T) {
	// Insert ATG that matches preceding bases → multi-base dup
	transcript := createDupTestTranscript()
	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "G", Alt: "GATG"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.3_5dup", hgvsc)
}

func TestHGVS_Dup_WithThreePrimeShift(t *testing.T) {
	// Insertion in poly-A: dup shifts to 3' end
	cds := "ATGAAACCCATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGA"
	transcript := &cache.Transcript{
		ID: "ENST_DUP_SHIFT", GeneName: "DUPSHIFT", Chrom: "1",
		Start: 990, End: 1210, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	// Insert A at pos 1003 (CDS 4='A'), shifts to last A → c.6dup
	v := &vcf.Variant{Chrom: "1", Pos: 1003, Ref: "A", Alt: "AA"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.6dup", hgvsc)
}

func TestHGVS_Dup_ReverseStrand(t *testing.T) {
	// KRAS dup on reverse strand
	transcript := createKRASTranscript()
	v := &vcf.Variant{Chrom: "12", Pos: 25245351, Ref: "C", Alt: "CC"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.35dup", hgvsc)
}

// === Delins tests (indel type) ===
// HGVS: deletion+insertion where ref and alt differ in length

func TestHGVS_Delins_Forward(t *testing.T) {
	// Delete 3 bases, insert 1: c.5_7delinsG
	cds := "ATGCCCAAAGGGGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
	transcript := &cache.Transcript{
		ID: "ENST_DI", GeneName: "DI", Chrom: "1",
		Start: 990, End: 1210, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	v := &vcf.Variant{Chrom: "1", Pos: 1003, Ref: "CCCA", Alt: "CG"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.5_7delinsG", hgvsc)
}

func TestHGVS_Delins_Reverse(t *testing.T) {
	// Delins on KRAS reverse strand
	transcript := createKRASTranscript()
	v := &vcf.Variant{Chrom: "12", Pos: 25245348, Ref: "GCCC", Alt: "GA"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.34_36delinsT", hgvsc)
}

func TestHGVS_Delins_SuffixClipping(t *testing.T) {
	// VCF representation with shared suffix that should be clipped
	// REF=ATGCA, ALT=ACA → after prefix/suffix clip: delete TG → c.2_3del
	cds := "ATGATGCAAGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGA"
	transcript := &cache.Transcript{
		ID: "ENST_SC", GeneName: "SC", Chrom: "1",
		Start: 990, End: 1210, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	v := &vcf.Variant{Chrom: "1", Pos: 1000, Ref: "ATGCA", Alt: "ACA"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.2_3del", hgvsc)
}

// === Intronic notation tests ===
// HGVS: c.87+1 (5' splice), c.88-1 (3' splice)

func TestHGVS_Intronic_FivePrimeSplice(t *testing.T) {
	// Position in intron near 5' splice site (after exon 1)
	// Exon 1 ends at 1020, intron starts at 1021
	transcript := createForwardTranscript()
	v := &vcf.Variant{Chrom: "1", Pos: 1021, Ref: "A", Alt: "G"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	// Should be c.21+1A>G (last CDS position of exon 1 is 21, then +1)
	assert.Contains(t, hgvsc, "+")
}

func TestHGVS_Intronic_ThreePrimeSplice(t *testing.T) {
	// Position in intron near 3' splice site (before exon 2)
	// Exon 2 starts at 1100, intron position 1099
	transcript := createForwardTranscript()
	v := &vcf.Variant{Chrom: "1", Pos: 1099, Ref: "A", Alt: "G"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	// Should be c.22-1A>G (first CDS position of exon 2 is 22, then -1)
	assert.Contains(t, hgvsc, "-")
}

// === Splice junction duplication tests ===

func TestHGVS_SpliceJunction_Dup_Forward(t *testing.T) {
	// Insertion at intronic position duplicating first exonic CDS base
	cds := "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
	transcript := &cache.Transcript{
		ID: "ENST_SJ", GeneName: "SJ", Chrom: "1",
		Start: 990, End: 1200, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1180,
		Exons: []cache.Exon{
			{Number: 1, Start: 990, End: 1020, CDSStart: 1000, CDSEnd: 1020, Frame: 0},
			{Number: 2, Start: 1100, End: 1200, CDSStart: 1100, CDSEnd: 1180, Frame: 0},
		},
		CDSSequence: cds,
	}

	// Insert A at intronic pos 1099 → matches CDS[21]='A' → c.22dup
	v := &vcf.Variant{Chrom: "1", Pos: 1099, Ref: "T", Alt: "TA"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.22dup", hgvsc)
}

func TestHGVS_SpliceJunction_Insert_NonMatching(t *testing.T) {
	// Insertion at splice junction that doesn't match → plain ins, not dup
	cds := "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
	transcript := &cache.Transcript{
		ID: "ENST_SJ2", GeneName: "SJ2", Chrom: "1",
		Start: 990, End: 1200, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1180,
		Exons: []cache.Exon{
			{Number: 1, Start: 990, End: 1020, CDSStart: 1000, CDSEnd: 1020, Frame: 0},
			{Number: 2, Start: 1100, End: 1200, CDSStart: 1100, CDSEnd: 1180, Frame: 0},
		},
		CDSSequence: cds,
	}

	// Insert C at intronic pos 1099, CDS[21]='A' → no match → ins
	v := &vcf.Variant{Chrom: "1", Pos: 1099, Ref: "T", Alt: "TC"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.NotContains(t, hgvsc, "dup")
	assert.Contains(t, hgvsc, "ins")
}

// === HGVSp tests derived from HGVS spec ===

func TestHGVSp_Missense(t *testing.T) {
	// p.Trp24Cys
	result := &ConsequenceResult{
		Consequence:     ConsequenceMissenseVariant,
		ProteinPosition: 24,
		RefAA:           'W',
		AltAA:           'C',
	}
	assert.Equal(t, "p.Trp24Cys", FormatHGVSp(result))
}

func TestHGVSp_Nonsense(t *testing.T) {
	// p.Trp24Ter (stop gained)
	result := &ConsequenceResult{
		Consequence:     ConsequenceStopGained,
		ProteinPosition: 24,
		RefAA:           'W',
		AltAA:           '*',
	}
	assert.Equal(t, "p.Trp24Ter", FormatHGVSp(result))
}

func TestHGVSp_Silent(t *testing.T) {
	// p.Cys188= (synonymous)
	result := &ConsequenceResult{
		Consequence:     ConsequenceSynonymousVariant,
		ProteinPosition: 188,
		RefAA:           'C',
		AltAA:           'C',
	}
	assert.Equal(t, "p.Cys188=", FormatHGVSp(result))
}

func TestHGVS_StartLost(t *testing.T) {
	// p.Met1? (start lost — always this format)
	result := &ConsequenceResult{
		Consequence: ConsequenceStartLost,
	}
	assert.Equal(t, "p.Met1?", FormatHGVSp(result))
}

func TestHGVSp_StopLost_WithExtension(t *testing.T) {
	// p.Ter110Glnext*17 (stop codon to Gln, new stop at +17)
	result := &ConsequenceResult{
		Consequence:     ConsequenceStopLost,
		ProteinPosition: 110,
		RefAA:           '*',
		AltAA:           'Q',
		StopLostExtDist: 17,
	}
	assert.Equal(t, "p.Ter110Glnext*17", FormatHGVSp(result))
}

func TestHGVSp_StopLost_UnknownExtension(t *testing.T) {
	// p.Ter327Argext*? (no new stop found)
	result := &ConsequenceResult{
		Consequence:     ConsequenceStopLost,
		ProteinPosition: 327,
		RefAA:           '*',
		AltAA:           'R',
	}
	assert.Equal(t, "p.Ter327Argext*?", FormatHGVSp(result))
}

func TestHGVSp_StopRetained(t *testing.T) {
	// p.Ter130= (synonymous stop)
	result := &ConsequenceResult{
		Consequence:     ConsequenceStopRetained,
		ProteinPosition: 130,
		RefAA:           '*',
		AltAA:           '*',
	}
	assert.Equal(t, "p.Ter130=", FormatHGVSp(result))
}

func TestHGVSp_Frameshift_WithStop(t *testing.T) {
	// p.Arg97ProfsTer23
	result := &ConsequenceResult{
		Consequence:        ConsequenceFrameshiftVariant,
		ProteinPosition:    97,
		RefAA:              'R',
		AltAA:              'P',
		FrameshiftStopDist: 23,
	}
	assert.Equal(t, "p.Arg97ProfsTer23", FormatHGVSp(result))
}

func TestHGVSp_Frameshift_ShortestPossible(t *testing.T) {
	// p.Glu5ValfsTer5 (new stop relatively close)
	result := &ConsequenceResult{
		Consequence:        ConsequenceFrameshiftVariant,
		ProteinPosition:    5,
		RefAA:              'E',
		AltAA:              'V',
		FrameshiftStopDist: 5,
	}
	assert.Equal(t, "p.Glu5ValfsTer5", FormatHGVSp(result))
}

func TestHGVSp_InframeDeletion_Single(t *testing.T) {
	// p.Val7del
	result := &ConsequenceResult{
		Consequence:     ConsequenceInframeDeletion,
		ProteinPosition: 7,
		RefAA:           'V',
	}
	assert.Equal(t, "p.Val7del", FormatHGVSp(result))
}

func TestHGVSp_InframeDeletion_Range(t *testing.T) {
	// p.Lys23_Val25del
	result := &ConsequenceResult{
		Consequence:        ConsequenceInframeDeletion,
		ProteinPosition:    23,
		ProteinEndPosition: 25,
		RefAA:              'K',
		EndAA:              'V',
	}
	assert.Equal(t, "p.Lys23_Val25del", FormatHGVSp(result))
}

func TestHGVSp_InframeInsertion_Single(t *testing.T) {
	// p.His4_Gln5insAla
	result := &ConsequenceResult{
		Consequence:     ConsequenceInframeInsertion,
		ProteinPosition: 4,
		RefAA:           'H',
		EndAA:           'Q',
		InsertedAAs:     "A",
	}
	assert.Equal(t, "p.His4_Gln5insAla", FormatHGVSp(result))
}

func TestHGVSp_InframeInsertion_Multiple(t *testing.T) {
	// p.Lys2_Gly3insGlnSerLys
	result := &ConsequenceResult{
		Consequence:     ConsequenceInframeInsertion,
		ProteinPosition: 2,
		RefAA:           'K',
		EndAA:           'G',
		InsertedAAs:     "QSK",
	}
	assert.Equal(t, "p.Lys2_Gly3insGlnSerLys", FormatHGVSp(result))
}

func TestHGVSp_Duplication_Single(t *testing.T) {
	// p.Val7dup
	result := &ConsequenceResult{
		Consequence:     ConsequenceInframeInsertion,
		ProteinPosition: 7,
		RefAA:           'V',
		InsertedAAs:     "V",
		IsDup:           true,
	}
	assert.Equal(t, "p.Val7dup", FormatHGVSp(result))
}

func TestHGVSp_Duplication_Range(t *testing.T) {
	// p.Lys23_Val25dup
	result := &ConsequenceResult{
		Consequence:        ConsequenceInframeInsertion,
		ProteinPosition:    23,
		ProteinEndPosition: 25,
		RefAA:              'K',
		EndAA:              'V',
		InsertedAAs:        "KAV",
		IsDup:              true,
	}
	assert.Equal(t, "p.Lys23_Val25dup", FormatHGVSp(result))
}

func TestHGVSp_Delins_SingleAA(t *testing.T) {
	// p.Cys28delinsTrpVal
	result := &ConsequenceResult{
		Consequence:     ConsequenceInframeDeletion,
		ProteinPosition: 28,
		RefAA:           'C',
		InsertedAAs:     "WV",
	}
	assert.Equal(t, "p.Cys28delinsTrpVal", FormatHGVSp(result))
}

func TestHGVSp_Delins_Range(t *testing.T) {
	// p.Cys28_Lys29delinsTrp
	result := &ConsequenceResult{
		Consequence:        ConsequenceInframeDeletion,
		ProteinPosition:    28,
		ProteinEndPosition: 29,
		RefAA:              'C',
		EndAA:              'K',
		InsertedAAs:        "W",
	}
	assert.Equal(t, "p.Cys28_Lys29delinsTrp", FormatHGVSp(result))
}

func TestHGVSp_MissenseDelins_MultiCodon(t *testing.T) {
	// Multi-codon MNV: p.Arg49Trp (from c.145_147delinsTGG)
	// When MNV spans multiple codons with multiple AA changes → delins
	result := &ConsequenceResult{
		Consequence:        ConsequenceMissenseVariant,
		ProteinPosition:    49,
		ProteinEndPosition: 50,
		RefAA:              'R',
		EndAA:              'G',
		InsertedAAs:        "WE",
		IsDelIns:           true,
	}
	assert.Equal(t, "p.Arg49_Gly50delinsTrpGlu", FormatHGVSp(result))
}

// === Upstream/Downstream should return empty HGVSc ===

func TestHGVS_Upstream_Empty(t *testing.T) {
	transcript := createForwardTranscript()
	v := &vcf.Variant{Chrom: "1", Pos: 500, Ref: "A", Alt: "G"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Empty(t, hgvsc)
}

func TestHGVS_Downstream_Empty(t *testing.T) {
	transcript := createForwardTranscript()
	v := &vcf.Variant{Chrom: "1", Pos: 5000, Ref: "A", Alt: "G"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Empty(t, hgvsc)
}

// === Nil/edge case tests ===

func TestHGVS_NilInputs(t *testing.T) {
	assert.Empty(t, FormatHGVSc(nil, nil, nil))
	assert.Empty(t, FormatHGVSc(&vcf.Variant{}, nil, nil))
	assert.Empty(t, FormatHGVSc(nil, &cache.Transcript{}, nil))
}

func TestHGVSp_NonCoding_Empty(t *testing.T) {
	result := &ConsequenceResult{
		Consequence: ConsequenceIntronVariant,
	}
	assert.Empty(t, FormatHGVSp(result))
}

// === No-change (ref==alt) tests ===
// HGVS rule: nucleotides tested and found not changed are described as c.123=

func TestHGVS_NoChange_SNV(t *testing.T) {
	transcript := createDupTestTranscript()
	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "G", Alt: "G"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.3=", hgvsc)
}

func TestHGVS_NoChange_UTR(t *testing.T) {
	transcript := createForwardTranscript()
	v := &vcf.Variant{Chrom: "1", Pos: 998, Ref: "A", Alt: "A"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.-2=", hgvsc)
}

// === Insertion at end of CDS ===
// HGVS: c.{lastCDS}_*1insX when insertion is between CDS and 3'UTR

func TestHGVS_Insertion_EndOfCDS(t *testing.T) {
	// Single exon, CDS = 9 bases. Insert after last CDS base.
	cds := "ATGATGATG"
	transcript := &cache.Transcript{
		ID: "ENST_END", GeneName: "END", Chrom: "1",
		Start: 990, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1008,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1020, CDSStart: 1000, CDSEnd: 1008, Frame: 0}},
		CDSSequence: cds,
	}
	transcript.BuildCDSIndex()

	// Insert C after last CDS base (CDS pos 9 = G at genomic 1008)
	v := &vcf.Variant{Chrom: "1", Pos: 1008, Ref: "G", Alt: "GC"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.9_*1insC", hgvsc)
}

func TestHGVS_Insertion_EndOfCDS_WithShift(t *testing.T) {
	// CDS ends with GG. Insert G at second-to-last position → shifts to end.
	cds := "ATGATGATGG"
	transcript := &cache.Transcript{
		ID: "ENST_END2", GeneName: "END2", Chrom: "1",
		Start: 990, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1009,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1020, CDSStart: 1000, CDSEnd: 1009, Frame: 0}},
		CDSSequence: cds,
	}
	transcript.BuildCDSIndex()

	// Insert G at CDS pos 9 (G), should shift to pos 10 (also G) → dup
	v := &vcf.Variant{Chrom: "1", Pos: 1008, Ref: "G", Alt: "GG"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	// CDS[8]='G', CDS[9]='G'. Insert G at anchor 8, shifts to 9.
	// Dup check: CDS[9]='G' == inserted 'G' → c.10dup
	assert.Equal(t, "c.10dup", hgvsc)
}
