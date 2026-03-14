package annotate

import (
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Edge case transcripts used across tests in this file.

// createMultiExonForwardTranscript creates a forward-strand transcript with
// 3 coding exons plus UTRs, suitable for testing exon/intron boundaries,
// splice sites, and split codons.
//
// Layout (genomic):
//
//	5'UTR     Exon1(CDS)      Intron1     Exon2(CDS)    Intron2     Exon3(CDS)     3'UTR
//	[990-999] [1000-1010]  [1011-1999]  [2000-2008]  [2009-2999]  [3000-3011]  [3012-3020]
//
// CDS = 11 + 9 + 12 = 32 bases (not divisible by 3 for incomplete terminal codon testing)
// CDS: ATG GCT AAA | GAA GGG CCC | GGG TAA CCC T (last codon incomplete)
// Protein: M A K E G P G * ...
func createMultiExonForwardTranscript() *cache.Transcript {
	return &cache.Transcript{
		ID:       "ENST_EDGE_FWD",
		GeneID:   "ENSG_EDGE_FWD",
		GeneName: "EDGEFWD",
		Chrom:    "1",
		Start:    990,
		End:      3020,
		Strand:   1,
		Biotype:  "protein_coding",
		CDSStart: 1000,
		CDSEnd:   3011,
		Exons: []cache.Exon{
			{Number: 1, Start: 990, End: 1010, CDSStart: 1000, CDSEnd: 1010, Frame: 0},
			{Number: 2, Start: 2000, End: 2008, CDSStart: 2000, CDSEnd: 2008, Frame: 1},
			{Number: 3, Start: 3000, End: 3020, CDSStart: 3000, CDSEnd: 3011, Frame: 1},
		},
		// CDS (11 + 9 + 12 = 32bp):
		// Exon1 CDS: ATGGCTAAAGA (11bp, codons: ATG GCT AAA GA|)
		// Exon2 CDS: AGGGCCCGG  (9bp, split codon: |A from exon1 + AG from exon2 = GAA, then GGC CCG G|)
		// Exon3 CDS: GGTAACCCTTTT (12bp, split codon: |G from exon2 + GG from exon3 = GGG, then TAA CCC TTT T)
		CDSSequence:  "ATGGCTAAAGAAGGGCCCGGGGTAACCCTTTT",
		UTR3Sequence: "GCATAA",
	}
}

// createNonCodingTranscript creates a lincRNA transcript for testing
// non_coding_transcript_exon_variant.
func createNonCodingTranscript() *cache.Transcript {
	return &cache.Transcript{
		ID:       "ENST_LINCRNA",
		GeneID:   "ENSG_LINCRNA",
		GeneName: "TESTLINC",
		Chrom:    "1",
		Start:    5000,
		End:      6000,
		Strand:   1,
		Biotype:  "lincRNA",
		CDSStart: 0, // No CDS
		CDSEnd:   0,
		Exons: []cache.Exon{
			{Number: 1, Start: 5000, End: 5200},
			{Number: 2, Start: 5500, End: 5700},
		},
	}
}

// createMiRNATranscript creates a miRNA transcript.
func createMiRNATranscript() *cache.Transcript {
	return &cache.Transcript{
		ID:       "ENST_MIRNA",
		GeneID:   "ENSG_MIRNA",
		GeneName: "MIR21",
		Chrom:    "17",
		Start:    59841266,
		End:      59841337,
		Strand:   1,
		Biotype:  "miRNA",
		CDSStart: 0,
		CDSEnd:   0,
		Exons: []cache.Exon{
			{Number: 1, Start: 59841266, End: 59841337},
		},
	}
}

// --- Test 1: Exon/intron boundary variants ---

func TestEdge_ExonBoundaryLastBase(t *testing.T) {
	tr := createMultiExonForwardTranscript()

	// Variant at the very last base of exon 1 (pos 1010)
	v := &vcf.Variant{Chrom: "1", Pos: 1010, Ref: "A", Alt: "G"}
	result := PredictConsequence(v, tr)

	// Should be a coding consequence (in exon) AND splice_region_variant
	// since it's within 3bp of the exon boundary
	assert.Contains(t, result.Consequence, "splice_region_variant",
		"last base of exon should be splice region")
	assert.NotEqual(t, ConsequenceIntronVariant, result.Consequence,
		"should not be intronic — it's in the exon")
}

func TestEdge_ExonBoundaryFirstBase(t *testing.T) {
	tr := createMultiExonForwardTranscript()

	// Variant at the first base of exon 2 (pos 2000)
	v := &vcf.Variant{Chrom: "1", Pos: 2000, Ref: "A", Alt: "T"}
	result := PredictConsequence(v, tr)

	assert.Contains(t, result.Consequence, "splice_region_variant",
		"first base of exon should be splice region")
	assert.NotEqual(t, ConsequenceIntronVariant, result.Consequence)
}

func TestEdge_IntronFirstBase(t *testing.T) {
	tr := createMultiExonForwardTranscript()

	// First base of intron 1 (exon1 End + 1 = 1011)
	v := &vcf.Variant{Chrom: "1", Pos: 1011, Ref: "A", Alt: "G"}
	result := PredictConsequence(v, tr)

	// Forward strand: after exon end = splice donor
	assert.Equal(t, ConsequenceSpliceDonor, result.Consequence)
	assert.Equal(t, ImpactHigh, result.Impact)
}

func TestEdge_IntronLastBase(t *testing.T) {
	tr := createMultiExonForwardTranscript()

	// Last base before exon 2 (exon2 Start - 1 = 1999)
	v := &vcf.Variant{Chrom: "1", Pos: 1999, Ref: "A", Alt: "G"}
	result := PredictConsequence(v, tr)

	// Forward strand: before exon start = splice acceptor
	assert.Equal(t, ConsequenceSpliceAcceptor, result.Consequence)
	assert.Equal(t, ImpactHigh, result.Impact)
}

// --- Test 2: Non-coding transcript exon variant ---

func TestEdge_NonCodingExonVariant(t *testing.T) {
	tr := createNonCodingTranscript()

	v := &vcf.Variant{Chrom: "1", Pos: 5100, Ref: "A", Alt: "G"}
	result := PredictConsequence(v, tr)

	assert.Equal(t, ConsequenceNonCodingExon, result.Consequence)
	assert.Equal(t, ImpactModifier, result.Impact)
}

func TestEdge_MiRNAVariant(t *testing.T) {
	tr := createMiRNATranscript()

	v := &vcf.Variant{Chrom: "17", Pos: 59841300, Ref: "C", Alt: "T"}
	result := PredictConsequence(v, tr)

	assert.Equal(t, ConsequenceMatureMiRNA, result.Consequence)
}

// --- Test 3: Compound consequences ---

func TestEdge_StopGainedPlusSpliceRegion(t *testing.T) {
	// Create a transcript where a stop_gained can occur at a splice region position.
	// Forward strand, single exon CDS for simplicity but with 2 exons to have boundary.
	// Exon 1: 1000-1010 (CDS), Exon 2: 2000-2010 (CDS)
	// Place a variant at the 3rd-to-last base of exon 1 (splice region = within 3bp of boundary)
	// that creates a stop codon.

	tr := &cache.Transcript{
		ID: "ENST_COMPOUND", GeneName: "COMPOUND", Chrom: "1",
		Start: 990, End: 2020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 2010,
		Exons: []cache.Exon{
			{Number: 1, Start: 990, End: 1010, CDSStart: 1000, CDSEnd: 1010, Frame: 0},
			{Number: 2, Start: 2000, End: 2020, CDSStart: 2000, CDSEnd: 2010, Frame: 2},
		},
		// CDS: 11bp (exon1) + 11bp (exon2) = 22bp
		// ATG GCT AAA | TAG GCT AAA CC
		// But we want stop at a splice region. Exon 1 ends at 1010.
		// Position 1008 (3rd from end of exon) is splice region.
		// CDS pos at 1008 = 1008-1000+1 = 9, codon 3, pos 2
		// If codon 3 = AAA and we change pos 9 (3rd base) to G → AAG (Lys, not stop)
		// We need TAG/TAA/TGA at a splice region position.
		// Let's adjust: put a TAA-able codon at the boundary.
		// Codon at CDS pos 7-9 (genomic 1006-1008): make it TAC (Tyr)
		// Mutate CDS pos 8 (genomic 1007, 2nd from boundary = splice region) A→A stays same...
		// Let's just make CDS pos 9 (genomic 1008) create a stop.
		// CDS: ATG GCT TAC → mutate C at pos 9 to G → TAG = stop
		CDSSequence: "ATGGCTTACGCTAGGCTAAACC",
	}

	// pos 1008 = CDS pos 9 = 3rd base of codon 3 (TAC)
	// Mutate C→G: TAC → TAG = stop
	v := &vcf.Variant{Chrom: "1", Pos: 1008, Ref: "C", Alt: "G"}
	result := PredictConsequence(v, tr)

	assert.Contains(t, result.Consequence, ConsequenceStopGained,
		"should detect stop_gained")
	assert.Contains(t, result.Consequence, ConsequenceSpliceRegion,
		"should also have splice_region_variant (within 3bp of exon boundary)")
}

// --- Test 4: Start codon edge cases ---

func TestEdge_DeletionFrom5UTRIntoStartCodon(t *testing.T) {
	// Forward strand. 5'UTR at 990-999, CDS starts at 1000 (ATG).
	tr := &cache.Transcript{
		ID: "ENST_START_DEL", GeneName: "STARTDEL", Chrom: "1",
		Start: 985, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:       []cache.Exon{{Number: 1, Start: 985, End: 1020, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence: "ATGGCTAAATAA",
	}

	// Delete from UTR pos 998 into start codon: ref=GGATG (5 bases), alt=G
	// Spans 998-1002, which includes start codon at 1000-1002
	v := &vcf.Variant{Chrom: "1", Pos: 998, Ref: "GGATG", Alt: "G"}
	result := PredictConsequence(v, tr)

	assert.Equal(t, ConsequenceStartLost, result.Consequence,
		"deletion spanning from 5'UTR into start codon should be start_lost")
}

func TestEdge_MNVSpanningStartCodon(t *testing.T) {
	// Forward strand, MNV that changes the start codon.
	tr := &cache.Transcript{
		ID: "ENST_START_MNV", GeneName: "STARTMNV", Chrom: "1",
		Start: 995, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:       []cache.Exon{{Number: 1, Start: 995, End: 1020, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence: "ATGGCTAAATAA",
	}

	// Change ATG → CTG at start codon (CDS pos 1-3)
	v := &vcf.Variant{Chrom: "1", Pos: 1000, Ref: "ATG", Alt: "CTG"}
	result := PredictConsequence(v, tr)

	assert.Equal(t, ConsequenceStartLost, result.Consequence,
		"MNV changing start codon should be start_lost")
}

// --- Test 5: Stop codon edge cases ---

func TestEdge_MNVSpanningStopCodon(t *testing.T) {
	// Forward strand: CDS = ATG GCT AAA TAA (12bp: M A K *)
	// MNV changing stop codon TAA → TCA
	tr := &cache.Transcript{
		ID: "ENST_STOP_MNV", GeneName: "STOPMNV", Chrom: "1",
		Start: 995, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:        []cache.Exon{{Number: 1, Start: 995, End: 1020, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence:  "ATGGCTAAATAA",
		UTR3Sequence: "GCATAA",
	}

	// Change TAA → TCA at stop codon (genomic 1009-1011)
	v := &vcf.Variant{Chrom: "1", Pos: 1009, Ref: "TAA", Alt: "TCA"}
	result := PredictConsequence(v, tr)

	assert.Equal(t, ConsequenceStopLost, result.Consequence,
		"MNV changing stop codon should be stop_lost")
}

// --- Test 6: Split codon at exon boundary ---

func TestEdge_SplitCodonSNV(t *testing.T) {
	// Multi-exon transcript where a codon spans an exon boundary.
	// The createMultiExonForwardTranscript has split codons at exon junctions.
	tr := createMultiExonForwardTranscript()

	// Exon 1 CDS: pos 1000-1010 (11bp: ATGGCTAAAGA)
	// Last 2 bases of exon 1: GA (CDS pos 10-11, part of split codon)
	// First base of exon 2: A (CDS pos 12, completes the split codon GAA = Glu)
	// So codon 4 spans the exon1/exon2 boundary: GA|A

	// SNV at the last base of exon 1 (pos 1010, CDS pos 11)
	// This is the 2nd base of codon 4 (pos-in-codon = 1)
	// Codon 4 = GAA (Glu). Changing 2nd base: GAA → GXA
	v := &vcf.Variant{Chrom: "1", Pos: 1010, Ref: "A", Alt: "C"}
	result := PredictConsequence(v, tr)

	// Should get a coding consequence (not intron), even though the codon
	// is split across exons
	assert.NotEqual(t, ConsequenceIntronVariant, result.Consequence)
	assert.Greater(t, result.ProteinPosition, int64(0),
		"should have a protein position for split codon variant")

	// The variant is also at a splice region (last base of exon)
	assert.Contains(t, result.Consequence, "splice_region_variant")
}

func TestEdge_SplitCodonFirstBaseNextExon(t *testing.T) {
	tr := createMultiExonForwardTranscript()

	// First base of exon 2 (pos 2000, CDS pos 12)
	// This completes codon 4 (GAA → Glu)
	// Changing it: GA_ where _ is the mutant
	v := &vcf.Variant{Chrom: "1", Pos: 2000, Ref: "A", Alt: "T"}
	result := PredictConsequence(v, tr)

	assert.NotEqual(t, ConsequenceIntronVariant, result.Consequence)
	assert.Greater(t, result.ProteinPosition, int64(0))
	assert.Contains(t, result.Consequence, "splice_region_variant")
}

// --- Test 7: Incomplete terminal codon ---

func TestEdge_IncompleteTerminalCodon(t *testing.T) {
	// CDS length 32bp is not divisible by 3 (32/3 = 10 remainder 2)
	// The last 2 bases form an incomplete codon.
	tr := createMultiExonForwardTranscript()

	// Variant at the very last base of CDS (pos 3011, CDS pos 32)
	// This is in the incomplete terminal codon (bases 31-32, no 33rd base)
	v := &vcf.Variant{Chrom: "1", Pos: 3011, Ref: "T", Alt: "A"}
	result := PredictConsequence(v, tr)

	// Should still get some consequence (coding_sequence_variant or similar),
	// not a crash or empty result
	assert.NotEmpty(t, result.Consequence,
		"variant in incomplete terminal codon should still get a consequence")
}

// --- Test 8: Reverse strand boundary tests ---

func TestEdge_ReverseStrandSpliceConsistency(t *testing.T) {
	// Verify that splice donor/acceptor assignment is correct for reverse strand.
	// For reverse strand:
	//   - Before exon.Start (lower genomic coord side) = DONOR
	//   - After exon.End (higher genomic coord side) = ACCEPTOR
	tr := createKRASTranscript()

	// Exon 2: Start=25245274, End=25245395

	// Donor: Start-1 = 25245273
	vDonor := &vcf.Variant{Chrom: "12", Pos: 25245273, Ref: "A", Alt: "G"}
	rDonor := PredictConsequence(vDonor, tr)
	assert.Equal(t, ConsequenceSpliceDonor, rDonor.Consequence)

	// Acceptor: End+1 = 25245396
	vAcceptor := &vcf.Variant{Chrom: "12", Pos: 25245396, Ref: "A", Alt: "G"}
	rAcceptor := PredictConsequence(vAcceptor, tr)
	assert.Equal(t, ConsequenceSpliceAcceptor, rAcceptor.Consequence)
}

// --- Test 9: Deletion spanning entire CDS ---

func TestEdge_DeletionSpanningEntireCDS(t *testing.T) {
	tr := &cache.Transcript{
		ID: "ENST_WHOLECDS", GeneName: "WHOLECDS", Chrom: "1",
		Start: 990, End: 1020, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1011,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1020, CDSStart: 1000, CDSEnd: 1011, Frame: 0}},
		CDSSequence: "ATGGCTAAATAA",
	}

	// Delete from UTR through entire CDS and beyond
	v := &vcf.Variant{Chrom: "1", Pos: 995, Ref: "AAAAAAAAAAAAAAAAAAAAA", Alt: "A"}
	result := PredictConsequence(v, tr)

	// Should detect start_lost (deletion spans start codon)
	assert.Equal(t, ConsequenceStartLost, result.Consequence)
}

// --- Test 10: Adjacent exon boundaries (very short intron) ---

func TestEdge_VeryShortIntron(t *testing.T) {
	// Intron of only 4bp between exons — splice sites overlap with splice region.
	tr := &cache.Transcript{
		ID: "ENST_SHORTINTRON", GeneName: "SHORTINT", Chrom: "1",
		Start: 990, End: 1030, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1020,
		Exons: []cache.Exon{
			{Number: 1, Start: 990, End: 1005, CDSStart: 1000, CDSEnd: 1005, Frame: 0},
			{Number: 2, Start: 1010, End: 1030, CDSStart: 1010, CDSEnd: 1020, Frame: 0},
		},
		CDSSequence: "ATGGCTAAAGAATAA",
	}

	// Intron is 1006-1009 (4bp). Donor at 1006,1007. Acceptor at 1008,1009.
	vDonor := &vcf.Variant{Chrom: "1", Pos: 1006, Ref: "A", Alt: "G"}
	rDonor := PredictConsequence(vDonor, tr)
	assert.Equal(t, ConsequenceSpliceDonor, rDonor.Consequence)

	vAcceptor := &vcf.Variant{Chrom: "1", Pos: 1009, Ref: "A", Alt: "G"}
	rAcceptor := PredictConsequence(vAcceptor, tr)
	assert.Equal(t, ConsequenceSpliceAcceptor, rAcceptor.Consequence)
}
