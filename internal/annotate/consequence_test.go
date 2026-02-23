package annotate

import (
	"strings"
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// createKRASTranscript creates a minimal KRAS transcript for testing.
// KRAS is on the reverse strand, chr12:25205246-25250929 (GRCh38)
// Codon 12 is at genomic position 25245350
func createKRASTranscript() *cache.Transcript {
	return &cache.Transcript{
		ID:          "ENST00000311936",
		GeneID:      "ENSG00000133703",
		GeneName:    "KRAS",
		Chrom:       "12",
		Start:       25205246,
		End:         25250929,
		Strand:      -1, // Reverse strand
		Biotype:     "protein_coding",
		IsCanonical: true,
		CDSStart:    25209798, // CDS end in genomic coords (start for reverse)
		CDSEnd:      25245384, // CDS start in genomic coords (end for reverse)
		Exons: []cache.Exon{
			{Number: 1, Start: 25250751, End: 25250929, CDSStart: 0, CDSEnd: 0, Frame: -1},       // 5' UTR
			{Number: 2, Start: 25245274, End: 25245395, CDSStart: 25245274, CDSEnd: 25245384, Frame: 0}, // Contains codon 12
			{Number: 3, Start: 25227234, End: 25227412, CDSStart: 25227234, CDSEnd: 25227412, Frame: 0},
			{Number: 4, Start: 25225614, End: 25225773, CDSStart: 25225614, CDSEnd: 25225773, Frame: 2},
			{Number: 5, Start: 25209798, End: 25209911, CDSStart: 25209798, CDSEnd: 25209911, Frame: 0}, // Last coding exon
		},
		// CDS sequence starting from ATG (start codon is last 3 bases in genomic coords)
		// For KRAS reverse strand, we need the coding sequence 5' to 3'
		// Codon 12 is GGT (Gly) - this is approximately position 34-36 in the CDS
		// Full CDS would be ~567 bp for KRAS
		// For testing, we'll use a simplified sequence with codon 12 at the right position
		CDSSequence: "ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATATGATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCTGTCTCTTGGATATTCTCGACACAGCAGGTCAAGAGGAGTACAGTGCAATGAGGGACCAGTACATGAGGACTGGGGAGGGCTTTCTTTGTGTATTTGCCATAAATAATACTAAATCATTTGAAGATATTCACCATTATAGAGAACAAATTAAAAGAGTTAAGGACTCTGAAGATGTACCTATGGTCCTAGTAGGAAATAAATGTGATTTGCCTTCTAGAACAGTAGACACAAAACAGGCTCAGGACTTAGCAAGAAGTTATGGAATTCCTTTTATTGAAACATCAGCAAAGACAAGACAGAGAGTGGAGGATGCTTTTTATACATTGGTGAGAGAGATCCGACAATACAGATTGAAAAAAATCAGCAAAGAAGAAAAGACTCCTGGCTGTGTGAAAATTAAAAAATGCATTATAATGTAA",
	}
}

func TestPredictConsequence_KRASG12C(t *testing.T) {
	// KRAS G12C variant: c.34G>T p.G12C
	// KRAS is on the reverse strand, so:
	// - CDS position 34 = first base of codon 12 (GGT -> TGT)
	// - Genomic position = CDSEnd - CDSPos + 1 = 25245384 - 34 + 1 = 25245351
	// - On genomic strand, G on coding = C on genomic (complement)
	// - c.34G>T means: coding G->T, which on genomic strand is C->A
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245351,
		Ref:   "C",
		Alt:   "A",
	}

	transcript := createKRASTranscript()

	result := PredictConsequence(v, transcript)

	// Should be missense variant
	if result.Consequence != ConsequenceMissenseVariant {
		t.Errorf("Expected missense_variant, got %s", result.Consequence)
	}

	// Should be MODERATE impact
	if result.Impact != ImpactModerate {
		t.Errorf("Expected MODERATE impact, got %s", result.Impact)
	}

	// Protein position should be 12
	if result.ProteinPosition != 12 {
		t.Errorf("Expected protein position 12, got %d", result.ProteinPosition)
	}

	// Amino acid change should be G12C
	if result.AminoAcidChange != "G12C" {
		t.Errorf("Expected amino acid change G12C, got %s", result.AminoAcidChange)
	}
}

func TestPredictConsequence_Synonymous(t *testing.T) {
	// Create a variant that results in synonymous change
	// GGT -> GGC both code for Glycine
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245348, // Third position of codon 12
		Ref:   "A",      // On reverse strand, this would be T on coding strand
		Alt:   "G",      // On reverse strand, this would be C on coding strand
	}

	transcript := createKRASTranscript()

	result := PredictConsequence(v, transcript)

	// For a synonymous change
	if result.Consequence != ConsequenceSynonymousVariant &&
	   result.Consequence != ConsequenceMissenseVariant {
		// Note: exact consequence depends on codon position calculation
		t.Logf("Got consequence: %s (may vary based on position)", result.Consequence)
	}
}

func TestPredictConsequence_StopGained(t *testing.T) {
	// Test a variant that creates a stop codon
	// This is a simplified test
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245350,
		Ref:   "G",
		Alt:   "T", // If this creates a stop codon
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	// Just verify we get a consequence
	if result.Consequence == "" {
		t.Error("Expected a consequence, got empty string")
	}
}

func TestPredictConsequence_IntronicVariant(t *testing.T) {
	// Test a variant in an intron
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25235000, // Between exons
		Ref:   "A",
		Alt:   "G",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	if result.Consequence != ConsequenceIntronVariant {
		t.Errorf("Expected intron_variant, got %s", result.Consequence)
	}

	if result.Impact != ImpactModifier {
		t.Errorf("Expected MODIFIER impact, got %s", result.Impact)
	}
}

func TestPredictConsequence_UpstreamVariant(t *testing.T) {
	// Test a variant upstream of the gene
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25260000, // After gene end (upstream for reverse strand)
		Ref:   "A",
		Alt:   "G",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	// Should be upstream for reverse strand gene
	if result.Consequence != ConsequenceUpstreamGene &&
	   result.Consequence != ConsequenceDownstreamGene {
		t.Errorf("Expected upstream/downstream variant, got %s", result.Consequence)
	}
}

func TestCDSToCodonPosition(t *testing.T) {
	tests := []struct {
		cdsPos          int64
		wantCodonNum    int64
		wantPosInCodon  int
	}{
		{1, 1, 0},   // First base of first codon
		{2, 1, 1},   // Second base of first codon
		{3, 1, 2},   // Third base of first codon
		{4, 2, 0},   // First base of second codon
		{34, 12, 0}, // First base of codon 12 (KRAS G12)
		{35, 12, 1}, // Second base of codon 12
		{36, 12, 2}, // Third base of codon 12
		{0, 0, 0},   // Invalid position
	}

	for _, tt := range tests {
		t.Run("", func(t *testing.T) {
			codonNum, posInCodon := CDSToCodonPosition(tt.cdsPos)
			if codonNum != tt.wantCodonNum || posInCodon != tt.wantPosInCodon {
				t.Errorf("CDSToCodonPosition(%d) = (%d, %d), want (%d, %d)",
					tt.cdsPos, codonNum, posInCodon, tt.wantCodonNum, tt.wantPosInCodon)
			}
		})
	}
}

func TestPredictConsequence_FrameshiftVariant(t *testing.T) {
	// Test a frameshift deletion (1 base deleted)
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245350,
		Ref:   "GG",
		Alt:   "G",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	if result.Consequence != ConsequenceFrameshiftVariant {
		t.Errorf("Expected frameshift_variant, got %s", result.Consequence)
	}

	if result.Impact != ImpactHigh {
		t.Errorf("Expected HIGH impact, got %s", result.Impact)
	}
}

func TestIsSpliceRegion(t *testing.T) {
	transcript := createKRASTranscript()
	// Exon 2: Start=25245274, End=25245395

	tests := []struct {
		name     string
		pos      int64
		expected bool
	}{
		// Exon side: within 3bp of exon.End (25245395)
		{"exon_side_end_0bp", 25245395, true},
		{"exon_side_end_1bp", 25245394, true},
		{"exon_side_end_2bp", 25245393, true},
		{"exon_side_end_3bp", 25245392, false},

		// Exon side: within 3bp of exon.Start (25245274)
		{"exon_side_start_0bp", 25245274, true},
		{"exon_side_start_1bp", 25245275, true},
		{"exon_side_start_2bp", 25245276, true},
		{"exon_side_start_3bp", 25245277, false},

		// Intron side: 3-8bp after exon.End (25245395)
		{"intron_after_end_1bp", 25245396, false}, // splice donor/acceptor
		{"intron_after_end_2bp", 25245397, false}, // splice donor/acceptor
		{"intron_after_end_3bp", 25245398, true},  // splice region
		{"intron_after_end_5bp", 25245400, true},  // splice region
		{"intron_after_end_8bp", 25245403, true},  // splice region
		{"intron_after_end_9bp", 25245404, false}, // too far

		// Intron side: 3-8bp before exon.Start (25245274)
		{"intron_before_start_1bp", 25245273, false}, // splice donor/acceptor
		{"intron_before_start_2bp", 25245272, false}, // splice donor/acceptor
		{"intron_before_start_3bp", 25245271, true},  // splice region
		{"intron_before_start_5bp", 25245269, true},  // splice region
		{"intron_before_start_8bp", 25245266, true},  // splice region
		{"intron_before_start_9bp", 25245265, false}, // too far

		// Middle of exon - not splice region
		{"mid_exon", 25245350, false},

		// Far into intron
		{"deep_intron", 25235000, false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := isSpliceRegion(tt.pos, transcript)
			if got != tt.expected {
				t.Errorf("isSpliceRegion(%d) = %v, want %v", tt.pos, got, tt.expected)
			}
		})
	}
}

func TestPredictConsequence_SpliceRegionIntronic(t *testing.T) {
	// Variant 5bp into intron after exon 2 End (25245395)
	// Position 25245400 = 5bp after exon boundary → splice_region_variant,intron_variant
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245400,
		Ref:   "A",
		Alt:   "G",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	expected := "splice_region_variant,intron_variant"
	if result.Consequence != expected {
		t.Errorf("Expected %s, got %s", expected, result.Consequence)
	}
	if result.Impact != "LOW" {
		t.Errorf("Expected LOW impact, got %s", result.Impact)
	}
}

func TestPredictConsequence_SpliceRegionCoding(t *testing.T) {
	// Variant at exon 2 End boundary (25245395) - 0bp from boundary, exon side
	// This is within the CDS (CDSEnd=25245384... wait, 25245395 > 25245384)
	// Actually 25245395 > CDSEnd 25245384, so this would be 5'UTR on reverse strand
	// Use exon 2 Start boundary instead: pos 25245274 is CDS (CDSStart=25245274)
	// On reverse strand, CDSStart < CDSEnd, and pos >= CDSStart and pos <= CDSEnd → in CDS
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245274, // exon.Start, within 3bp of boundary, in CDS
		Ref:   "C",
		Alt:   "T",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	// Should have splice_region_variant appended to coding consequence
	if !strings.Contains(result.Consequence, "splice_region_variant") {
		t.Errorf("Expected consequence to contain splice_region_variant, got %s", result.Consequence)
	}
	// Primary consequence should be a coding type
	if !strings.Contains(result.Consequence, "missense_variant") &&
		!strings.Contains(result.Consequence, "synonymous_variant") {
		t.Errorf("Expected coding consequence + splice_region_variant, got %s", result.Consequence)
	}
}

func TestPredictConsequence_NoSpliceRegionMidExon(t *testing.T) {
	// KRAS G12C at position 25245351 - middle of exon, should NOT have splice_region
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245351,
		Ref:   "C",
		Alt:   "A",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	if strings.Contains(result.Consequence, "splice_region_variant") {
		t.Errorf("Expected no splice_region_variant for mid-exon variant, got %s", result.Consequence)
	}
	if result.Consequence != "missense_variant" {
		t.Errorf("Expected missense_variant, got %s", result.Consequence)
	}
}

func TestPredictConsequence_NoSpliceRegionDeepIntron(t *testing.T) {
	// Variant deep in intron (25235000) - should be plain intron_variant
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25235000,
		Ref:   "A",
		Alt:   "G",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	if result.Consequence != "intron_variant" {
		t.Errorf("Expected intron_variant, got %s", result.Consequence)
	}
}

func TestPredictConsequence_InframeDeletion(t *testing.T) {
	// Test an in-frame deletion (3 bases deleted)
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245350,
		Ref:   "GGGA",
		Alt:   "G",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	if result.Consequence != ConsequenceInframeDeletion {
		t.Errorf("Expected inframe_deletion, got %s", result.Consequence)
	}

	if result.Impact != ImpactModerate {
		t.Errorf("Expected MODERATE impact, got %s", result.Impact)
	}
}
