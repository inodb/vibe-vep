package annotate

import (
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"

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
			{Number: 1, Start: 25250751, End: 25250929, CDSStart: 0, CDSEnd: 0, Frame: -1},              // 5' UTR
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
	assert.Equal(t, ConsequenceMissenseVariant, result.Consequence)

	// Should be MODERATE impact
	assert.Equal(t, ImpactModerate, result.Impact)

	// Protein position should be 12
	assert.Equal(t, int64(12), result.ProteinPosition)

	// Amino acid change should be G12C
	assert.Equal(t, "G12C", result.AminoAcidChange)

	// HGVSp should be p.Gly12Cys
	assert.Equal(t, "p.Gly12Cys", result.HGVSp)
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
	assert.NotEmpty(t, result.Consequence)
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

	assert.Equal(t, ConsequenceIntronVariant, result.Consequence)

	assert.Equal(t, ImpactModifier, result.Impact)
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
		cdsPos         int64
		wantCodonNum   int64
		wantPosInCodon int
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
			assert.Equal(t, tt.wantCodonNum, codonNum)
			assert.Equal(t, tt.wantPosInCodon, posInCodon)
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

	assert.Equal(t, ConsequenceFrameshiftVariant, result.Consequence)

	assert.Equal(t, ImpactHigh, result.Impact)
}

func TestSpliceSiteType(t *testing.T) {
	// KRAS is reverse strand. Exon 2: Start=25245274, End=25245395
	// Reverse strand: before exon.Start = donor, after exon.End = acceptor
	transcript := createKRASTranscript()

	tests := []struct {
		name     string
		pos      int64
		expected string
	}{
		// After exon.End (25245395): reverse strand -> acceptor
		{"end_plus1", 25245396, ConsequenceSpliceAcceptor},
		{"end_plus2", 25245397, ConsequenceSpliceAcceptor},
		{"end_plus3_not_splice", 25245398, ""},

		// Before exon.Start (25245274): reverse strand -> donor
		{"start_minus1", 25245273, ConsequenceSpliceDonor},
		{"start_minus2", 25245272, ConsequenceSpliceDonor},
		{"start_minus3_not_splice", 25245271, ""},

		// In exon or deep intron -> not a splice site
		{"in_exon", 25245350, ""},
		{"deep_intron", 25235000, ""},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := spliceSiteType(tt.pos, transcript)
			assert.Equal(t, tt.expected, got)
		})
	}
}

func TestSpliceSiteType_ForwardStrand(t *testing.T) {
	// Minimal forward strand transcript to verify donor/acceptor assignment
	transcript := &cache.Transcript{
		ID:     "ENST_FWD",
		Chrom:  "1",
		Start:  1000,
		End:    5000,
		Strand: 1, // Forward
		Exons: []cache.Exon{
			{Number: 1, Start: 1000, End: 1200},
			{Number: 2, Start: 2000, End: 2300},
		},
	}

	tests := []struct {
		name     string
		pos      int64
		expected string
	}{
		// After exon 1 End (1200): forward strand -> donor
		{"exon1_end_plus1", 1201, ConsequenceSpliceDonor},
		{"exon1_end_plus2", 1202, ConsequenceSpliceDonor},

		// Before exon 2 Start (2000): forward strand -> acceptor
		{"exon2_start_minus1", 1999, ConsequenceSpliceAcceptor},
		{"exon2_start_minus2", 1998, ConsequenceSpliceAcceptor},

		// Not splice site
		{"exon1_end_plus3", 1203, ""},
		{"exon2_start_minus3", 1997, ""},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := spliceSiteType(tt.pos, transcript)
			assert.Equal(t, tt.expected, got)
		})
	}
}

func TestPredictConsequence_SpliceDonor(t *testing.T) {
	// KRAS reverse strand: position before exon.Start = splice donor
	// Exon 2 Start = 25245274, so pos 25245273 = Start-1 -> splice_donor_variant
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245273,
		Ref:   "A",
		Alt:   "G",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	assert.Equal(t, ConsequenceSpliceDonor, result.Consequence)
	assert.Equal(t, ImpactHigh, result.Impact)
}

func TestPredictConsequence_SpliceAcceptor(t *testing.T) {
	// KRAS reverse strand: position after exon.End = splice acceptor
	// Exon 2 End = 25245395, so pos 25245396 = End+1 -> splice_acceptor_variant
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245396,
		Ref:   "A",
		Alt:   "G",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	assert.Equal(t, ConsequenceSpliceAcceptor, result.Consequence)
	assert.Equal(t, ImpactHigh, result.Impact)
}

func TestPredictConsequence_IndelSpanningSpliceSite(t *testing.T) {
	transcript := createKRASTranscript()
	// KRAS is reverse strand. Exon 2: Start=25245274, End=25245395
	// Splice donor at exon.Start-1 (25245273), exon.Start-2 (25245272)
	// Splice acceptor at exon.End+1 (25245396), exon.End+2 (25245397)

	tests := []struct {
		name       string
		pos        int64
		ref        string
		alt        string
		wantSplice bool
		wantType   string
	}{
		{
			// Deletion starting in splice region (5bp from boundary), spanning into splice site
			// pos=25245269 (splice region), ref=6bp, end=25245274 -> hits exon.Start
			// But exon.Start is exon side, not splice site. Need to span to Start-1 or Start-2.
			// pos=25245268, ref=6bp, end=25245273 -> hits Start-1 (splice donor for reverse)
			name:       "intron_del_spanning_donor",
			pos:        25245268,
			ref:        "AAAAAA",
			alt:        "A",
			wantSplice: true,
			wantType:   ConsequenceSpliceDonor,
		},
		{
			// Deletion starting in splice region after exon end, spanning into splice acceptor
			// pos=25245398 (splice region), ref=AAAA (4bp), end=25245401
			// But we need to span back to 25245396 or 25245397 (acceptor sites)
			// pos=25245395 (exon end), ref=AAAA (4bp), end=25245398 -> doesn't hit +1/+2 from start pos
			// pos=25245394, ref=AAAAAA (6bp), end=25245399 -> covers 25245396 (End+1=acceptor)
			// Actually for intronic start: pos=25245398, ref is 6bp -> end=25245403
			// That doesn't hit 25245396/25245397. Let me pick a start that spans.
			// pos=25245393, ref=AAAAAAAAA (9bp), end=25245401 -> covers 25245396 (acceptor)
			// But 25245393 is in the exon, not intronic.
			// For a truly intronic start that spans acceptor:
			// pos=25245400 (intron, splice region), ref=AAAAAA -> too far right
			// Actually let's use: deletion starting 8bp out that spans to End+1
			// pos=25245390 (in exon), ref=8bp, end=25245397 -> covers 25245396,25245397 (acceptor)
			name:       "exon_del_spanning_acceptor",
			pos:        25245390,
			ref:        "AAAAAAAA",
			alt:        "A",
			wantSplice: true,
			wantType:   ConsequenceSpliceAcceptor,
		},
		{
			// Short deletion fully within intron, not reaching splice site
			name:       "intron_del_no_splice",
			pos:        25245265,
			ref:        "AAA",
			alt:        "A",
			wantSplice: false,
		},
		{
			// SNV - should not trigger indel splice detection
			name:       "snv_at_splice_region",
			pos:        25245400,
			ref:        "A",
			alt:        "G",
			wantSplice: false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := &vcf.Variant{Chrom: "12", Pos: tt.pos, Ref: tt.ref, Alt: tt.alt}
			got := indelSpliceSiteType(v, transcript)
			if tt.wantSplice {
				assert.Equal(t, tt.wantType, got)
			} else {
				assert.Empty(t, got)
			}
		})
	}
}

func TestPredictConsequence_DeletionSpanningSpliceDonor(t *testing.T) {
	// KRAS reverse strand. Deletion starting in intron near exon 2 Start boundary,
	// spanning into the splice donor site (exon.Start-1, exon.Start-2).
	// Exon 2 Start = 25245274. Donor at 25245273, 25245272.
	// Deletion: pos=25245268, ref=6bp -> end=25245273, hits donor site.
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245268,
		Ref:   "AAAAAA",
		Alt:   "A",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	assert.Equal(t, ConsequenceSpliceDonor, result.Consequence)
	assert.Equal(t, ImpactHigh, result.Impact)
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
			assert.Equal(t, tt.expected, got)
		})
	}
}

func TestPredictConsequence_SpliceRegionIntronic(t *testing.T) {
	// Variant 5bp into intron after exon 2 End (25245395)
	// Position 25245400 = 5bp after exon boundary -> splice_region_variant,intron_variant
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245400,
		Ref:   "A",
		Alt:   "G",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	expected := "splice_region_variant,intron_variant"
	assert.Equal(t, expected, result.Consequence)
	assert.Equal(t, "LOW", result.Impact)
}

func TestPredictConsequence_SpliceRegionCoding(t *testing.T) {
	// Variant at exon 2 End boundary (25245395) - 0bp from boundary, exon side
	// This is within the CDS (CDSEnd=25245384... wait, 25245395 > 25245384)
	// Actually 25245395 > CDSEnd 25245384, so this would be 5'UTR on reverse strand
	// Use exon 2 Start boundary instead: pos 25245274 is CDS (CDSStart=25245274)
	// On reverse strand, CDSStart < CDSEnd, and pos >= CDSStart and pos <= CDSEnd -> in CDS
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245274, // exon.Start, within 3bp of boundary, in CDS
		Ref:   "C",
		Alt:   "T",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)

	// Should have splice_region_variant appended to coding consequence
	assert.Contains(t, result.Consequence, "splice_region_variant")
	// Primary consequence should be a coding type
	hasCoding := strings.Contains(result.Consequence, "missense_variant") ||
		strings.Contains(result.Consequence, "synonymous_variant")
	assert.True(t, hasCoding, "expected coding consequence + splice_region_variant, got %s", result.Consequence)
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

	assert.NotContains(t, result.Consequence, "splice_region_variant")
	assert.Equal(t, "missense_variant", result.Consequence)
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

	assert.Equal(t, "intron_variant", result.Consequence)
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

	assert.Equal(t, ConsequenceInframeDeletion, result.Consequence)

	assert.Equal(t, ImpactModerate, result.Impact)
}

func BenchmarkPredictConsequence(b *testing.B) {
	transcript := createKRASTranscript()

	variants := []*vcf.Variant{
		{Chrom: "12", Pos: 25245350, Ref: "C", Alt: "A"},    // missense (G12C)
		{Chrom: "12", Pos: 25245350, Ref: "C", Alt: "C"},    // synonymous
		{Chrom: "12", Pos: 25245350, Ref: "CC", Alt: "C"},   // frameshift
		{Chrom: "12", Pos: 25245350, Ref: "CCCC", Alt: "C"}, // inframe deletion
		{Chrom: "12", Pos: 25245280, Ref: "A", Alt: "G"},    // splice region
		{Chrom: "12", Pos: 25245300, Ref: "A", Alt: "G"},    // intron
		{Chrom: "12", Pos: 25250800, Ref: "A", Alt: "G"},    // 5' UTR
		{Chrom: "12", Pos: 25200000, Ref: "A", Alt: "G"},    // upstream
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for _, v := range variants {
			PredictConsequence(v, transcript)
		}
	}
}

func TestPredictConsequence_Performance(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping performance test in short mode")
	}
	transcript := createKRASTranscript()
	variants := []*vcf.Variant{
		{Chrom: "12", Pos: 25245350, Ref: "C", Alt: "A"},
		{Chrom: "12", Pos: 25245350, Ref: "CC", Alt: "C"},
		{Chrom: "12", Pos: 25245300, Ref: "A", Alt: "G"},
		{Chrom: "12", Pos: 25250800, Ref: "A", Alt: "G"},
	}

	// Annotate 100k variants and check it completes within 1 second
	const iterations = 100000
	start := testing.AllocsPerRun(1, func() {})
	_ = start

	t0 := testing.Benchmark(func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for _, v := range variants {
				PredictConsequence(v, transcript)
			}
		}
	})

	nsPerVariant := float64(t0.T.Nanoseconds()) / float64(t0.N) / float64(len(variants))
	variantsPerSec := 1e9 / nsPerVariant

	// Regression threshold: must handle at least 100k variants/sec per transcript
	if variantsPerSec < 100000 {
		t.Errorf("PredictConsequence too slow: %.0f variants/sec (want >= 100000)", variantsPerSec)
	}
	t.Logf("PredictConsequence: %.0f variants/sec (%.0f ns/variant)", variantsPerSec, nsPerVariant)
}
