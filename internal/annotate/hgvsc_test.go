package annotate

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

func TestFormatHGVSc_KRAS_SNV(t *testing.T) {
	// KRAS G12C: c.34G>T on reverse strand
	// Genomic pos 25245351, ref=C (genomic), alt=A (genomic)
	// On coding strand: ref=G, alt=T → c.34G>T
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245351,
		Ref:   "C",
		Alt:   "A",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	if hgvsc != "c.34G>T" {
		t.Errorf("Expected c.34G>T, got %q", hgvsc)
	}
}

func TestFormatHGVSc_KRAS_Intronic(t *testing.T) {
	// Test an intronic position between exon 2 and exon 3 of KRAS (reverse strand)
	// Exon 2: 25245274-25245395 (CDS: 25245274-25245384)
	// Exon 3: 25227234-25227412
	// An intronic position just after exon 2 (downstream in genomic = upstream in transcript)
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245272, // 2bp before exon 2 start (intronic, on the 3' side for reverse strand)
		Ref:   "A",
		Alt:   "G",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	// For reverse strand, pos 25245272 is 2bp before exon 2 start (25245274)
	// In transcript terms, this is after exon 2's last CDS base
	// Exon 2 CDS ends at position 25245274 (the last coding base at the 3' end genomically)
	// GenomicToCDS(25245274) gives the CDS position of that boundary
	// The offset is 25245274 - 25245272 = 2
	// Since reverse strand and pos < exon.Start, this is "after" the exon in transcript terms
	// → c.{CDSpos}+{offset}
	t.Logf("HGVSc for intronic position: %s (consequence: %s)", hgvsc, result.Consequence)
	if hgvsc == "" {
		t.Error("Expected non-empty HGVSc for intronic variant")
	}
}

func TestFormatHGVSc_ForwardStrand_SNV(t *testing.T) {
	// Create a simple forward-strand transcript for testing
	transcript := createForwardTranscript()

	// SNV in CDS at position 1005 → CDS pos depends on exon layout
	v := &vcf.Variant{
		Chrom: "1",
		Pos:   1005,
		Ref:   "A",
		Alt:   "G",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("Forward strand SNV HGVSc: %s", hgvsc)
	if hgvsc == "" {
		t.Error("Expected non-empty HGVSc for coding SNV")
	}
	// Should be c.{pos}A>G (forward strand, no complement needed)
	if len(hgvsc) > 2 && hgvsc[:2] != "c." {
		t.Errorf("Expected HGVSc to start with 'c.', got %q", hgvsc)
	}
}

func TestFormatHGVSc_Deletion(t *testing.T) {
	transcript := createForwardTranscript()

	// Single base deletion in CDS
	v := &vcf.Variant{
		Chrom: "1",
		Pos:   1005,
		Ref:   "AG",
		Alt:   "A",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("Deletion HGVSc: %s", hgvsc)
	if hgvsc == "" {
		t.Error("Expected non-empty HGVSc for deletion")
	}
	if len(hgvsc) < 4 || hgvsc[len(hgvsc)-3:] != "del" {
		t.Errorf("Expected HGVSc to end with 'del', got %q", hgvsc)
	}
}

func TestFormatHGVSc_Insertion(t *testing.T) {
	transcript := createForwardTranscript()

	// Insertion in CDS
	v := &vcf.Variant{
		Chrom: "1",
		Pos:   1005,
		Ref:   "A",
		Alt:   "AGG",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("Insertion HGVSc: %s", hgvsc)
	if hgvsc == "" {
		t.Error("Expected non-empty HGVSc for insertion")
	}
}

func TestFormatHGVSc_UTR(t *testing.T) {
	transcript := createForwardTranscript()

	// 5'UTR variant (before CDS start)
	v := &vcf.Variant{
		Chrom: "1",
		Pos:   998, // In exon 1 but before CDSStart=1000
		Ref:   "A",
		Alt:   "G",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("5'UTR HGVSc: %s (consequence: %s)", hgvsc, result.Consequence)
	if hgvsc == "" {
		t.Error("Expected non-empty HGVSc for 5'UTR variant")
	}
	// Should start with "c.-"
	if len(hgvsc) > 2 && hgvsc[:3] != "c.-" {
		t.Errorf("Expected HGVSc to start with 'c.-', got %q", hgvsc)
	}
}

func TestFormatHGVSc_UpstreamDownstream(t *testing.T) {
	transcript := createForwardTranscript()

	// Upstream variant - should return empty
	v := &vcf.Variant{
		Chrom: "1",
		Pos:   500,
		Ref:   "A",
		Alt:   "G",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	if hgvsc != "" {
		t.Errorf("Expected empty HGVSc for upstream variant, got %q", hgvsc)
	}
}

func TestGenomicToHGVScPos_CDS(t *testing.T) {
	transcript := createForwardTranscript()

	// Position in CDS
	pos := genomicToHGVScPos(1005, transcript)
	if pos != "6" {
		t.Errorf("Expected CDS position '6', got %q", pos)
	}
}

func TestGenomicToHGVScPos_FivePrimeUTR(t *testing.T) {
	transcript := createForwardTranscript()

	// Position in 5'UTR
	pos := genomicToHGVScPos(998, transcript)
	if pos != "-2" {
		t.Errorf("Expected 5'UTR position '-2', got %q", pos)
	}
}

// createForwardTranscript creates a simple forward-strand transcript for testing.
func createForwardTranscript() *cache.Transcript {
	// Simple 2-exon forward strand transcript
	// Exon 1: 990-1020 (CDS starts at 1000)
	// Intron: 1021-1099
	// Exon 2: 1100-1200 (CDS ends at 1180)
	//
	// CDS: 1000-1020 (21bp from exon 1) + 1100-1180 (81bp from exon 2) = 102bp = 34 codons
	return &cache.Transcript{
		ID:          "ENST00000000001",
		GeneID:      "ENSG00000000001",
		GeneName:    "TEST",
		Chrom:       "1",
		Start:       990,
		End:         1200,
		Strand:      1,
		Biotype:     "protein_coding",
		IsCanonical: true,
		CDSStart:    1000,
		CDSEnd:      1180,
		Exons: []cache.Exon{
			{Number: 1, Start: 990, End: 1020, CDSStart: 1000, CDSEnd: 1020, Frame: 0},
			{Number: 2, Start: 1100, End: 1200, CDSStart: 1100, CDSEnd: 1180, Frame: 0},
		},
		CDSSequence: "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG",
	}
}
