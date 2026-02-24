package annotate

import (
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"

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

	assert.Equal(t, "c.34G>T", hgvsc)
}

func TestFormatHGVSc_KRAS_Intronic(t *testing.T) {
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245272, // 2bp before exon 2 start (intronic, on the 3' side for reverse strand)
		Ref:   "A",
		Alt:   "G",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("HGVSc for intronic position: %s (consequence: %s)", hgvsc, result.Consequence)
	assert.NotEmpty(t, hgvsc)
}

func TestFormatHGVSc_ForwardStrand_SNV(t *testing.T) {
	transcript := createForwardTranscript()

	v := &vcf.Variant{
		Chrom: "1",
		Pos:   1005,
		Ref:   "A",
		Alt:   "G",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("Forward strand SNV HGVSc: %s", hgvsc)
	assert.NotEmpty(t, hgvsc)
	if len(hgvsc) > 2 {
		assert.True(t, strings.HasPrefix(hgvsc, "c."), "expected HGVSc to start with 'c.', got %q", hgvsc)
	}
}

func TestFormatHGVSc_Deletion(t *testing.T) {
	transcript := createForwardTranscript()

	v := &vcf.Variant{
		Chrom: "1",
		Pos:   1005,
		Ref:   "AG",
		Alt:   "A",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("Deletion HGVSc: %s", hgvsc)
	assert.NotEmpty(t, hgvsc)
	assert.True(t, strings.HasSuffix(hgvsc, "del"), "expected HGVSc to end with 'del', got %q", hgvsc)
}

func TestFormatHGVSc_Insertion(t *testing.T) {
	transcript := createForwardTranscript()

	v := &vcf.Variant{
		Chrom: "1",
		Pos:   1005,
		Ref:   "A",
		Alt:   "AGG",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("Insertion HGVSc: %s", hgvsc)
	assert.NotEmpty(t, hgvsc)
}

func TestFormatHGVSc_UTR(t *testing.T) {
	transcript := createForwardTranscript()

	v := &vcf.Variant{
		Chrom: "1",
		Pos:   998, // In exon 1 but before CDSStart=1000
		Ref:   "A",
		Alt:   "G",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("5'UTR HGVSc: %s (consequence: %s)", hgvsc, result.Consequence)
	assert.NotEmpty(t, hgvsc)
	if len(hgvsc) > 2 {
		assert.True(t, strings.HasPrefix(hgvsc, "c.-"), "expected HGVSc to start with 'c.-', got %q", hgvsc)
	}
}

func TestFormatHGVSc_UpstreamDownstream(t *testing.T) {
	transcript := createForwardTranscript()

	v := &vcf.Variant{
		Chrom: "1",
		Pos:   500,
		Ref:   "A",
		Alt:   "G",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Empty(t, hgvsc)
}

func TestGenomicToHGVScPos_CDS(t *testing.T) {
	transcript := createForwardTranscript()

	pos := genomicToHGVScPos(1005, transcript)
	assert.Equal(t, "6", pos)
}

func TestGenomicToHGVScPos_FivePrimeUTR(t *testing.T) {
	transcript := createForwardTranscript()

	pos := genomicToHGVScPos(998, transcript)
	assert.Equal(t, "-2", pos)
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
