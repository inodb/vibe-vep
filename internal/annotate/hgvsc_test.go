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

func TestFormatHGVSc_DupPrecedingBase(t *testing.T) {
	// Test: insertion that duplicates the preceding (anchor) base.
	// Forward strand transcript, CDS = "ATGATCGATCG..."
	// Insert at pos 1002 (CDS pos 3 = 'G'), ref=G, alt=GG
	// Inserted 'G' matches CDS pos 3 → c.3dup
	transcript := createDupTestTranscript()

	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "G", Alt: "GG"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.3dup", hgvsc)
}

func TestFormatHGVSc_DupFollowingBase(t *testing.T) {
	// Test: insertion that duplicates the following base (not the anchor).
	// Forward strand transcript, CDS = "ATGATCGATCG..."
	// CDS pos 6='C', pos 7='G'
	// Insert at pos 1005 (CDS pos 6 = 'C'), ref=C, alt=CG
	// Inserted 'G' doesn't match CDS pos 6 ('C'), but matches CDS pos 7 ('G') → c.7dup
	transcript := createDupTestTranscript()

	v := &vcf.Variant{Chrom: "1", Pos: 1005, Ref: "C", Alt: "CG"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.7dup", hgvsc)
}

func TestFormatHGVSc_DupMultiBase(t *testing.T) {
	// Test: multi-base duplication.
	// CDS = "ATGATCGATCG..."
	// CDS pos 1-3 = "ATG"
	// Insert at pos 1002 (CDS 3), ref=G, alt=GATG
	// Inserted "ATG" matches CDS pos 1-3 → c.1_3dup
	transcript := createDupTestTranscript()

	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "G", Alt: "GATG"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.1_3dup", hgvsc)
}

func TestFormatHGVSc_InsertionNotDup(t *testing.T) {
	// Test: insertion that does NOT match adjacent bases → plain insertion.
	// CDS = "ATGATCGATCG..."
	// Insert at pos 1002 (CDS 3 = 'G'), ref=G, alt=GCC
	// Inserted 'CC' doesn't match "TG" (preceding) or "AT" (following) → insertion
	transcript := createDupTestTranscript()

	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "G", Alt: "GCC"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Contains(t, hgvsc, "ins")
	assert.NotContains(t, hgvsc, "dup")
}

func TestFormatHGVSc_DupReverseStrand(t *testing.T) {
	// Test: duplication on reverse strand.
	// KRAS reverse strand. CDS[0]='A', CDS[1]='T', CDS[2]='G', CDS[3]='A', CDS[4]='C', CDS[5]='T'
	// CDS position 1 = genomic 25245384 (CDSEnd for reverse strand)
	// CDS position 34 = genomic 25245351 (KRAS G12 position)
	// CDS positions go: 25245384=1, 25245383=2, 25245382=3, ...
	// For a simple dup test, use a position where the base is known:
	// CDS[33] (pos 34) is first base of codon 12 = 'G' (GGT -> Gly)
	// CDS[34] (pos 35) is second base = 'G'
	// Insert 'G' at CDS pos 34: should dup pos 34
	// Genomic pos for CDS 34 = 25245384 - 34 + 1 = 25245351
	// On reverse strand, insert at 25245351: ref = genomic C (rev comp = G), alt = CC (rev comp = GG)
	transcript := createKRASTranscript()

	// VCF: pos=25245351, ref=C, alt=CC → on coding strand: ref=G, alt=GG → ins G
	// CDS pos 34 = 'G', so inserted G matches preceding base → c.34dup
	v := &vcf.Variant{Chrom: "12", Pos: 25245351, Ref: "C", Alt: "CC"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.34dup", hgvsc)
}

// createDupTestTranscript creates a forward-strand transcript with a known CDS for dup testing.
func createDupTestTranscript() *cache.Transcript {
	// Single exon forward strand transcript for simplicity.
	// Exon: 990-1200, CDS: 1000-1101 (102bp = 34 codons)
	// CDS: pos 1=1000, pos 2=1001, ..., pos 102=1101
	// CDSSequence[0..101] = "ATGATCGATCG..." (repeating pattern)
	cds := "ATGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
	return &cache.Transcript{
		ID:       "ENST00000DUP001",
		GeneName: "DUPTEST",
		Chrom:    "1",
		Start:    990,
		End:      1210,
		Strand:   1,
		Biotype:  "protein_coding",
		CDSStart: 1000,
		CDSEnd:   1101,
		Exons: []cache.Exon{
			{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0},
		},
		CDSSequence: cds,
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
