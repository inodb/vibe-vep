package annotate

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestAnnotator_KRASG12C_Integration is an end-to-end integration test
// that annotates the KRAS G12C variant and verifies the complete output.
func TestAnnotator_KRASG12C_Integration(t *testing.T) {
	// Load test cache
	testCacheDir := findTestCacheDir(t)
	c := cache.New()
	loader := cache.NewLoader(testCacheDir, "homo_sapiens", "GRCh38")

	require.NoError(t, loader.Load(c, "12"), "loading cache")
	require.NotZero(t, c.TranscriptCount(), "no transcripts loaded")

	// Create annotator
	ann := NewAnnotator(c)

	// Create KRAS G12C variant (c.34G>T p.G12C)
	// Genomic coordinates: chr12:25245351 C>A (reverse strand)
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245351,
		ID:    ".",
		Ref:   "C",
		Alt:   "A",
	}

	// Annotate
	annotations, err := ann.Annotate(v)
	require.NoError(t, err, "annotation failed")
	require.NotEmpty(t, annotations)

	// Find the canonical transcript annotation
	var canonicalAnn *Annotation
	for _, a := range annotations {
		if a.IsCanonical {
			canonicalAnn = a
			break
		}
	}

	require.NotNil(t, canonicalAnn, "expected canonical transcript annotation")

	// Verify annotation fields
	tests := []struct {
		name     string
		got      interface{}
		expected interface{}
	}{
		{"GeneName", canonicalAnn.GeneName, "KRAS"},
		{"TranscriptID", canonicalAnn.TranscriptID, "ENST00000311936"},
		{"Consequence", canonicalAnn.Consequence, ConsequenceMissenseVariant},
		{"Impact", canonicalAnn.Impact, ImpactModerate},
		{"CDSPosition", canonicalAnn.CDSPosition, int64(34)},
		{"ProteinPosition", canonicalAnn.ProteinPosition, int64(12)},
		{"AminoAcidChange", canonicalAnn.AminoAcidChange, "G12C"},
		{"IsCanonical", canonicalAnn.IsCanonical, true},
		{"Biotype", canonicalAnn.Biotype, "protein_coding"},
	}

	for _, tt := range tests {
		assert.Equal(t, tt.expected, tt.got, tt.name)
	}

	// Verify codon change contains expected codons
	assert.Contains(t, strings.ToUpper(canonicalAnn.CodonChange), "GGT")
	assert.Contains(t, strings.ToUpper(canonicalAnn.CodonChange), "TGT")
}

func TestAnnotator_IntergenicVariant(t *testing.T) {
	// Load test cache
	testCacheDir := findTestCacheDir(t)
	c := cache.New()
	loader := cache.NewLoader(testCacheDir, "homo_sapiens", "GRCh38")

	require.NoError(t, loader.Load(c, "12"), "loading cache")

	ann := NewAnnotator(c)

	// Variant at position not overlapping any transcript
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   1000000, // Far from KRAS
		Ref:   "A",
		Alt:   "G",
	}

	annotations, err := ann.Annotate(v)
	require.NoError(t, err, "annotation failed")
	require.Len(t, annotations, 1)

	assert.Equal(t, ConsequenceIntergenicVariant, annotations[0].Consequence)
	assert.Equal(t, ImpactModifier, annotations[0].Impact)
}

func TestAnnotator_CanonicalOnly(t *testing.T) {
	testCacheDir := findTestCacheDir(t)
	c := cache.New()
	loader := cache.NewLoader(testCacheDir, "homo_sapiens", "GRCh38")

	require.NoError(t, loader.Load(c, "12"), "loading cache")

	ann := NewAnnotator(c)
	ann.SetCanonicalOnly(true)

	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245351,
		Ref:   "C",
		Alt:   "A",
	}

	annotations, err := ann.Annotate(v)
	require.NoError(t, err, "annotation failed")

	// Should only have canonical transcript
	assert.Len(t, annotations, 1)

	if len(annotations) > 0 {
		assert.True(t, annotations[0].IsCanonical)
	}
}

func TestAnnotator_AnnotateAll(t *testing.T) {
	testCacheDir := findTestCacheDir(t)
	c := cache.New()
	loader := cache.NewLoader(testCacheDir, "homo_sapiens", "GRCh38")

	require.NoError(t, loader.Load(c, "12"), "loading cache")

	ann := NewAnnotator(c)
	ann.SetCanonicalOnly(true)

	// Parse the test VCF
	vcfPath := findTestVCF(t, "kras_g12c.vcf")
	parser, err := vcf.NewParser(vcfPath)
	require.NoError(t, err, "creating parser")
	defer parser.Close()

	// Use a mock writer to capture output
	var buf bytes.Buffer
	writer := &mockWriter{buf: &buf}

	require.NoError(t, ann.AnnotateAll(parser, writer), "AnnotateAll failed")

	// Verify we got KRAS annotation
	output := buf.String()
	assert.Contains(t, output, "KRAS")
	assert.Contains(t, output, "missense_variant")
}

// mockWriter implements AnnotationWriter for testing.
type mockWriter struct {
	buf *bytes.Buffer
}

func (w *mockWriter) WriteHeader() error {
	return nil
}

func (w *mockWriter) Write(v *vcf.Variant, ann *Annotation) error {
	w.buf.WriteString(ann.GeneName + "\t" + ann.Consequence + "\n")
	return nil
}

func (w *mockWriter) Flush() error {
	return nil
}

// TestAnnotator_VEP_MRPL39 validates consequence predictions against the
// ensembl-vep test suite expectations for the MRPL39 gene on chr21.
// VEP test variant rs142513484 (chr21:25585733 C>T) expects:
//   - ENST00000352957: missense_variant (A/T, codons Gca/Aca, CDS pos 991)
//   - ENST00000307301: 3_prime_UTR_variant
func TestAnnotator_VEP_MRPL39(t *testing.T) {
	testCacheDir := findTestCacheDir(t)
	c := cache.New()
	loader := cache.NewLoader(testCacheDir, "homo_sapiens", "GRCh38")

	require.NoError(t, loader.Load(c, "21"), "loading chr21 cache")
	require.NotZero(t, c.TranscriptCount(), "no transcripts loaded")

	ann := NewAnnotator(c)

	tests := []struct {
		name         string
		variant      *vcf.Variant
		transcriptID string
		wantConseq   string
		wantImpact   string
		wantAA       string // amino acid change (empty = don't check)
	}{
		{
			name:         "rs142513484 canonical missense",
			variant:      &vcf.Variant{Chrom: "21", Pos: 25585733, Ref: "C", Alt: "T"},
			transcriptID: "ENST00000352957",
			wantConseq:   ConsequenceMissenseVariant,
			wantImpact:   ImpactModerate,
			wantAA:       "A331T",
		},
		{
			name:         "rs142513484 3'UTR on ENST00000307301",
			variant:      &vcf.Variant{Chrom: "21", Pos: 25585733, Ref: "C", Alt: "T"},
			transcriptID: "ENST00000307301",
			wantConseq:   Consequence3PrimeUTR,
			wantImpact:   ImpactModifier,
		},
		{
			name:         "rs199510789 synonymous shared",
			variant:      &vcf.Variant{Chrom: "21", Pos: 25588859, Ref: "C", Alt: "T"},
			transcriptID: "ENST00000352957",
			wantConseq:   ConsequenceSynonymousVariant,
			wantImpact:   ImpactLow,
		},
		{
			name:         "rs115023280 missense ENST00000352957",
			variant:      &vcf.Variant{Chrom: "21", Pos: 25592826, Ref: "G", Alt: "A"},
			transcriptID: "ENST00000352957",
			wantConseq:   ConsequenceMissenseVariant,
			wantImpact:   ImpactModerate,
		},
		{
			name:         "rs1135638 synonymous",
			variant:      &vcf.Variant{Chrom: "21", Pos: 25592836, Ref: "G", Alt: "A"},
			transcriptID: "ENST00000352957",
			wantConseq:   ConsequenceSynonymousVariant,
			wantImpact:   ImpactLow,
		},
		{
			name:         "rs187353664 intron in canonical",
			variant:      &vcf.Variant{Chrom: "21", Pos: 25587701, Ref: "T", Alt: "C"},
			transcriptID: "ENST00000352957",
			wantConseq:   ConsequenceIntronVariant,
			wantImpact:   ImpactModifier,
		},
		{
			name:         "rs7278168 synonymous",
			variant:      &vcf.Variant{Chrom: "21", Pos: 25603910, Ref: "C", Alt: "T"},
			transcriptID: "ENST00000352957",
			wantConseq:   ConsequenceSynonymousVariant,
			wantImpact:   ImpactLow,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			annotations, err := ann.Annotate(tt.variant)
			require.NoError(t, err, "annotation failed")
			require.NotEmpty(t, annotations, "no annotations returned")

			// Find annotation for expected transcript
			var found *Annotation
			for _, a := range annotations {
				if strings.HasPrefix(a.TranscriptID, tt.transcriptID) {
					found = a
					break
				}
			}
			require.NotNilf(t, found, "no annotation for transcript %s", tt.transcriptID)

			assert.Equal(t, tt.wantConseq, strings.Split(found.Consequence, ",")[0],
				"consequence mismatch")
			assert.Equal(t, tt.wantImpact, found.Impact, "impact mismatch")

			if tt.wantAA != "" {
				assert.Equal(t, tt.wantAA, found.AminoAcidChange, "amino acid change")
			}
		})
	}
}

// findTestCacheDir locates the test cache directory.
func findTestCacheDir(t *testing.T) string {
	t.Helper()

	paths := []string{
		filepath.Join("testdata", "cache"),
		filepath.Join("..", "..", "testdata", "cache"),
	}

	for _, p := range paths {
		if _, err := os.Stat(p); err == nil {
			return p
		}
	}

	t.Fatal("Test cache directory not found")
	return ""
}

// findTestVCF locates a test VCF file.
func findTestVCF(t *testing.T, name string) string {
	t.Helper()

	paths := []string{
		filepath.Join("testdata", name),
		filepath.Join("..", "..", "testdata", name),
	}

	for _, p := range paths {
		if _, err := os.Stat(p); err == nil {
			return p
		}
	}

	t.Fatalf("Test VCF file not found: %s", name)
	return ""
}
