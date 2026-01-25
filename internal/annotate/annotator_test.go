package annotate

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// TestAnnotator_KRASG12C_Integration is an end-to-end integration test
// that annotates the KRAS G12C variant and verifies the complete output.
func TestAnnotator_KRASG12C_Integration(t *testing.T) {
	// Load test cache
	testCacheDir := findTestCacheDir(t)
	c := cache.New()
	loader := cache.NewLoader(testCacheDir, "homo_sapiens", "GRCh38")

	if err := loader.Load(c, "12"); err != nil {
		t.Fatalf("Failed to load cache: %v", err)
	}

	if c.TranscriptCount() == 0 {
		t.Fatal("No transcripts loaded")
	}

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
	if err != nil {
		t.Fatalf("Annotation failed: %v", err)
	}

	if len(annotations) == 0 {
		t.Fatal("Expected at least one annotation")
	}

	// Find the canonical transcript annotation
	var canonicalAnn *Annotation
	for _, a := range annotations {
		if a.IsCanonical {
			canonicalAnn = a
			break
		}
	}

	if canonicalAnn == nil {
		t.Fatal("Expected to find canonical transcript annotation")
	}

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
		if tt.got != tt.expected {
			t.Errorf("%s: got %v, expected %v", tt.name, tt.got, tt.expected)
		}
	}

	// Verify codon change contains expected codons
	if !strings.Contains(strings.ToUpper(canonicalAnn.CodonChange), "GGT") {
		t.Errorf("CodonChange should contain GGT (ref codon), got %s", canonicalAnn.CodonChange)
	}
	if !strings.Contains(strings.ToUpper(canonicalAnn.CodonChange), "TGT") {
		t.Errorf("CodonChange should contain TGT (alt codon), got %s", canonicalAnn.CodonChange)
	}
}

func TestAnnotator_IntergenicVariant(t *testing.T) {
	// Load test cache
	testCacheDir := findTestCacheDir(t)
	c := cache.New()
	loader := cache.NewLoader(testCacheDir, "homo_sapiens", "GRCh38")

	if err := loader.Load(c, "12"); err != nil {
		t.Fatalf("Failed to load cache: %v", err)
	}

	ann := NewAnnotator(c)

	// Variant at position not overlapping any transcript
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   1000000, // Far from KRAS
		Ref:   "A",
		Alt:   "G",
	}

	annotations, err := ann.Annotate(v)
	if err != nil {
		t.Fatalf("Annotation failed: %v", err)
	}

	if len(annotations) != 1 {
		t.Fatalf("Expected 1 annotation for intergenic, got %d", len(annotations))
	}

	if annotations[0].Consequence != ConsequenceIntergenicVariant {
		t.Errorf("Expected intergenic_variant, got %s", annotations[0].Consequence)
	}

	if annotations[0].Impact != ImpactModifier {
		t.Errorf("Expected MODIFIER impact, got %s", annotations[0].Impact)
	}
}

func TestAnnotator_CanonicalOnly(t *testing.T) {
	testCacheDir := findTestCacheDir(t)
	c := cache.New()
	loader := cache.NewLoader(testCacheDir, "homo_sapiens", "GRCh38")

	if err := loader.Load(c, "12"); err != nil {
		t.Fatalf("Failed to load cache: %v", err)
	}

	ann := NewAnnotator(c)
	ann.SetCanonicalOnly(true)

	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245351,
		Ref:   "C",
		Alt:   "A",
	}

	annotations, err := ann.Annotate(v)
	if err != nil {
		t.Fatalf("Annotation failed: %v", err)
	}

	// Should only have canonical transcript
	if len(annotations) != 1 {
		t.Errorf("Expected 1 annotation (canonical only), got %d", len(annotations))
	}

	if len(annotations) > 0 && !annotations[0].IsCanonical {
		t.Error("Expected only canonical transcript annotation")
	}
}

func TestAnnotator_AnnotateAll(t *testing.T) {
	testCacheDir := findTestCacheDir(t)
	c := cache.New()
	loader := cache.NewLoader(testCacheDir, "homo_sapiens", "GRCh38")

	if err := loader.Load(c, "12"); err != nil {
		t.Fatalf("Failed to load cache: %v", err)
	}

	ann := NewAnnotator(c)
	ann.SetCanonicalOnly(true)

	// Parse the test VCF
	vcfPath := findTestVCF(t, "kras_g12c.vcf")
	parser, err := vcf.NewParser(vcfPath)
	if err != nil {
		t.Fatalf("Failed to create parser: %v", err)
	}
	defer parser.Close()

	// Use a mock writer to capture output
	var buf bytes.Buffer
	writer := &mockWriter{buf: &buf}

	if err := ann.AnnotateAll(parser, writer); err != nil {
		t.Fatalf("AnnotateAll failed: %v", err)
	}

	// Verify we got KRAS annotation
	output := buf.String()
	if !strings.Contains(output, "KRAS") {
		t.Error("Expected output to contain KRAS")
	}
	if !strings.Contains(output, "missense_variant") {
		t.Error("Expected output to contain missense_variant")
	}
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
