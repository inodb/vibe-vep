package vcf

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestParser_SingleVariant(t *testing.T) {
	// Find testdata directory
	testFile := findTestFile(t, "kras_g12c.vcf")

	parser, err := NewParser(testFile)
	if err != nil {
		t.Fatalf("Failed to create parser: %v", err)
	}
	defer parser.Close()

	// Read the first (and only) variant
	v, err := parser.Next()
	if err != nil {
		t.Fatalf("Failed to read variant: %v", err)
	}

	if v == nil {
		t.Fatal("Expected a variant, got nil")
	}

	// Verify KRAS G12C variant (c.34G>T p.G12C)
	// On reverse strand: coding G->T = genomic C->A
	if v.Chrom != "12" {
		t.Errorf("Expected chrom 12, got %s", v.Chrom)
	}
	if v.Pos != 25245351 {
		t.Errorf("Expected pos 25245351, got %d", v.Pos)
	}
	if v.Ref != "C" {
		t.Errorf("Expected ref C, got %s", v.Ref)
	}
	if v.Alt != "A" {
		t.Errorf("Expected alt A, got %s", v.Alt)
	}

	// Should be a SNV
	if !v.IsSNV() {
		t.Error("KRAS G12C should be classified as SNV")
	}

	// No more variants
	v2, err := parser.Next()
	if err != nil {
		t.Fatalf("Error checking for more variants: %v", err)
	}
	if v2 != nil {
		t.Error("Expected no more variants")
	}
}

func TestParser_MultipleVariants(t *testing.T) {
	testFile := findTestFile(t, "multi_variant.vcf")

	parser, err := NewParser(testFile)
	if err != nil {
		t.Fatalf("Failed to create parser: %v", err)
	}
	defer parser.Close()

	// Count variants
	count := 0
	for {
		v, err := parser.Next()
		if err != nil {
			t.Fatalf("Error reading variant: %v", err)
		}
		if v == nil {
			break
		}
		count++
	}

	if count != 5 {
		t.Errorf("Expected 5 variants, got %d", count)
	}
}

func TestParser_Header(t *testing.T) {
	testFile := findTestFile(t, "kras_g12c.vcf")

	parser, err := NewParser(testFile)
	if err != nil {
		t.Fatalf("Failed to create parser: %v", err)
	}
	defer parser.Close()

	header := parser.Header()
	if len(header) == 0 {
		t.Error("Expected header lines")
	}

	// Check for required header elements
	hasFileformat := false
	hasChromLine := false
	for _, line := range header {
		if line == "##fileformat=VCFv4.2" {
			hasFileformat = true
		}
		if line[:6] == "#CHROM" {
			hasChromLine = true
		}
	}

	if !hasFileformat {
		t.Error("Missing ##fileformat header")
	}
	if !hasChromLine {
		t.Error("Missing #CHROM header line")
	}
}

func TestSplitMultiAllelic(t *testing.T) {
	tests := []struct {
		name     string
		alt      string
		expected int
	}{
		{"single allele", "C", 1},
		{"two alleles", "C,T", 2},
		{"three alleles", "C,T,G", 3},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := &Variant{
				Chrom: "12",
				Pos:   100,
				Ref:   "A",
				Alt:   tt.alt,
			}

			variants := SplitMultiAllelic(v)
			if len(variants) != tt.expected {
				t.Errorf("Expected %d variants, got %d", tt.expected, len(variants))
			}

			// Each variant should have only one alt allele
			for _, split := range variants {
				if strings.Contains(split.Alt, ",") {
					t.Errorf("Split variant should not contain comma in alt: %s", split.Alt)
				}
			}
		})
	}
}

func TestParseError(t *testing.T) {
	err := &ParseError{
		Line:    42,
		Message: "expected 8 columns, found 7",
	}

	expected := "vcf parse error at line 42: expected 8 columns, found 7"
	if err.Error() != expected {
		t.Errorf("Error message mismatch: got %q, want %q", err.Error(), expected)
	}
}

// findTestFile locates a test file in the testdata directory.
func findTestFile(t *testing.T, name string) string {
	t.Helper()

	// Try different relative paths
	paths := []string{
		filepath.Join("testdata", name),
		filepath.Join("..", "..", "testdata", name),
	}

	for _, p := range paths {
		if _, err := os.Stat(p); err == nil {
			return p
		}
	}

	t.Fatalf("Test file not found: %s", name)
	return ""
}

