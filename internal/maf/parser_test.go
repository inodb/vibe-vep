package maf

import (
	"os"
	"path/filepath"
	"testing"
)

func TestParser_ParseVariants(t *testing.T) {
	testFile := findTestFile(t, "sample.maf")

	parser, err := NewParser(testFile)
	if err != nil {
		t.Fatalf("Failed to create parser: %v", err)
	}
	defer parser.Close()

	// Verify column indices were parsed correctly
	cols := parser.Columns()
	if cols.Chromosome != 4 {
		t.Errorf("Expected Chromosome column at index 4, got %d", cols.Chromosome)
	}
	if cols.StartPosition != 5 {
		t.Errorf("Expected Start_Position column at index 5, got %d", cols.StartPosition)
	}
	if cols.ReferenceAllele != 11 {
		t.Errorf("Expected Reference_Allele column at index 11, got %d", cols.ReferenceAllele)
	}
	if cols.TumorSeqAllele2 != 13 {
		t.Errorf("Expected Tumor_Seq_Allele2 column at index 13, got %d", cols.TumorSeqAllele2)
	}

	// Read first variant (TRUB1)
	v, err := parser.Next()
	if err != nil {
		t.Fatalf("Failed to read variant: %v", err)
	}
	if v == nil {
		t.Fatal("Expected a variant, got nil")
	}

	if v.Chrom != "10" {
		t.Errorf("Expected chrom 10, got %s", v.Chrom)
	}
	if v.Pos != 116734973 {
		t.Errorf("Expected pos 116734973, got %d", v.Pos)
	}
	if v.Ref != "G" {
		t.Errorf("Expected ref G, got %s", v.Ref)
	}
	if v.Alt != "A" {
		t.Errorf("Expected alt A, got %s", v.Alt)
	}

	// Read second variant (RRM1)
	v, err = parser.Next()
	if err != nil {
		t.Fatalf("Failed to read variant: %v", err)
	}
	if v == nil {
		t.Fatal("Expected a variant, got nil")
	}

	if v.Chrom != "11" {
		t.Errorf("Expected chrom 11, got %s", v.Chrom)
	}
	if v.Pos != 4148284 {
		t.Errorf("Expected pos 4148284, got %d", v.Pos)
	}

	// Read third variant (KRAS G12C)
	v, err = parser.Next()
	if err != nil {
		t.Fatalf("Failed to read variant: %v", err)
	}
	if v == nil {
		t.Fatal("Expected a variant, got nil")
	}

	if v.Chrom != "12" {
		t.Errorf("Expected chrom 12, got %s", v.Chrom)
	}
	if v.Pos != 25398285 {
		t.Errorf("Expected pos 25398285, got %d", v.Pos)
	}
	if v.Ref != "G" {
		t.Errorf("Expected ref G, got %s", v.Ref)
	}
	if v.Alt != "T" {
		t.Errorf("Expected alt T, got %s", v.Alt)
	}

	// Count remaining variants
	count := 3 // Already read 3
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

func TestParser_WithAnnotation(t *testing.T) {
	testFile := findTestFile(t, "sample.maf")

	parser, err := NewParser(testFile)
	if err != nil {
		t.Fatalf("Failed to create parser: %v", err)
	}
	defer parser.Close()

	// Read first variant with annotation
	v, ann, err := parser.NextWithAnnotation()
	if err != nil {
		t.Fatalf("Failed to read variant: %v", err)
	}
	if v == nil || ann == nil {
		t.Fatal("Expected variant and annotation, got nil")
	}

	// Verify annotation data
	if ann.HugoSymbol != "TRUB1" {
		t.Errorf("Expected HugoSymbol TRUB1, got %s", ann.HugoSymbol)
	}
	if ann.Consequence != "stop_gained" {
		t.Errorf("Expected Consequence stop_gained, got %s", ann.Consequence)
	}
	if ann.HGVSpShort != "p.W295*" {
		t.Errorf("Expected HGVSpShort p.W295*, got %s", ann.HGVSpShort)
	}
	if ann.TranscriptID != "ENST00000298746" {
		t.Errorf("Expected TranscriptID ENST00000298746, got %s", ann.TranscriptID)
	}
	if ann.NCBIBuild != "GRCh37" {
		t.Errorf("Expected NCBIBuild GRCh37, got %s", ann.NCBIBuild)
	}

	// Read KRAS G12C variant
	parser.NextWithAnnotation() // skip RRM1
	v, ann, err = parser.NextWithAnnotation()
	if err != nil {
		t.Fatalf("Failed to read variant: %v", err)
	}

	if ann.HugoSymbol != "KRAS" {
		t.Errorf("Expected HugoSymbol KRAS, got %s", ann.HugoSymbol)
	}
	if ann.HGVSpShort != "p.G12C" {
		t.Errorf("Expected HGVSpShort p.G12C, got %s", ann.HGVSpShort)
	}
}

func TestParser_Header(t *testing.T) {
	testFile := findTestFile(t, "sample.maf")

	parser, err := NewParser(testFile)
	if err != nil {
		t.Fatalf("Failed to create parser: %v", err)
	}
	defer parser.Close()

	header := parser.Header()
	if header == "" {
		t.Error("Expected header line")
	}

	// Check for expected columns in header
	if len(header) < 50 {
		t.Error("Header seems too short")
	}
}

func TestParseError(t *testing.T) {
	err := &ParseError{
		Line:    42,
		Message: "required column not found",
	}

	expected := "maf parse error at line 42: required column not found"
	if err.Error() != expected {
		t.Errorf("Error message mismatch: got %q, want %q", err.Error(), expected)
	}
}

func TestParser_ImplementsVariantParser(t *testing.T) {
	testFile := findTestFile(t, "sample.maf")

	parser, err := NewParser(testFile)
	if err != nil {
		t.Fatalf("Failed to create parser: %v", err)
	}
	defer parser.Close()

	// Verify the parser has the expected interface methods
	_ = parser.LineNumber()
	v, err := parser.Next()
	if err != nil {
		t.Fatalf("Next() failed: %v", err)
	}
	if v == nil {
		t.Fatal("Expected variant from Next()")
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
