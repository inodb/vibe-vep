package cache

import (
	"os"
	"path/filepath"
	"testing"
)

func TestLoader_LoadJSONFile(t *testing.T) {
	// Find testdata directory
	testCacheDir := findTestCacheDir(t)

	loader := NewLoader(testCacheDir, "homo_sapiens", "GRCh38")
	c := New()

	// Load chromosome 12
	err := loader.Load(c, "12")
	if err != nil {
		t.Fatalf("Failed to load chromosome 12: %v", err)
	}

	// Should have loaded KRAS transcripts
	if c.TranscriptCount() == 0 {
		t.Fatal("Expected at least one transcript, got 0")
	}

	// Find KRAS canonical transcript
	transcript := c.GetTranscript("ENST00000311936")
	if transcript == nil {
		t.Fatal("Expected to find ENST00000311936 (KRAS canonical)")
	}

	// Verify transcript properties
	if transcript.GeneName != "KRAS" {
		t.Errorf("Expected gene name KRAS, got %s", transcript.GeneName)
	}

	if transcript.Strand != -1 {
		t.Errorf("Expected strand -1 (reverse), got %d", transcript.Strand)
	}

	if !transcript.IsCanonical {
		t.Error("Expected ENST00000311936 to be canonical")
	}

	if !transcript.IsProteinCoding() {
		t.Error("Expected KRAS to be protein coding")
	}

	if len(transcript.Exons) != 5 {
		t.Errorf("Expected 5 exons, got %d", len(transcript.Exons))
	}
}

func TestLoader_FindTranscripts(t *testing.T) {
	testCacheDir := findTestCacheDir(t)

	loader := NewLoader(testCacheDir, "homo_sapiens", "GRCh38")
	c := New()

	err := loader.Load(c, "12")
	if err != nil {
		t.Fatalf("Failed to load chromosome 12: %v", err)
	}

	// Position 25245351 should overlap KRAS transcripts
	transcripts := c.FindTranscripts("12", 25245351)
	if len(transcripts) == 0 {
		t.Fatal("Expected to find transcripts at position 25245351")
	}

	// Should find KRAS
	foundKRAS := false
	for _, tr := range transcripts {
		if tr.GeneName == "KRAS" {
			foundKRAS = true
			break
		}
	}

	if !foundKRAS {
		t.Error("Expected to find KRAS transcript at position 25245351")
	}
}

func TestLoader_FindTranscripts_Intergenic(t *testing.T) {
	testCacheDir := findTestCacheDir(t)

	loader := NewLoader(testCacheDir, "homo_sapiens", "GRCh38")
	c := New()

	err := loader.Load(c, "12")
	if err != nil {
		t.Fatalf("Failed to load chromosome 12: %v", err)
	}

	// Position 1000000 should not overlap any transcripts in our test data
	transcripts := c.FindTranscripts("12", 1000000)
	if len(transcripts) != 0 {
		t.Errorf("Expected no transcripts at position 1000000, got %d", len(transcripts))
	}
}

func TestLoader_LoadNonexistentChromosome(t *testing.T) {
	testCacheDir := findTestCacheDir(t)

	loader := NewLoader(testCacheDir, "homo_sapiens", "GRCh38")
	c := New()

	// Loading a chromosome that doesn't exist should not error
	err := loader.Load(c, "99")
	if err != nil {
		t.Errorf("Expected no error for nonexistent chromosome, got: %v", err)
	}

	// Cache should be empty for that chromosome
	transcripts := c.FindTranscripts("99", 100000)
	if len(transcripts) != 0 {
		t.Errorf("Expected no transcripts for chromosome 99, got %d", len(transcripts))
	}
}

func TestLoader_TranscriptCDSSequence(t *testing.T) {
	testCacheDir := findTestCacheDir(t)

	loader := NewLoader(testCacheDir, "homo_sapiens", "GRCh38")
	c := New()

	err := loader.Load(c, "12")
	if err != nil {
		t.Fatalf("Failed to load chromosome 12: %v", err)
	}

	transcript := c.GetTranscript("ENST00000311936")
	if transcript == nil {
		t.Fatal("Expected to find ENST00000311936")
	}

	// Verify CDS sequence starts with ATG (start codon)
	if len(transcript.CDSSequence) < 3 {
		t.Fatal("CDS sequence too short")
	}

	startCodon := transcript.CDSSequence[:3]
	if startCodon != "ATG" {
		t.Errorf("Expected CDS to start with ATG, got %s", startCodon)
	}

	// Verify codon 12 is GGT (Glycine)
	// CDS positions 34-36 correspond to codon 12
	if len(transcript.CDSSequence) < 36 {
		t.Fatal("CDS sequence too short for codon 12")
	}

	codon12 := transcript.CDSSequence[33:36] // 0-indexed: positions 33-35
	if codon12 != "GGT" {
		t.Errorf("Expected codon 12 to be GGT, got %s", codon12)
	}
}

// findTestCacheDir locates the test cache directory.
func findTestCacheDir(t *testing.T) string {
	t.Helper()

	// Try different relative paths
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
