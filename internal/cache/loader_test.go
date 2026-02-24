package cache

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestLoader_LoadJSONFile(t *testing.T) {
	// Find testdata directory
	testCacheDir := findTestCacheDir(t)

	loader := NewLoader(testCacheDir, "homo_sapiens", "GRCh38")
	c := New()

	// Load chromosome 12
	err := loader.Load(c, "12")
	require.NoError(t, err)

	// Should have loaded KRAS transcripts
	require.NotZero(t, c.TranscriptCount(), "Expected at least one transcript")

	// Find KRAS canonical transcript
	transcript := c.GetTranscript("ENST00000311936")
	require.NotNil(t, transcript, "Expected to find ENST00000311936 (KRAS canonical)")

	// Verify transcript properties
	assert.Equal(t, "KRAS", transcript.GeneName)
	assert.Equal(t, int8(-1), transcript.Strand)
	assert.True(t, transcript.IsCanonical)
	assert.True(t, transcript.IsProteinCoding())
	assert.Len(t, transcript.Exons, 5)
}

func TestLoader_FindTranscripts(t *testing.T) {
	testCacheDir := findTestCacheDir(t)

	loader := NewLoader(testCacheDir, "homo_sapiens", "GRCh38")
	c := New()

	err := loader.Load(c, "12")
	require.NoError(t, err)

	// Position 25245351 should overlap KRAS transcripts
	transcripts := c.FindTranscripts("12", 25245351)
	require.NotEmpty(t, transcripts, "Expected to find transcripts at position 25245351")

	// Should find KRAS
	foundKRAS := false
	for _, tr := range transcripts {
		if tr.GeneName == "KRAS" {
			foundKRAS = true
			break
		}
	}

	assert.True(t, foundKRAS, "Expected to find KRAS transcript at position 25245351")
}

func TestLoader_FindTranscripts_Intergenic(t *testing.T) {
	testCacheDir := findTestCacheDir(t)

	loader := NewLoader(testCacheDir, "homo_sapiens", "GRCh38")
	c := New()

	err := loader.Load(c, "12")
	require.NoError(t, err)

	// Position 1000000 should not overlap any transcripts in our test data
	transcripts := c.FindTranscripts("12", 1000000)
	assert.Empty(t, transcripts)
}

func TestLoader_LoadNonexistentChromosome(t *testing.T) {
	testCacheDir := findTestCacheDir(t)

	loader := NewLoader(testCacheDir, "homo_sapiens", "GRCh38")
	c := New()

	// Loading a chromosome that doesn't exist should not error
	err := loader.Load(c, "99")
	assert.NoError(t, err)

	// Cache should be empty for that chromosome
	transcripts := c.FindTranscripts("99", 100000)
	assert.Empty(t, transcripts)
}

func TestLoader_TranscriptCDSSequence(t *testing.T) {
	testCacheDir := findTestCacheDir(t)

	loader := NewLoader(testCacheDir, "homo_sapiens", "GRCh38")
	c := New()

	err := loader.Load(c, "12")
	require.NoError(t, err)

	transcript := c.GetTranscript("ENST00000311936")
	require.NotNil(t, transcript, "Expected to find ENST00000311936")

	// Verify CDS sequence starts with ATG (start codon)
	require.True(t, len(transcript.CDSSequence) >= 3, "CDS sequence too short")

	startCodon := transcript.CDSSequence[:3]
	assert.Equal(t, "ATG", startCodon)

	// Verify codon 12 is GGT (Glycine)
	// CDS positions 34-36 correspond to codon 12
	require.True(t, len(transcript.CDSSequence) >= 36, "CDS sequence too short for codon 12")

	codon12 := transcript.CDSSequence[33:36] // 0-indexed: positions 33-35
	assert.Equal(t, "GGT", codon12)
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
