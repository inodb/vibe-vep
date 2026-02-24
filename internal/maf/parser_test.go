package maf

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestParser_ParseVariants(t *testing.T) {
	testFile := findTestFile(t, "sample.maf")

	parser, err := NewParser(testFile)
	require.NoError(t, err)
	defer parser.Close()

	// Verify column indices were parsed correctly
	cols := parser.Columns()
	assert.Equal(t, 4, cols.Chromosome)
	assert.Equal(t, 5, cols.StartPosition)
	assert.Equal(t, 11, cols.ReferenceAllele)
	assert.Equal(t, 13, cols.TumorSeqAllele2)

	// Read first variant (TRUB1)
	v, err := parser.Next()
	require.NoError(t, err)
	require.NotNil(t, v)

	assert.Equal(t, "10", v.Chrom)
	assert.Equal(t, int64(116734973), v.Pos)
	assert.Equal(t, "G", v.Ref)
	assert.Equal(t, "A", v.Alt)

	// Read second variant (RRM1)
	v, err = parser.Next()
	require.NoError(t, err)
	require.NotNil(t, v)

	assert.Equal(t, "11", v.Chrom)
	assert.Equal(t, int64(4148284), v.Pos)

	// Read third variant (KRAS G12C)
	v, err = parser.Next()
	require.NoError(t, err)
	require.NotNil(t, v)

	assert.Equal(t, "12", v.Chrom)
	assert.Equal(t, int64(25398285), v.Pos)
	assert.Equal(t, "G", v.Ref)
	assert.Equal(t, "T", v.Alt)

	// Count remaining variants
	count := 3 // Already read 3
	for {
		v, err := parser.Next()
		require.NoError(t, err)
		if v == nil {
			break
		}
		count++
	}

	assert.Equal(t, 5, count)
}

func TestParser_WithAnnotation(t *testing.T) {
	testFile := findTestFile(t, "sample.maf")

	parser, err := NewParser(testFile)
	require.NoError(t, err)
	defer parser.Close()

	// Read first variant with annotation
	v, ann, err := parser.NextWithAnnotation()
	require.NoError(t, err)
	require.NotNil(t, v)
	require.NotNil(t, ann)

	// Verify annotation data
	assert.Equal(t, "TRUB1", ann.HugoSymbol)
	assert.Equal(t, "stop_gained", ann.Consequence)
	assert.Equal(t, "p.W295*", ann.HGVSpShort)
	assert.Equal(t, "ENST00000298746", ann.TranscriptID)
	assert.Equal(t, "GRCh37", ann.NCBIBuild)

	// Read KRAS G12C variant
	parser.NextWithAnnotation() // skip RRM1
	v, ann, err = parser.NextWithAnnotation()
	require.NoError(t, err)

	assert.Equal(t, "KRAS", ann.HugoSymbol)
	assert.Equal(t, "p.G12C", ann.HGVSpShort)
}

func TestParser_Header(t *testing.T) {
	testFile := findTestFile(t, "sample.maf")

	parser, err := NewParser(testFile)
	require.NoError(t, err)
	defer parser.Close()

	header := parser.Header()
	assert.NotEmpty(t, header)

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
	assert.Equal(t, expected, err.Error())
}

func TestParser_ImplementsVariantParser(t *testing.T) {
	testFile := findTestFile(t, "sample.maf")

	parser, err := NewParser(testFile)
	require.NoError(t, err)
	defer parser.Close()

	// Verify the parser has the expected interface methods
	_ = parser.LineNumber()
	v, err := parser.Next()
	require.NoError(t, err)
	require.NotNil(t, v)
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
