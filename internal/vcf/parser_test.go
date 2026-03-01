package vcf

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestParser_SingleVariant(t *testing.T) {
	// Find testdata directory
	testFile := findTestFile(t, "kras_g12c.vcf")

	parser, err := NewParser(testFile)
	require.NoError(t, err)
	defer parser.Close()

	// Read the first (and only) variant
	v, err := parser.Next()
	require.NoError(t, err)
	require.NotNil(t, v)

	// Verify KRAS G12C variant (c.34G>T p.G12C)
	// On reverse strand: coding G->T = genomic C->A
	assert.Equal(t, "12", v.Chrom)
	assert.Equal(t, int64(25245351), v.Pos)
	assert.Equal(t, "C", v.Ref)
	assert.Equal(t, "A", v.Alt)

	// Should be a SNV
	assert.True(t, v.IsSNV(), "KRAS G12C should be classified as SNV")

	// No more variants
	v2, err := parser.Next()
	require.NoError(t, err)
	assert.Nil(t, v2)
}

func TestParser_MultipleVariants(t *testing.T) {
	testFile := findTestFile(t, "multi_variant.vcf")

	parser, err := NewParser(testFile)
	require.NoError(t, err)
	defer parser.Close()

	// Count variants
	count := 0
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

func TestParser_Header(t *testing.T) {
	testFile := findTestFile(t, "kras_g12c.vcf")

	parser, err := NewParser(testFile)
	require.NoError(t, err)
	defer parser.Close()

	header := parser.Header()
	require.NotEmpty(t, header)

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

	assert.True(t, hasFileformat, "Missing ##fileformat header")
	assert.True(t, hasChromLine, "Missing #CHROM header line")
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
			require.Len(t, variants, tt.expected)

			// Each variant should have only one alt allele
			for _, split := range variants {
				assert.NotContains(t, split.Alt, ",")
			}
		})
	}
}

func TestParser_SampleColumns(t *testing.T) {
	input := "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n" +
		"12\t25245351\t.\tC\tA\t100\tPASS\tDP=50\tGT:DP\t0/1:30\t0/0:20\n"

	parser, err := NewParserFromReader(strings.NewReader(input))
	require.NoError(t, err)

	// Check sample names
	assert.Equal(t, []string{"TUMOR", "NORMAL"}, parser.SampleNames())

	v, err := parser.Next()
	require.NoError(t, err)
	require.NotNil(t, v)

	assert.Equal(t, "GT:DP\t0/1:30\t0/0:20", v.SampleColumns)
}

func TestParser_NoSampleColumns(t *testing.T) {
	input := "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"12\t25245351\t.\tC\tA\t100\tPASS\tDP=50\n"

	parser, err := NewParserFromReader(strings.NewReader(input))
	require.NoError(t, err)

	assert.Nil(t, parser.SampleNames())

	v, err := parser.Next()
	require.NoError(t, err)
	require.NotNil(t, v)

	assert.Equal(t, "", v.SampleColumns)
}

func TestSplitMultiAllelic_PreservesSampleColumns(t *testing.T) {
	v := &Variant{
		Chrom:         "12",
		Pos:           100,
		Ref:           "A",
		Alt:           "C,T",
		SampleColumns: "GT:DP\t0/1:30",
	}

	variants := SplitMultiAllelic(v)
	require.Len(t, variants, 2)
	assert.Equal(t, "GT:DP\t0/1:30", variants[0].SampleColumns)
	assert.Equal(t, "GT:DP\t0/1:30", variants[1].SampleColumns)
}

func TestParseError(t *testing.T) {
	err := &ParseError{
		Line:    42,
		Message: "expected 8 columns, found 7",
	}

	expected := "vcf parse error at line 42: expected 8 columns, found 7"
	assert.Equal(t, expected, err.Error())
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
