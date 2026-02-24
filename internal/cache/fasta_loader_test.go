package cache

import (
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestFASTALoader_ParseHeader(t *testing.T) {
	loader := &FASTALoader{}

	tests := []struct {
		header   string
		expected string
	}{
		// GENCODE format with pipe delimiters
		{">ENST00000311936.8|ENSG00000133703.14|OTTHUMG|KRAS-201|KRAS|567|", "ENST00000311936"},
		// Simple Ensembl format with space
		{">ENST00000311936.8 cds chromosome:GRCh38", "ENST00000311936"},
		// Just the ID
		{">ENST00000311936", "ENST00000311936"},
		// With version, no delimiter
		{">ENST00000311936.8", "ENST00000311936"},
	}

	for _, tt := range tests {
		t.Run(tt.header, func(t *testing.T) {
			got := loader.parseHeader(tt.header)
			assert.Equal(t, tt.expected, got)
		})
	}
}

func TestFASTALoader_ParseFASTA(t *testing.T) {
	fastaContent := `>ENST00000311936.8|ENSG00000133703.14|KRAS-201|KRAS
ATGACTGAATATAAACTTGTGGTAGTTGGAGCT
GGTGGCGTAGGCAAGAGTGCCTTGACGATACAG
>ENST00000000001.1|ENSG00000000001|TEST
ATGCGATCGATCGATCGATCG
`

	loader := NewFASTALoader("")
	err := loader.parseFASTA(strings.NewReader(fastaContent))
	require.NoError(t, err)

	// Check sequence count
	assert.Equal(t, 2, loader.SequenceCount())

	// Check KRAS sequence (concatenated without newlines)
	seq := loader.GetSequence("ENST00000311936")
	expectedSeq := "ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAG"
	assert.Equal(t, expectedSeq, seq)

	// Check HasSequence
	assert.True(t, loader.HasSequence("ENST00000311936"))
	assert.False(t, loader.HasSequence("ENST99999999999"))
}

func TestFASTALoader_LoadFile(t *testing.T) {
	loader := NewFASTALoader("../../testdata/sample_cds.fa")

	require.NoError(t, loader.Load())

	// Check KRAS sequence was loaded
	assert.True(t, loader.HasSequence("ENST00000311936"))

	seq := loader.GetSequence("ENST00000311936")
	assert.NotEmpty(t, seq)

	// KRAS CDS should start with ATG (start codon)
	assert.True(t, strings.HasPrefix(seq, "ATG"), "KRAS CDS should start with ATG, got %q", seq[:min(3, len(seq))])

	// KRAS CDS should end with TAA (stop codon)
	assert.True(t, strings.HasSuffix(seq, "TAA"), "KRAS CDS should end with TAA, got %q", seq[max(0, len(seq)-3):])
}

func TestFASTALoader_GetSequenceWithVersion(t *testing.T) {
	fastaContent := `>ENST00000311936.8|KRAS
ATGACTGAA
`

	loader := NewFASTALoader("")
	require.NoError(t, loader.parseFASTA(strings.NewReader(fastaContent)))

	// Should find sequence with or without version
	tests := []string{
		"ENST00000311936",
		"ENST00000311936.8",
	}

	for _, id := range tests {
		assert.True(t, loader.HasSequence(id), "HasSequence(%q) should be true", id)
		assert.Equal(t, "ATGACTGAA", loader.GetSequence(id), "GetSequence(%q)", id)
	}
}

func TestFASTALoader_CDSExtraction(t *testing.T) {
	// Full mRNA with UTR5 + CDS + UTR3; CDS is positions 11-19 (1-based)
	// UTR5: AAAAAAAAAA (10 bases), CDS: ATGCCCGAA (9 bases), UTR3: TTTTTTTTTT (10 bases)
	fastaContent := `>ENST00000999999.1|ENSG00000999999.1|GENE-201|GENE|29|UTR5:1-10|CDS:11-19|UTR3:20-29|
AAAAAAAAAAATGCCCGAATTTTTTTTTT
`

	loader := NewFASTALoader("")
	require.NoError(t, loader.parseFASTA(strings.NewReader(fastaContent)))

	seq := loader.GetSequence("ENST00000999999")
	assert.Equal(t, "ATGCCCGAA", seq)
}

func TestFASTALoader_CDSExtractionNoCDSAnnotation(t *testing.T) {
	// Header without CDS annotation — should return the full sequence
	fastaContent := `>ENST00000888888.1|ENSG00000888888.1|GENE2-201|GENE2|30|
ATGCCCGAAATGCCCGAAATGCCCGAATTT
`

	loader := NewFASTALoader("")
	require.NoError(t, loader.parseFASTA(strings.NewReader(fastaContent)))

	seq := loader.GetSequence("ENST00000888888")
	assert.Equal(t, "ATGCCCGAAATGCCCGAAATGCCCGAATTT", seq)
}

func TestFASTALoader_CDSExtractionNoUTR(t *testing.T) {
	// CDS:1-N (no UTR) — should return the entire sequence unchanged
	fastaContent := `>ENST00000777777.1|ENSG00000777777.1|GENE3-201|GENE3|12|CDS:1-12|
ATGCCCGAATAA
`

	loader := NewFASTALoader("")
	require.NoError(t, loader.parseFASTA(strings.NewReader(fastaContent)))

	seq := loader.GetSequence("ENST00000777777")
	assert.Equal(t, "ATGCCCGAATAA", seq)
}

func TestParseCDSRange(t *testing.T) {
	tests := []struct {
		header    string
		wantStart int
		wantEnd   int
		wantOK    bool
	}{
		{">ENST00000456328.2|ENSG001|GENE|459|UTR5:1-200|CDS:201-459|UTR3:460-1657|", 201, 459, true},
		{">ENST00000311936.8|ENSG002|KRAS|567|CDS:1-567|", 1, 567, true},
		{">ENST00000311936.8|ENSG002|KRAS|567|", 0, 0, false},
		{">ENST00000311936.8 simple header", 0, 0, false},
	}

	for _, tt := range tests {
		t.Run(tt.header, func(t *testing.T) {
			start, end, ok := parseCDSRange(tt.header)
			assert.Equal(t, tt.wantOK, ok)
			assert.Equal(t, tt.wantStart, start)
			assert.Equal(t, tt.wantEnd, end)
		})
	}
}
