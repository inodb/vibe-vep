package cache

import (
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestParseAttributes(t *testing.T) {
	tests := []struct {
		name     string
		input    string
		expected map[string]string
	}{
		{
			name:  "basic attributes",
			input: `gene_id "ENSG00000133703"; transcript_id "ENST00000311936"; gene_name "KRAS";`,
			expected: map[string]string{
				"gene_id":       "ENSG00000133703",
				"transcript_id": "ENST00000311936",
				"gene_name":     "KRAS",
			},
		},
		{
			name:  "with tags",
			input: `gene_id "ENSG00000133703"; tag "Ensembl_canonical"; tag "MANE_Select";`,
			expected: map[string]string{
				"gene_id": "ENSG00000133703",
				"tag":     "MANE_Select", // Last value wins
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := parseAttributes(tt.input)
			for key, want := range tt.expected {
				assert.Equal(t, want, result[key], "parseAttributes()[%q]", key)
			}
		})
	}
}

func TestStripVersion(t *testing.T) {
	tests := []struct {
		input    string
		expected string
	}{
		{"ENST00000311936.8", "ENST00000311936"},
		{"ENSG00000133703.14", "ENSG00000133703"},
		{"ENST00000311936", "ENST00000311936"},
		{"", ""},
	}

	for _, tt := range tests {
		assert.Equal(t, tt.expected, stripVersion(tt.input), "stripVersion(%q)", tt.input)
	}
}

func TestParseStrand(t *testing.T) {
	assert.Equal(t, int8(1), parseStrand("+"))
	assert.Equal(t, int8(-1), parseStrand("-"))
}

func TestGTFLoader_ParseGTF(t *testing.T) {
	gtfContent := `##description: Test GTF
chr12	HAVANA	gene	25205246	25250929	.	-	.	gene_id "ENSG00000133703"; gene_type "protein_coding"; gene_name "KRAS";
chr12	HAVANA	transcript	25205246	25250929	.	-	.	gene_id "ENSG00000133703"; transcript_id "ENST00000311936"; gene_type "protein_coding"; gene_name "KRAS"; transcript_type "protein_coding"; tag "Ensembl_canonical";
chr12	HAVANA	exon	25250751	25250929	.	-	.	gene_id "ENSG00000133703"; transcript_id "ENST00000311936"; gene_name "KRAS"; exon_number "1";
chr12	HAVANA	exon	25245274	25245395	.	-	.	gene_id "ENSG00000133703"; transcript_id "ENST00000311936"; gene_name "KRAS"; exon_number "2";
chr12	HAVANA	CDS	25250751	25250808	.	-	0	gene_id "ENSG00000133703"; transcript_id "ENST00000311936"; gene_name "KRAS"; exon_number "1";
chr12	HAVANA	CDS	25245274	25245395	.	-	2	gene_id "ENSG00000133703"; transcript_id "ENST00000311936"; gene_name "KRAS"; exon_number "2";
chr12	HAVANA	start_codon	25250806	25250808	.	-	0	gene_id "ENSG00000133703"; transcript_id "ENST00000311936"; gene_name "KRAS";
chr12	HAVANA	stop_codon	25245274	25245276	.	-	0	gene_id "ENSG00000133703"; transcript_id "ENST00000311936"; gene_name "KRAS";
`

	loader := &GTFLoader{}
	transcripts, err := loader.parseGTF(strings.NewReader(gtfContent), "")
	require.NoError(t, err)

	require.Len(t, transcripts, 1)

	tr := transcripts["ENST00000311936"]
	require.NotNil(t, tr, "parseGTF() did not return ENST00000311936")

	// Check transcript fields
	assert.Equal(t, "KRAS", tr.GeneName)
	assert.Equal(t, "12", tr.Chrom)
	assert.Equal(t, int8(-1), tr.Strand)
	assert.True(t, tr.IsCanonical)
	assert.Equal(t, "protein_coding", tr.Biotype)

	// Check exons
	require.Len(t, tr.Exons, 2)

	// Check CDS boundaries
	assert.NotZero(t, tr.CDSStart, "CDSStart should be set")
	assert.NotZero(t, tr.CDSEnd, "CDSEnd should be set")
}

func TestGTFLoader_LoadFile(t *testing.T) {
	loader := NewGTFLoader("../../testdata/sample.gtf")
	c := New()

	require.NoError(t, loader.Load(c))

	// Should have loaded KRAS transcript
	tr := c.GetTranscript("ENST00000311936")
	require.NotNil(t, tr, "GetTranscript(ENST00000311936) returned nil")

	assert.Equal(t, "KRAS", tr.GeneName)

	// Check exon count - KRAS has 6 exons
	assert.Len(t, tr.Exons, 6)

	// Check CDS is protein coding
	assert.True(t, tr.IsProteinCoding())
}

func TestGTFLoader_FilterChromosome(t *testing.T) {
	gtfContent := `chr12	HAVANA	transcript	25205246	25250929	.	-	.	gene_id "ENSG00000133703"; transcript_id "ENST00000311936"; gene_name "KRAS"; transcript_type "protein_coding";
chr12	HAVANA	exon	25250751	25250929	.	-	.	gene_id "ENSG00000133703"; transcript_id "ENST00000311936"; exon_number "1";
chr1	HAVANA	transcript	100000	200000	.	+	.	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; gene_name "TEST"; transcript_type "protein_coding";
chr1	HAVANA	exon	100000	100100	.	+	.	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "1";
`

	loader := &GTFLoader{}

	// Filter to chr12 only
	transcripts, err := loader.parseGTF(strings.NewReader(gtfContent), "chr12")
	require.NoError(t, err)

	require.Len(t, transcripts, 1)

	assert.Contains(t, transcripts, "ENST00000311936")
}
