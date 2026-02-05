package cache

import (
	"strings"
	"testing"
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
				if got := result[key]; got != want {
					t.Errorf("parseAttributes()[%q] = %q, want %q", key, got, want)
				}
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
		if got := stripVersion(tt.input); got != tt.expected {
			t.Errorf("stripVersion(%q) = %q, want %q", tt.input, got, tt.expected)
		}
	}
}

func TestParseStrand(t *testing.T) {
	if got := parseStrand("+"); got != 1 {
		t.Errorf("parseStrand(+) = %d, want 1", got)
	}
	if got := parseStrand("-"); got != -1 {
		t.Errorf("parseStrand(-) = %d, want -1", got)
	}
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
	if err != nil {
		t.Fatalf("parseGTF() error = %v", err)
	}

	if len(transcripts) != 1 {
		t.Fatalf("parseGTF() returned %d transcripts, want 1", len(transcripts))
	}

	tr := transcripts["ENST00000311936"]
	if tr == nil {
		t.Fatal("parseGTF() did not return ENST00000311936")
	}

	// Check transcript fields
	if tr.GeneName != "KRAS" {
		t.Errorf("GeneName = %q, want KRAS", tr.GeneName)
	}
	if tr.Chrom != "12" {
		t.Errorf("Chrom = %q, want 12", tr.Chrom)
	}
	if tr.Strand != -1 {
		t.Errorf("Strand = %d, want -1", tr.Strand)
	}
	if !tr.IsCanonical {
		t.Error("IsCanonical = false, want true")
	}
	if tr.Biotype != "protein_coding" {
		t.Errorf("Biotype = %q, want protein_coding", tr.Biotype)
	}

	// Check exons
	if len(tr.Exons) != 2 {
		t.Fatalf("len(Exons) = %d, want 2", len(tr.Exons))
	}

	// Check CDS boundaries
	if tr.CDSStart == 0 || tr.CDSEnd == 0 {
		t.Errorf("CDS not set: CDSStart=%d, CDSEnd=%d", tr.CDSStart, tr.CDSEnd)
	}
}

func TestGTFLoader_LoadFile(t *testing.T) {
	loader := NewGTFLoader("../../testdata/sample.gtf")
	c := New()

	if err := loader.Load(c); err != nil {
		t.Fatalf("Load() error = %v", err)
	}

	// Should have loaded KRAS transcript
	tr := c.GetTranscript("ENST00000311936")
	if tr == nil {
		t.Fatal("GetTranscript(ENST00000311936) returned nil")
	}

	if tr.GeneName != "KRAS" {
		t.Errorf("GeneName = %q, want KRAS", tr.GeneName)
	}

	// Check exon count - KRAS has 6 exons
	if len(tr.Exons) != 6 {
		t.Errorf("len(Exons) = %d, want 6", len(tr.Exons))
	}

	// Check CDS is protein coding
	if !tr.IsProteinCoding() {
		t.Error("IsProteinCoding() = false, want true")
	}
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
	if err != nil {
		t.Fatalf("parseGTF() error = %v", err)
	}

	if len(transcripts) != 1 {
		t.Fatalf("parseGTF() with filter returned %d transcripts, want 1", len(transcripts))
	}

	if _, ok := transcripts["ENST00000311936"]; !ok {
		t.Error("parseGTF() did not return ENST00000311936 for chr12 filter")
	}
}
