package cache

import (
	"strings"
	"testing"
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
			if got != tt.expected {
				t.Errorf("parseHeader(%q) = %q, want %q", tt.header, got, tt.expected)
			}
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
	if err != nil {
		t.Fatalf("parseFASTA() error = %v", err)
	}

	// Check sequence count
	if loader.SequenceCount() != 2 {
		t.Errorf("SequenceCount() = %d, want 2", loader.SequenceCount())
	}

	// Check KRAS sequence (concatenated without newlines)
	seq := loader.GetSequence("ENST00000311936")
	expectedSeq := "ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAG"
	if seq != expectedSeq {
		t.Errorf("GetSequence(ENST00000311936) = %q, want %q", seq, expectedSeq)
	}

	// Check HasSequence
	if !loader.HasSequence("ENST00000311936") {
		t.Error("HasSequence(ENST00000311936) = false, want true")
	}
	if loader.HasSequence("ENST99999999999") {
		t.Error("HasSequence(ENST99999999999) = true, want false")
	}
}

func TestFASTALoader_LoadFile(t *testing.T) {
	loader := NewFASTALoader("../../testdata/sample_cds.fa")

	if err := loader.Load(); err != nil {
		t.Fatalf("Load() error = %v", err)
	}

	// Check KRAS sequence was loaded
	if !loader.HasSequence("ENST00000311936") {
		t.Error("HasSequence(ENST00000311936) = false after loading sample_cds.fa")
	}

	seq := loader.GetSequence("ENST00000311936")
	if seq == "" {
		t.Error("GetSequence(ENST00000311936) returned empty string")
	}

	// KRAS CDS should start with ATG (start codon)
	if !strings.HasPrefix(seq, "ATG") {
		t.Errorf("KRAS CDS should start with ATG, got %q", seq[:min(3, len(seq))])
	}

	// KRAS CDS should end with TAA (stop codon)
	if !strings.HasSuffix(seq, "TAA") {
		t.Errorf("KRAS CDS should end with TAA, got %q", seq[max(0, len(seq)-3):])
	}
}

func TestFASTALoader_GetSequenceWithVersion(t *testing.T) {
	fastaContent := `>ENST00000311936.8|KRAS
ATGACTGAA
`

	loader := NewFASTALoader("")
	if err := loader.parseFASTA(strings.NewReader(fastaContent)); err != nil {
		t.Fatalf("parseFASTA() error = %v", err)
	}

	// Should find sequence with or without version
	tests := []string{
		"ENST00000311936",
		"ENST00000311936.8",
	}

	for _, id := range tests {
		if !loader.HasSequence(id) {
			t.Errorf("HasSequence(%q) = false, want true", id)
		}
		if seq := loader.GetSequence(id); seq != "ATGACTGAA" {
			t.Errorf("GetSequence(%q) = %q, want ATGACTGAA", id, seq)
		}
	}
}

func TestFASTALoader_CDSExtraction(t *testing.T) {
	// Full mRNA with UTR5 + CDS + UTR3; CDS is positions 11-19 (1-based)
	// UTR5: AAAAAAAAAA (10 bases), CDS: ATGCCCGAA (9 bases), UTR3: TTTTTTTTTT (10 bases)
	fastaContent := `>ENST00000999999.1|ENSG00000999999.1|GENE-201|GENE|29|UTR5:1-10|CDS:11-19|UTR3:20-29|
AAAAAAAAAAATGCCCGAATTTTTTTTTT
`

	loader := NewFASTALoader("")
	if err := loader.parseFASTA(strings.NewReader(fastaContent)); err != nil {
		t.Fatalf("parseFASTA() error = %v", err)
	}

	seq := loader.GetSequence("ENST00000999999")
	expected := "ATGCCCGAA"
	if seq != expected {
		t.Errorf("GetSequence(ENST00000999999) = %q, want %q (CDS only)", seq, expected)
	}
}

func TestFASTALoader_CDSExtractionNoCDSAnnotation(t *testing.T) {
	// Header without CDS annotation — should return the full sequence
	fastaContent := `>ENST00000888888.1|ENSG00000888888.1|GENE2-201|GENE2|30|
ATGCCCGAAATGCCCGAAATGCCCGAATTT
`

	loader := NewFASTALoader("")
	if err := loader.parseFASTA(strings.NewReader(fastaContent)); err != nil {
		t.Fatalf("parseFASTA() error = %v", err)
	}

	seq := loader.GetSequence("ENST00000888888")
	expected := "ATGCCCGAAATGCCCGAAATGCCCGAATTT"
	if seq != expected {
		t.Errorf("GetSequence(ENST00000888888) = %q, want %q (full sequence)", seq, expected)
	}
}

func TestFASTALoader_CDSExtractionNoUTR(t *testing.T) {
	// CDS:1-N (no UTR) — should return the entire sequence unchanged
	fastaContent := `>ENST00000777777.1|ENSG00000777777.1|GENE3-201|GENE3|12|CDS:1-12|
ATGCCCGAATAA
`

	loader := NewFASTALoader("")
	if err := loader.parseFASTA(strings.NewReader(fastaContent)); err != nil {
		t.Fatalf("parseFASTA() error = %v", err)
	}

	seq := loader.GetSequence("ENST00000777777")
	expected := "ATGCCCGAATAA"
	if seq != expected {
		t.Errorf("GetSequence(ENST00000777777) = %q, want %q", seq, expected)
	}
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
			if ok != tt.wantOK || start != tt.wantStart || end != tt.wantEnd {
				t.Errorf("parseCDSRange(%q) = (%d, %d, %v), want (%d, %d, %v)",
					tt.header, start, end, ok, tt.wantStart, tt.wantEnd, tt.wantOK)
			}
		})
	}
}
