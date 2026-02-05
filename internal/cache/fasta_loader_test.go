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
