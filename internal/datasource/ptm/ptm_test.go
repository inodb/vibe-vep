package ptm

import (
	"compress/gzip"
	"os"
	"path/filepath"
	"testing"
)

func TestStripVersion(t *testing.T) {
	tests := []struct {
		input, want string
	}{
		{"ENST00000357654.8", "ENST00000357654"},
		{"ENST00000357654", "ENST00000357654"},
		{"", ""},
	}
	for _, tt := range tests {
		got := stripVersion(tt.input)
		if got != tt.want {
			t.Errorf("stripVersion(%q) = %q, want %q", tt.input, got, tt.want)
		}
	}
}

func TestLoadAndLookup(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "ptm.json.gz")

	// Write a minimal ptm.json.gz with JSONL content.
	f, err := os.Create(path)
	if err != nil {
		t.Fatal(err)
	}
	gz := gzip.NewWriter(f)

	lines := []string{
		`{"uniprot_entry":"BRCA1_HUMAN","uniprot_accession":"P38398","position":109,"type":"Sumoylation","pubmed_ids":["28112733"],"sequence":"LEYANSYNFAKKENNSPEHLK","ensembl_transcript_ids":["ENST00000468300.5","ENST00000357654.8"]}`,
		`{"uniprot_entry":"KRAS_HUMAN","uniprot_accession":"P01116","position":12,"type":"Methylation","pubmed_ids":["12345"],"sequence":"MTEYKLVVVGAGGVGKS","ensembl_transcript_ids":["ENST00000311936.8"]}`,
		// Entry with empty transcript IDs — should be skipped.
		`{"uniprot_entry":"SKIP_HUMAN","uniprot_accession":"Q99999","position":1,"type":"Other","pubmed_ids":[],"sequence":"SKIP","ensembl_transcript_ids":[]}`,
	}
	for _, line := range lines {
		gz.Write([]byte(line + "\n"))
	}
	gz.Close()
	f.Close()

	store, err := Load(path)
	if err != nil {
		t.Fatal(err)
	}

	// Lookup by unversioned ID.
	ptms := store.LookupByTranscript("ENST00000357654")
	if len(ptms) != 1 {
		t.Fatalf("got %d PTMs for ENST00000357654, want 1", len(ptms))
	}
	if ptms[0].Type != "Sumoylation" {
		t.Errorf("type = %q, want Sumoylation", ptms[0].Type)
	}
	if ptms[0].UniprotAccession != "P38398" {
		t.Errorf("uniprot_accession = %q, want P38398", ptms[0].UniprotAccession)
	}

	// Lookup by versioned ID.
	ptms2 := store.LookupByTranscript("ENST00000357654.8")
	if len(ptms2) != 1 {
		t.Fatalf("got %d PTMs for versioned input, want 1", len(ptms2))
	}

	// ENST00000468300 should also have the BRCA1 PTM.
	ptms3 := store.LookupByTranscript("ENST00000468300")
	if len(ptms3) != 1 {
		t.Fatalf("got %d PTMs for ENST00000468300, want 1", len(ptms3))
	}

	// KRAS transcript.
	ptms4 := store.LookupByTranscript("ENST00000311936")
	if len(ptms4) != 1 {
		t.Fatalf("got %d PTMs for ENST00000311936, want 1", len(ptms4))
	}
	if ptms4[0].Type != "Methylation" {
		t.Errorf("type = %q, want Methylation", ptms4[0].Type)
	}

	// Unknown transcript.
	ptms5 := store.LookupByTranscript("ENST99999999999")
	if len(ptms5) != 0 {
		t.Errorf("got %d PTMs for unknown transcript, want 0", len(ptms5))
	}

	// Transcript count: ENST00000468300, ENST00000357654, ENST00000311936 = 3 transcripts.
	if store.TranscriptCount() != 3 {
		t.Errorf("transcript count = %d, want 3", store.TranscriptCount())
	}
}

func TestLookupByTranscriptNilStore(t *testing.T) {
	var s *Store
	ptms := s.LookupByTranscript("ENST00000357654")
	if ptms != nil {
		t.Errorf("expected nil from nil store, got %v", ptms)
	}
}
