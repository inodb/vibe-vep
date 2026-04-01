package uniprot

import (
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
	path := filepath.Join(dir, "mapping.txt")

	content := "enst_id\tfinal_uniprot_id\n" +
		"ENST00000357654\tP38398\n" +
		"ENST00000311936\tP01116\n" +
		"ENST00000000233\tP84085\n"

	if err := os.WriteFile(path, []byte(content), 0644); err != nil {
		t.Fatal(err)
	}

	store, err := Load(path)
	if err != nil {
		t.Fatal(err)
	}

	// Lookup existing transcript.
	if got := store.LookupByTranscript("ENST00000357654"); got != "P38398" {
		t.Errorf("got %q, want P38398", got)
	}

	// Lookup with version suffix.
	if got := store.LookupByTranscript("ENST00000357654.8"); got != "P38398" {
		t.Errorf("got %q for versioned input, want P38398", got)
	}

	// Lookup unknown.
	if got := store.LookupByTranscript("ENST99999999999"); got != "" {
		t.Errorf("got %q for unknown, want empty", got)
	}

	// Count.
	if store.Count() != 3 {
		t.Errorf("count = %d, want 3", store.Count())
	}
}

func TestLookupByTranscriptNilStore(t *testing.T) {
	var s *Store
	if got := s.LookupByTranscript("ENST00000357654"); got != "" {
		t.Errorf("expected empty from nil store, got %q", got)
	}
}
