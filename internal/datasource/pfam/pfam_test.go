package pfam

import (
	"os"
	"path/filepath"
	"testing"
)

func TestStripVersion(t *testing.T) {
	tests := []struct {
		input, want string
	}{
		{"ENST00000288602.11", "ENST00000288602"},
		{"ENST00000288602", "ENST00000288602"},
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

	// Write a minimal pfamA.txt
	pfamAPath := filepath.Join(dir, "pfamA.txt")
	if err := os.WriteFile(pfamAPath, []byte(
		"pfamA_acc\tpfamA_id\tdescription\n"+
			"PF07714\tPkinase_Tyr\tProtein tyrosine kinase\n"+
			"PF00130\tC1_1\tPhorbol esters/diacylglycerol binding domain\n",
	), 0644); err != nil {
		t.Fatal(err)
	}

	// Write a minimal ensembl_biomart_pfam.txt
	biomartPath := filepath.Join(dir, "ensembl_biomart_pfam.txt")
	if err := os.WriteFile(biomartPath, []byte(
		"Gene stable ID\tTranscript stable ID\tGene name\tPfam domain ID\tPfam domain start\tPfam domain end\n"+
			"ENSG00000157764\tENST00000288602\tBRAF\tPF07714\t457\t712\n"+
			"ENSG00000157764\tENST00000288602\tBRAF\tPF00130\t235\t281\n"+
			"ENSG00000157764\tENST00000288602.5\tBRAF\tPF02196\t157\t225\n",
	), 0644); err != nil {
		t.Fatal(err)
	}

	store, err := Load(pfamAPath, biomartPath)
	if err != nil {
		t.Fatal(err)
	}

	// Test domain lookup.
	d, ok := store.LookupDomain("PF07714")
	if !ok {
		t.Fatal("expected to find PF07714")
	}
	if d.Name != "Pkinase_Tyr" {
		t.Errorf("got name %q, want Pkinase_Tyr", d.Name)
	}

	_, ok = store.LookupDomain("PF99999")
	if ok {
		t.Error("expected not to find PF99999")
	}

	// Test transcript lookup (should merge versioned + unversioned).
	ranges := store.LookupTranscript("ENST00000288602")
	if len(ranges) != 3 {
		t.Errorf("got %d ranges for ENST00000288602, want 3", len(ranges))
	}

	// Test with versioned input.
	ranges2 := store.LookupTranscript("ENST00000288602.11")
	if len(ranges2) != 3 {
		t.Errorf("got %d ranges for versioned input, want 3", len(ranges2))
	}

	// Test counts.
	if store.DomainCount() != 2 {
		t.Errorf("domain count = %d, want 2", store.DomainCount())
	}
	if store.TranscriptCount() != 1 {
		t.Errorf("transcript count = %d, want 1", store.TranscriptCount())
	}
}

func TestLookupDomains(t *testing.T) {
	s := NewStore()
	s.domains["PF00001"] = Domain{Accession: "PF00001", Name: "7tm_1", Description: "7 transmembrane receptor"}
	s.domains["PF00002"] = Domain{Accession: "PF00002", Name: "7tm_2", Description: "7 transmembrane receptor 2"}

	results := s.LookupDomains([]string{"PF00001", "PF99999", "PF00002"})
	if len(results) != 2 {
		t.Fatalf("got %d results, want 2", len(results))
	}
	if results[0].Accession != "PF00001" || results[1].Accession != "PF00002" {
		t.Errorf("unexpected accessions: %v", results)
	}
}
