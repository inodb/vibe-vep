package cache

import (
	"os"
	"path/filepath"
	"testing"
)

func TestDuckDBLoader(t *testing.T) {
	// Create temp directory for test database
	tmpDir := t.TempDir()
	dbPath := filepath.Join(tmpDir, "test.duckdb")

	// Create loader and schema
	loader, err := NewDuckDBLoader(dbPath)
	if err != nil {
		t.Fatalf("NewDuckDBLoader: %v", err)
	}
	defer loader.Close()

	if err := loader.CreateSchema(); err != nil {
		t.Fatalf("CreateSchema: %v", err)
	}

	// Insert test transcript (KRAS)
	kras := &Transcript{
		ID:              "ENST00000311936",
		GeneID:          "ENSG00000133703",
		GeneName:        "KRAS",
		Chrom:           "12",
		Start:           25205246,
		End:             25250929,
		Strand:          -1,
		Biotype:         "protein_coding",
		IsCanonical:     true,
		CDSStart:        25209798,
		CDSEnd:          25245384,
		CDSSequence:     "ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGT",
		ProteinSequence: "MTEYKLVVVGAG",
		Exons: []Exon{
			{Number: 1, Start: 25250751, End: 25250929, CDSStart: 0, CDSEnd: 0, Frame: -1},
			{Number: 2, Start: 25245274, End: 25245395, CDSStart: 25245274, CDSEnd: 25245384, Frame: 0},
			{Number: 3, Start: 25227234, End: 25227412, CDSStart: 25227234, CDSEnd: 25227412, Frame: 2},
			{Number: 4, Start: 25225614, End: 25225773, CDSStart: 25225614, CDSEnd: 25225773, Frame: 1},
			{Number: 5, Start: 25209798, End: 25209911, CDSStart: 25209798, CDSEnd: 25209911, Frame: 0},
			{Number: 6, Start: 25205246, End: 25205332, CDSStart: 0, CDSEnd: 0, Frame: -1},
		},
	}

	if err := loader.InsertTranscript(kras); err != nil {
		t.Fatalf("InsertTranscript: %v", err)
	}

	// Test TranscriptCount
	count, err := loader.TranscriptCount()
	if err != nil {
		t.Fatalf("TranscriptCount: %v", err)
	}
	if count != 1 {
		t.Errorf("TranscriptCount = %d, want 1", count)
	}

	// Test GetTranscript
	got, err := loader.GetTranscript("ENST00000311936")
	if err != nil {
		t.Fatalf("GetTranscript: %v", err)
	}
	if got == nil {
		t.Fatal("GetTranscript returned nil")
	}
	if got.GeneName != "KRAS" {
		t.Errorf("GeneName = %q, want KRAS", got.GeneName)
	}
	if len(got.Exons) != 6 {
		t.Errorf("len(Exons) = %d, want 6", len(got.Exons))
	}

	// Test FindTranscripts
	found, err := loader.FindTranscripts("12", 25245351)
	if err != nil {
		t.Fatalf("FindTranscripts: %v", err)
	}
	if len(found) != 1 {
		t.Errorf("FindTranscripts len = %d, want 1", len(found))
	}
	if found[0].ID != "ENST00000311936" {
		t.Errorf("FindTranscripts ID = %q, want ENST00000311936", found[0].ID)
	}

	// Test FindTranscripts with no results
	empty, err := loader.FindTranscripts("12", 1)
	if err != nil {
		t.Fatalf("FindTranscripts: %v", err)
	}
	if len(empty) != 0 {
		t.Errorf("FindTranscripts len = %d, want 0", len(empty))
	}

	// Test GetTranscript not found
	notFound, err := loader.GetTranscript("ENST99999999999")
	if err != nil {
		t.Fatalf("GetTranscript not found: %v", err)
	}
	if notFound != nil {
		t.Error("GetTranscript should return nil for missing ID")
	}

	// Test Chromosomes
	chroms, err := loader.Chromosomes()
	if err != nil {
		t.Fatalf("Chromosomes: %v", err)
	}
	if len(chroms) != 1 || chroms[0] != "12" {
		t.Errorf("Chromosomes = %v, want [12]", chroms)
	}
}

func TestDuckDBLoader_LoadToCache(t *testing.T) {
	tmpDir := t.TempDir()
	dbPath := filepath.Join(tmpDir, "test.duckdb")

	loader, err := NewDuckDBLoader(dbPath)
	if err != nil {
		t.Fatalf("NewDuckDBLoader: %v", err)
	}
	defer loader.Close()

	if err := loader.CreateSchema(); err != nil {
		t.Fatalf("CreateSchema: %v", err)
	}

	// Insert multiple transcripts
	transcripts := []*Transcript{
		{
			ID:       "ENST00000001",
			GeneID:   "ENSG00000001",
			GeneName: "GENE1",
			Chrom:    "1",
			Start:    1000,
			End:      2000,
			Strand:   1,
			Biotype:  "protein_coding",
		},
		{
			ID:       "ENST00000002",
			GeneID:   "ENSG00000002",
			GeneName: "GENE2",
			Chrom:    "1",
			Start:    3000,
			End:      4000,
			Strand:   1,
			Biotype:  "protein_coding",
		},
		{
			ID:       "ENST00000003",
			GeneID:   "ENSG00000003",
			GeneName: "GENE3",
			Chrom:    "2",
			Start:    1000,
			End:      2000,
			Strand:   -1,
			Biotype:  "protein_coding",
		},
	}

	for _, tr := range transcripts {
		if err := loader.InsertTranscript(tr); err != nil {
			t.Fatalf("InsertTranscript: %v", err)
		}
	}

	// Test Load single chromosome
	c := New()
	if err := loader.Load(c, "1"); err != nil {
		t.Fatalf("Load: %v", err)
	}
	if c.TranscriptCount() != 2 {
		t.Errorf("TranscriptCount = %d, want 2", c.TranscriptCount())
	}

	// Test LoadAll
	c2 := New()
	if err := loader.LoadAll(c2); err != nil {
		t.Fatalf("LoadAll: %v", err)
	}
	if c2.TranscriptCount() != 3 {
		t.Errorf("TranscriptCount = %d, want 3", c2.TranscriptCount())
	}
}

func TestIsDuckDB(t *testing.T) {
	tests := []struct {
		path string
		want bool
	}{
		{"transcripts.duckdb", true},
		{"cache.db", true},
		{"s3://bucket/cache.duckdb", true},
		{"/home/user/.vep", false},
		{"cache.json", false},
		{"data.gz", false},
	}

	for _, tc := range tests {
		got := IsDuckDB(tc.path)
		if got != tc.want {
			t.Errorf("IsDuckDB(%q) = %v, want %v", tc.path, got, tc.want)
		}
	}
}

func TestDuckDBLoader_ReopenDB(t *testing.T) {
	tmpDir := t.TempDir()
	dbPath := filepath.Join(tmpDir, "persist.duckdb")

	// Create and populate database
	loader1, err := NewDuckDBLoader(dbPath)
	if err != nil {
		t.Fatalf("NewDuckDBLoader: %v", err)
	}

	if err := loader1.CreateSchema(); err != nil {
		t.Fatalf("CreateSchema: %v", err)
	}

	tr := &Transcript{
		ID:       "ENST00000001",
		GeneID:   "ENSG00000001",
		GeneName: "TEST",
		Chrom:    "1",
		Start:    1000,
		End:      2000,
		Strand:   1,
		Biotype:  "protein_coding",
	}
	if err := loader1.InsertTranscript(tr); err != nil {
		t.Fatalf("InsertTranscript: %v", err)
	}
	loader1.Close()

	// Reopen and verify data persists
	loader2, err := NewDuckDBLoader(dbPath)
	if err != nil {
		t.Fatalf("Reopen NewDuckDBLoader: %v", err)
	}
	defer loader2.Close()

	count, err := loader2.TranscriptCount()
	if err != nil {
		t.Fatalf("TranscriptCount: %v", err)
	}
	if count != 1 {
		t.Errorf("After reopen TranscriptCount = %d, want 1", count)
	}

	got, err := loader2.GetTranscript("ENST00000001")
	if err != nil {
		t.Fatalf("GetTranscript: %v", err)
	}
	if got == nil || got.GeneName != "TEST" {
		t.Error("Data not persisted after reopen")
	}
}

func TestDuckDBLoader_FileExists(t *testing.T) {
	tmpDir := t.TempDir()
	dbPath := filepath.Join(tmpDir, "test.duckdb")

	loader, err := NewDuckDBLoader(dbPath)
	if err != nil {
		t.Fatalf("NewDuckDBLoader: %v", err)
	}
	loader.Close()

	// Check file was created
	if _, err := os.Stat(dbPath); os.IsNotExist(err) {
		t.Error("DuckDB file was not created")
	}
}
