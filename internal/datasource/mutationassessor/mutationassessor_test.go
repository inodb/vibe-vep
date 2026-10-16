package mutationassessor

import (
	"bytes"
	"compress/gzip"
	"os"
	"path/filepath"
	"testing"
)

// rows captured from mutationassessor_v4_1.tsv.gz.
const sample = "uniprotId\tsv\thgvspShort\tf_score\tf_impact\tmsa\tmav\n" +
	"P00519\t4\tp.L25A\t3.479676644421987\tlow\tP00519.4_full_b0.6_rg0.3\t4\n" +
	"P00519\t4\tp.L25C\t4.799205309229616\tlow\tP00519.4_full_b0.6_rg0.3\t4\n"

// a second file repeating one substitution, to pin the duplicate policy.
const overlap = "uniprotId\tsv\thgvspShort\tf_score\tf_impact\tmsa\tmav\n" +
	"P00519\t4\tp.L25A\t9.999\thigh\tsomething-else\t9\n" +
	"Q9Y6K9\t2\tp.R123W\t2.5\tmedium\tQ9Y6K9.2_full\t2\n"

func writeGz(t *testing.T, dir, name, content string) string {
	t.Helper()
	path := filepath.Join(dir, name)
	var buf bytes.Buffer
	zw := gzip.NewWriter(&buf)
	if _, err := zw.Write([]byte(content)); err != nil {
		t.Fatal(err)
	}
	if err := zw.Close(); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, buf.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

func TestBuildAndLookup(t *testing.T) {
	dir := t.TempDir()
	a := writeGz(t, dir, "ma_1.tsv.gz", sample)
	b := writeGz(t, dir, "ma_2.tsv.gz", overlap)
	db := filepath.Join(dir, "mutation_assessor.sqlite")

	if err := Build(db, []string{a, b}); err != nil {
		t.Fatalf("build: %v", err)
	}
	if !Ready(db) {
		t.Fatal("index should be ready after a build")
	}

	store, err := Open(db)
	if err != nil {
		t.Fatalf("open: %v", err)
	}
	defer store.Close()

	got, ok := store.Lookup("P00519", "p.L25C")
	if !ok {
		t.Fatal("expected a prediction for P00519 p.L25C")
	}
	if got.SV != 4 || got.MAV != 4 {
		t.Errorf("sv/mav = %d/%d, want 4/4", got.SV, got.MAV)
	}
	if got.FunctionalImpactPrediction != "low" {
		t.Errorf("prediction = %q, want low", got.FunctionalImpactPrediction)
	}
	if got.FunctionalImpactScore < 4.79 || got.FunctionalImpactScore > 4.80 {
		t.Errorf("score = %v, want ~4.799", got.FunctionalImpactScore)
	}

	// A substitution present in two files keeps the first prediction rather
	// than failing the load or silently taking the last.
	first, ok := store.Lookup("P00519", "p.L25A")
	if !ok {
		t.Fatal("expected a prediction for the duplicated substitution")
	}
	if first.FunctionalImpactPrediction != "low" {
		t.Errorf("duplicate resolved to %q, want the first file's low", first.FunctionalImpactPrediction)
	}

	// The second file's own rows still load.
	if _, ok := store.Lookup("Q9Y6K9", "p.R123W"); !ok {
		t.Error("expected the second file's rows to load")
	}

	// Misses are misses, not errors.
	if _, ok := store.Lookup("P00519", "p.Z99Z"); ok {
		t.Error("expected no prediction for an absent substitution")
	}
	if _, ok := store.Lookup("", ""); ok {
		t.Error("expected no prediction for empty keys")
	}
}

func TestReadyOnMissingIndex(t *testing.T) {
	if Ready(filepath.Join(t.TempDir(), "absent.sqlite")) {
		t.Error("a missing index should not be ready")
	}
}

// A nil store answers without being checked for first.
func TestNilStore(t *testing.T) {
	var s *Store
	if _, ok := s.Lookup("P00519", "p.L25A"); ok {
		t.Error("nil store should report no prediction")
	}
	if err := s.Close(); err != nil {
		t.Errorf("closing a nil store: %v", err)
	}
}
