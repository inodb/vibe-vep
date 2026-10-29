package genomeref

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// streamStrings stands in for the FASTA reader, so the test does not need a
// file on disk to exercise Build.
func streamStrings(seqs map[string]string) func(string, func(string) bool, func(string, []byte)) error {
	order := []string{"1", "2", "GL000207.1"}
	return func(_ string, want func(string) bool, emit func(string, []byte)) error {
		for _, name := range order {
			seq, ok := seqs[name]
			if !ok || !want(name) {
				continue
			}
			emit(name, []byte(seq))
		}
		return nil
	}
}

func buildTestStore(t *testing.T, seqs map[string]string) *Store {
	t.Helper()
	path := filepath.Join(t.TempDir(), FileName)
	want := func(chrom string) bool { return chrom == "1" || chrom == "2" }
	if err := Build("ignored.fa", path, want, streamStrings(seqs)); err != nil {
		t.Fatalf("Build: %v", err)
	}
	s, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	t.Cleanup(func() { s.Close() })
	return s
}

func TestRangeReadsBackWhatWasPacked(t *testing.T) {
	// Deliberately not a multiple of four, so the final partial byte is read.
	const chr1 = "ACGTTTGGCCAAGTACGTNNNNNACGTA"
	const chr2 = "TTTTAAAACCCCGGGGA"
	s := buildTestStore(t, map[string]string{"1": chr1, "2": chr2})

	for _, tc := range []struct {
		chrom      string
		start, end int64
		want       string
	}{
		{"1", 1, int64(len(chr1)), chr1},
		{"1", 1, 1, "A"},
		{"1", 5, 9, "TTGGC"},
		{"1", 19, 23, "NNNNN"},
		{"1", 17, 25, "GTNNNNNAC"},
		{"2", 1, int64(len(chr2)), chr2},
		{"2", 14, 17, "GGGA"},
	} {
		got, err := s.Range(tc.chrom, tc.start, tc.end)
		if err != nil {
			t.Fatalf("Range(%s,%d,%d): %v", tc.chrom, tc.start, tc.end, err)
		}
		if string(got) != tc.want {
			t.Errorf("Range(%s,%d,%d) = %q, want %q", tc.chrom, tc.start, tc.end, got, tc.want)
		}
	}
}

func TestRangeClipsAndRefusesWhatItDoesNotHold(t *testing.T) {
	s := buildTestStore(t, map[string]string{"1": "ACGTACGT", "2": "TTTT", "GL000207.1": "AAAA"})

	if got, _ := s.Range("1", 6, 100); string(got) != "CGT" {
		t.Errorf("window past the end = %q, want %q", got, "CGT")
	}
	if got, _ := s.Range("1", -5, 2); string(got) != "AC" {
		t.Errorf("window before the start = %q, want %q", got, "AC")
	}
	if got, _ := s.Range("1", 20, 30); len(got) != 0 {
		t.Errorf("window entirely past the end = %q, want empty", got)
	}
	// "chr1" and "1" name the same chromosome.
	if got, _ := s.Range("chr1", 1, 4); string(got) != "ACGT" {
		t.Errorf("chr-prefixed name = %q, want %q", got, "ACGT")
	}
	// The scaffold was filtered out at build time, so asking for it says so.
	if _, err := s.Range("GL000207.1", 1, 4); err == nil {
		t.Error("a chromosome the index does not hold should be an error")
	}
}

func TestSoftMaskedBasesComeBackUpperCase(t *testing.T) {
	s := buildTestStore(t, map[string]string{"1": "acgtACGT", "2": "T"})
	got, err := s.Range("1", 1, 8)
	if err != nil {
		t.Fatal(err)
	}
	if string(got) != "ACGTACGT" {
		t.Errorf("Range = %q, want %q", got, "ACGTACGT")
	}
}

func TestOpenRejectsAFileThatIsNotAnIndex(t *testing.T) {
	path := filepath.Join(t.TempDir(), FileName)
	if err := os.WriteFile(path, []byte(strings.Repeat("x", 64)), 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := Open(path); err == nil {
		t.Error("Open should refuse a file without the magic")
	}
}
