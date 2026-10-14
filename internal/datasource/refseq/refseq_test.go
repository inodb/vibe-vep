package refseq

import (
	"compress/gzip"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// GENCODE metadata.RefSeq has NO header line — the first line is data.
const sampleTSV = "" +
	"ENST00000269305.4\tNM_000546.5\tNP_000537.3\n" +
	"ENST00000269305.4\tNM_001126112.2\tNP_001119584.1\n" +
	"ENST00000269305.4\tNM_001126118.1\tNP_001119590.1\n" +
	"ENST00000311936.3\tNM_004985.3\tNP_004976.2\n" +
	"ENST00000569605.1\tNR_015445.1\t\n" + // non-coding: no protein column value
	"ENST00000441802.2\tNM_005245.3\tNP_005236.2\n"

func writeTemp(t *testing.T, name, content string, gz bool) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), name)
	f, err := os.Create(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	if gz {
		w := gzip.NewWriter(f)
		if _, err := w.Write([]byte(content)); err != nil {
			t.Fatal(err)
		}
		if err := w.Close(); err != nil {
			t.Fatal(err)
		}
		return path
	}
	if _, err := f.WriteString(content); err != nil {
		t.Fatal(err)
	}
	return path
}

func TestLoadPlain(t *testing.T) {
	s, err := Load(writeTemp(t, "refseq.tsv", sampleTSV, false))
	if err != nil {
		t.Fatalf("Load: %v", err)
	}
	if got := s.Count(); got != 4 {
		t.Errorf("Count()=%d, want 4", got)
	}
}

func TestLoadGzip(t *testing.T) {
	s, err := Load(writeTemp(t, "refseq.tsv.gz", sampleTSV, true))
	if err != nil {
		t.Fatalf("Load: %v", err)
	}
	if got := s.LookupByTranscript("ENST00000311936.3"); len(got) != 1 || got[0] != "NM_004985.3" {
		t.Errorf("LookupByTranscript = %v, want [NM_004985.3]", got)
	}
}

// The first line of the file is data, not a header. A copy-pasted "skip header"
// would silently drop TP53's primary accession, which is exactly the value
// genome-nexus reports.
func TestFirstLineIsData(t *testing.T) {
	s, err := Load(writeTemp(t, "refseq.tsv", sampleTSV, false))
	if err != nil {
		t.Fatal(err)
	}
	got := s.LookupByTranscript("ENST00000269305.4")
	if len(got) == 0 || got[0] != "NM_000546.5" {
		t.Fatalf("first accession = %v, want NM_000546.5 first", got)
	}
}

// Order is load-bearing: genome-nexus's RefSeqResolver reports element 0 only.
func TestMultiAccessionOrderPreserved(t *testing.T) {
	s, err := Load(writeTemp(t, "refseq.tsv", sampleTSV, false))
	if err != nil {
		t.Fatal(err)
	}
	want := []string{"NM_000546.5", "NM_001126112.2", "NM_001126118.1"}
	got := s.LookupByTranscript("ENST00000269305.4")
	if len(got) != len(want) {
		t.Fatalf("got %d accessions %v, want %d", len(got), got, len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("accession[%d]=%q, want %q", i, got[i], want[i])
		}
	}
}

func TestStrippedVersionFallback(t *testing.T) {
	s, err := Load(writeTemp(t, "refseq.tsv", sampleTSV, false))
	if err != nil {
		t.Fatal(err)
	}
	// A different transcript version still resolves, so a GENCODE release skew
	// between the GTF and this file degrades gracefully instead of to zero.
	if got := s.LookupByTranscript("ENST00000311936.8"); len(got) != 1 || got[0] != "NM_004985.3" {
		t.Errorf("stripped-version lookup = %v, want [NM_004985.3]", got)
	}
	// Map() is what the loader joins on, and is keyed unversioned.
	if got := s.Map()["ENST00000311936"]; len(got) != 1 {
		t.Errorf("Map()[unversioned] = %v, want 1 accession", got)
	}
}

func TestUnknownTranscriptReturnsNil(t *testing.T) {
	s, err := Load(writeTemp(t, "refseq.tsv", sampleTSV, false))
	if err != nil {
		t.Fatal(err)
	}
	if got := s.LookupByTranscript("ENST99999999.1"); len(got) != 0 {
		t.Errorf("unknown transcript = %v, want empty", got)
	}
}

func TestMalformedLinesSkipped(t *testing.T) {
	input := "\n" +
		"ENST00000000001.1\n" + // single field
		"ENST00000000002.1\t\n" + // empty accession
		"\t\t\n" + // no transcript
		"# comment line\n" +
		"ENST00000311936.3\tNM_004985.3\tNP_004976.2\n"
	s, err := Load(writeTemp(t, "refseq.tsv", input, false))
	if err != nil {
		t.Fatalf("Load: %v", err)
	}
	if got := s.Count(); got != 1 {
		t.Errorf("Count()=%d, want 1 (malformed rows skipped)", got)
	}
}

func TestDuplicateAccessionDeduped(t *testing.T) {
	input := "ENST00000311936.3\tNM_004985.3\t\nENST00000311936.3\tNM_004985.3\t\n"
	s, err := Load(writeTemp(t, "refseq.tsv", input, false))
	if err != nil {
		t.Fatal(err)
	}
	if got := s.LookupByTranscript("ENST00000311936.3"); len(got) != 1 {
		t.Errorf("got %v, want a single deduped accession", got)
	}
}

func TestLoadMissingFile(t *testing.T) {
	if _, err := Load(filepath.Join(t.TempDir(), "nope.tsv")); err == nil {
		t.Error("expected an error for a missing file")
	}
}

func TestNilStoreIsSafe(t *testing.T) {
	var s *Store
	if got := s.LookupByTranscript("ENST00000311936.3"); got != nil {
		t.Errorf("nil store lookup = %v, want nil", got)
	}
	if s.Count() != 0 || s.Map() != nil {
		t.Error("nil store should report zero/nil")
	}
}

func TestPARSuffixStripsToBaseID(t *testing.T) {
	// GRCh38 PAR_Y transcripts carry a suffix after the version.
	input := "ENST00000362079.3_PAR_Y\tNM_001185183.2\t\n"
	s, err := Load(writeTemp(t, "refseq.tsv", input, false))
	if err != nil {
		t.Fatal(err)
	}
	for _, key := range []string{"ENST00000362079.3_PAR_Y", "ENST00000362079"} {
		if got := s.LookupByTranscript(key); len(got) != 1 {
			t.Errorf("lookup(%q) = %v, want 1 accession", key, got)
		}
	}
	if !strings.HasPrefix(s.Map()["ENST00000362079"][0], "NM_") {
		t.Error("expected an NM_ accession under the stripped key")
	}
}
