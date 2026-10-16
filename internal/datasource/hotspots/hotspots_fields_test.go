package hotspots

import (
	"os"
	"path/filepath"
	"testing"
)

// A row from the real hotspots_v2_and_3d.txt, trimmed to the columns that
// matter, to pin that every field the API reports is read.
func TestLoadReadsEveryReportedField(t *testing.T) {
	header := "hugo_symbol\tresidue\tamino_acid_position\tq_value\ttranscript_id\ttranscript_id_version\ttype\tversion\ttumor_count\tmissense_count\ttrunc_count\tinframe_count\tsplice_count\n"
	row := "NRAS\tQ61\t61\t0\tENST00000369535\tENST00000369535.4\tsingle residue\tv3\t422\t420\t1\t0\t1\n"

	dir := t.TempDir()
	path := filepath.Join(dir, "hotspots.txt")
	if err := os.WriteFile(path, []byte(header+row), 0o644); err != nil {
		t.Fatal(err)
	}

	store, err := Load(path)
	if err != nil {
		t.Fatalf("load: %v", err)
	}
	h, ok := store.Lookup("ENST00000369535", 61)
	if !ok {
		t.Fatal("expected NRAS Q61 to be a hotspot")
	}
	for _, c := range []struct {
		name string
		got  any
		want any
	}{
		{"HugoSymbol", h.HugoSymbol, "NRAS"},
		{"Residue", h.Residue, "Q61"},
		{"TranscriptIDVersion", h.TranscriptIDVersion, "ENST00000369535.4"},
		{"Version", h.Version, "v3"},
		{"Type", h.Type, "single residue"},
		{"TumorCount", h.TumorCount, 422},
		{"MissenseCount", h.MissenseCount, 420},
		{"TruncatingCount", h.TruncatingCount, 1},
		{"InframeCount", h.InframeCount, 0},
		{"SpliceCount", h.SpliceCount, 1},
	} {
		if c.got != c.want {
			t.Errorf("%s = %v, want %v", c.name, c.got, c.want)
		}
	}
}
