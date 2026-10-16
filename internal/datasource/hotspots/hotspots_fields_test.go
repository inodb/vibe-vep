package hotspots

import (
	"os"
	"path/filepath"
	"testing"
)

const testHeader = "hugo_symbol\tresidue\tamino_acid_position\tq_value\ttranscript_id\ttranscript_id_version\ttype\tversion\ttumor_count\tmissense_count\ttrunc_count\tinframe_count\tsplice_count\n"

func loadRows(t *testing.T, rows ...string) *Store {
	t.Helper()
	path := filepath.Join(t.TempDir(), "hotspots.txt")
	body := testHeader
	for _, r := range rows {
		body += r + "\n"
	}
	if err := os.WriteFile(path, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}
	store, err := Load(path)
	if err != nil {
		t.Fatalf("load: %v", err)
	}
	return store
}

func deref(t *testing.T, name string, p *int) int {
	t.Helper()
	if p == nil {
		t.Fatalf("%s is nil, want a count", name)
	}
	return *p
}

// A row from the real hotspots_v2_and_3d.txt, trimmed to the columns that
// matter, to pin that every field the API reports is read. The counts are
// written as floats there, which is why they are not read with Atoi.
func TestLoadReadsEveryReportedField(t *testing.T) {
	store := loadRows(t, "NRAS\tQ61\t61\t0\tENST00000369535\tENST00000369535.4\tsingle residue\tv3\t422\t420.0\t1.0\t0.0\t1.0")

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
		{"MissenseCount", deref(t, "MissenseCount", h.MissenseCount), 420},
		{"TruncatingCount", deref(t, "TruncatingCount", h.TruncatingCount), 1},
		{"InframeCount", deref(t, "InframeCount", h.InframeCount), 0},
		{"SpliceCount", deref(t, "SpliceCount", h.SpliceCount), 1},
	} {
		if c.got != c.want {
			t.Errorf("%s = %v, want %v", c.name, c.got, c.want)
		}
	}
}

// In-frame indel hotspots span a range of residues. Their position column
// reads "486-494", which must not be mistaken for an unparseable row.
func TestLoadKeepsRangeHotspots(t *testing.T) {
	store := loadRows(t,
		"BRAF\t486-494\t486-494\t3.95e-09\tENST00000288602\t6\tin-frame indel\tv2\t7\t0.0\t0.0\t7.0\t0.0",
		"BRAF\tV600\t600\t0\tENST00000288602\t6\tsingle residue\tv2\t897\t897.0\t0.0\t0.0\t0.0",
	)

	all := store.LookupAll("ENST00000288602")
	if len(all) != 2 {
		t.Fatalf("LookupAll returned %d hotspots, want 2", len(all))
	}
	rng := all[0]
	if rng.Residue != "486-494" || rng.Position != 486 || rng.PositionEnd != 494 {
		t.Errorf("range hotspot = %q %d-%d, want 486-494", rng.Residue, rng.Position, rng.PositionEnd)
	}
	if got := deref(t, "InframeCount", rng.InframeCount); got != 7 {
		t.Errorf("InframeCount = %d, want 7", got)
	}
	// Annotation asks about one residue, so a span is not an answer to it.
	if _, ok := store.Lookup("ENST00000288602", 486); ok {
		t.Error("Lookup matched a range hotspot; it should consider single residues only")
	}
}

// 3d hotspots leave the per-class columns empty. genome-nexus reports them as
// null, so they must stay distinguishable from a genuine zero.
func TestLoadLeavesAbsentCountsNil(t *testing.T) {
	store := loadRows(t, "BRAF\tG464\t464\t\tENST00000288602\t6\t3d\tv2\t2\t\t\t\t")

	all := store.LookupAll("ENST00000288602")
	if len(all) != 1 {
		t.Fatalf("LookupAll returned %d hotspots, want 1", len(all))
	}
	h := all[0]
	if h.TumorCount != 2 {
		t.Errorf("TumorCount = %d, want 2", h.TumorCount)
	}
	if h.MissenseCount != nil || h.TruncatingCount != nil || h.InframeCount != nil || h.SpliceCount != nil {
		t.Errorf("3d hotspot reported per-class counts: %v %v %v %v",
			h.MissenseCount, h.TruncatingCount, h.InframeCount, h.SpliceCount)
	}
}

// genome-nexus-importer loads this file with
// `--mode upsert --upsertFields hugo_symbol,residue,type,tumor_count`, so rows
// agreeing on those four collapse to one. 3d rows repeat per structural
// cluster and are the usual case; a differing tumour count is a separate row.
func TestLoadDedupsOnTheImportersUpsertKey(t *testing.T) {
	store := loadRows(t,
		"BRAF\tG466\t466\t\tENST00000288602\t6\t3d\tv2\t17\t\t\t\t",
		"BRAF\tG466\t466\t\tENST00000288602\t6\t3d\tv2\t17\t\t\t\t",
		"BRAF\tG466\t466\t\tENST00000288602\t6\t3d\tv2\t17\t\t\t\t",
		"BRAF\tG466\t466\t6.7e-35\tENST00000288602\t6\tsingle residue\tv2\t34\t34.0\t0.0\t0.0\t0.0",
	)

	all := store.LookupAll("ENST00000288602")
	if len(all) != 2 {
		t.Fatalf("LookupAll returned %d hotspots, want 2", len(all))
	}
	seen := map[string]int{}
	for _, h := range all {
		seen[h.Type] = h.TumorCount
	}
	if seen["3d"] != 17 || seen["single residue"] != 34 {
		t.Errorf("got %v, want 3d=17 and single residue=34", seen)
	}
}
