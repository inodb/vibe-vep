package server

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/datasource/hotspots"
)

func TestHotspotVersionMatching(t *testing.T) {
	// genome-nexus's rule is deliberately asymmetric: v3 takes everything,
	// v2 takes only rows that say v2 or say nothing at all.
	cases := []struct {
		requested, row string
		want           bool
	}{
		{"v3", "v3", true},
		{"v3", "v2", true},
		{"v3", "", true},
		{"v2", "v2", true},
		{"v2", "", true},
		{"v2", "v3", false},
	}
	for _, c := range cases {
		if got := versionMatches(c.requested, c.row); got != c.want {
			t.Errorf("versionMatches(%q, %q) = %v, want %v", c.requested, c.row, got, c.want)
		}
	}
}

func TestHotspotTypeMatches(t *testing.T) {
	cases := []struct {
		name            string
		hotspotType     string
		classifications []string
		want            bool
	}{
		{"single residue needs a missense", "single residue", []string{"Missense_Mutation"}, true},
		{"single residue rejects a frameshift", "single residue", []string{"Frame_Shift_Del"}, false},
		{"3d also needs a missense", "3d", []string{"Missense_Mutation"}, true},
		{"3d rejects a nonsense", "3d", []string{"Nonsense_Mutation"}, false},
		{"splice site needs a splice", "splice site", []string{"Splice_Site"}, true},
		{"splice site rejects a missense", "splice site", []string{"Missense_Mutation"}, false},
		{"in-frame needs an in-frame indel", "in-frame indel", []string{"In_Frame_Del"}, true},
		{"in-frame rejects a missense", "in-frame indel", []string{"Missense_Mutation"}, false},
		{"one match of several is enough", "single residue",
			[]string{"Splice_Region", "Missense_Mutation"}, true},
		{"an unrecognised type matches anything", "something new", []string{"Silent"}, true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if got := hotspotTypeMatches(c.hotspotType, c.classifications); got != c.want {
				t.Errorf("got %v, want %v", got, c.want)
			}
		})
	}
}

func TestResidueRangeAndOverlap(t *testing.T) {
	for _, c := range []struct {
		residue            string
		wantStart, wantEnd int
		ok                 bool
	}{
		{"Q61", 61, 61, true},
		{"V600_K601", 600, 601, true},
		{"", 0, 0, false},
		{"noDigitsHere", 0, 0, false},
	} {
		start, end, ok := residueRange(c.residue)
		if ok != c.ok || start != c.wantStart || end != c.wantEnd {
			t.Errorf("residueRange(%q) = %d,%d,%v want %d,%d,%v",
				c.residue, start, end, ok, c.wantStart, c.wantEnd, c.ok)
		}
	}

	// A single residue inside a span counts as overlapping, in either order.
	if !overlaps(600, 601, 600, 600) || !overlaps(600, 600, 600, 601) {
		t.Error("a residue inside a span should overlap it")
	}
	if overlaps(600, 601, 602, 603) {
		t.Error("adjacent spans should not overlap")
	}
}

func TestRenderHotspotsFiltersByVersion(t *testing.T) {
	spots := []hotspots.Hotspot{
		{TranscriptID: "ENST1", Residue: "Q61", Type: "single residue", Version: "v2", TumorCount: 10},
		{TranscriptID: "ENST1", Residue: "G12", Type: "single residue", Version: "v3", TumorCount: 20},
	}
	if got := renderHotspots(spots, "v3"); len(got) != 2 {
		t.Errorf("v3 should take every row, got %d", len(got))
	}
	got := renderHotspots(spots, "v2")
	if len(got) != 1 || got[0].Residue != "Q61" {
		t.Errorf("v2 should take only the v2 row, got %v", got)
	}
	// Empty renders as an empty slice, never nil, so it serialises as [].
	if got := renderHotspots(nil, "v3"); got == nil {
		t.Error("expected an empty slice, not nil")
	}
}
