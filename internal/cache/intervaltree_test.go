package cache

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestBuildIntervalTree_Empty(t *testing.T) {
	tree := BuildIntervalTree(nil)
	assert.Empty(t, tree.FindOverlaps(100))
}

func TestIntervalTree_SingleTranscript(t *testing.T) {
	tx := &Transcript{ID: "ENST001", Start: 100, End: 200}
	tree := BuildIntervalTree([]*Transcript{tx})

	assert.Len(t, tree.FindOverlaps(150), 1)
	assert.Equal(t, "ENST001", tree.FindOverlaps(150)[0].ID)

	assert.Len(t, tree.FindOverlaps(100), 1, "start boundary inclusive")
	assert.Len(t, tree.FindOverlaps(200), 1, "end boundary inclusive")
	assert.Empty(t, tree.FindOverlaps(99), "before start")
	assert.Empty(t, tree.FindOverlaps(201), "after end")
}

func TestIntervalTree_Overlapping(t *testing.T) {
	transcripts := []*Transcript{
		{ID: "A", Start: 100, End: 300},
		{ID: "B", Start: 150, End: 250},
		{ID: "C", Start: 200, End: 400},
	}
	tree := BuildIntervalTree(transcripts)

	results := tree.FindOverlaps(175)
	assert.Len(t, results, 2, "pos 175 overlaps A and B")
	ids := map[string]bool{}
	for _, r := range results {
		ids[r.ID] = true
	}
	assert.True(t, ids["A"])
	assert.True(t, ids["B"])

	results = tree.FindOverlaps(250)
	assert.Len(t, results, 3, "pos 250 overlaps A, B, C")

	results = tree.FindOverlaps(350)
	assert.Len(t, results, 1, "pos 350 overlaps only C")
	assert.Equal(t, "C", results[0].ID)
}

func TestIntervalTree_NonOverlapping(t *testing.T) {
	transcripts := []*Transcript{
		{ID: "A", Start: 100, End: 200},
		{ID: "B", Start: 300, End: 400},
		{ID: "C", Start: 500, End: 600},
	}
	tree := BuildIntervalTree(transcripts)

	assert.Len(t, tree.FindOverlaps(150), 1)
	assert.Equal(t, "A", tree.FindOverlaps(150)[0].ID)

	assert.Empty(t, tree.FindOverlaps(250), "gap between A and B")

	assert.Len(t, tree.FindOverlaps(350), 1)
	assert.Equal(t, "B", tree.FindOverlaps(350)[0].ID)
}

func TestIntervalTree_MaxEndPruning(t *testing.T) {
	// A short interval followed by a long one — maxEnd should allow finding the long one
	transcripts := []*Transcript{
		{ID: "short", Start: 100, End: 110},
		{ID: "long", Start: 105, End: 500},
	}
	tree := BuildIntervalTree(transcripts)

	results := tree.FindOverlaps(400)
	assert.Len(t, results, 1)
	assert.Equal(t, "long", results[0].ID)
}

func TestIntervalTree_FindOverlapsWithinDistance(t *testing.T) {
	transcripts := []*Transcript{
		{ID: "A", Start: 1000, End: 2000},
		{ID: "B", Start: 5000, End: 6000},
	}
	tree := BuildIntervalTree(transcripts)

	// distance=0 behaves like FindOverlaps.
	assert.Empty(t, tree.FindOverlapsWithinDistance(999, 0))
	assert.Len(t, tree.FindOverlapsWithinDistance(1500, 0), 1)

	// 5 bp upstream of A: within a 10 bp padding.
	results := tree.FindOverlapsWithinDistance(995, 10)
	assert.Len(t, results, 1)
	assert.Equal(t, "A", results[0].ID)

	// 5 bp downstream of A: within a 10 bp padding.
	results = tree.FindOverlapsWithinDistance(2005, 10)
	assert.Len(t, results, 1)
	assert.Equal(t, "A", results[0].ID)

	// Just outside the padding on both sides — no hits.
	assert.Empty(t, tree.FindOverlapsWithinDistance(989, 10))
	assert.Empty(t, tree.FindOverlapsWithinDistance(2011, 10))

	// Large padding overlapping both transcripts from between them.
	results = tree.FindOverlapsWithinDistance(3500, 2000)
	ids := map[string]bool{}
	for _, r := range results {
		ids[r.ID] = true
	}
	assert.True(t, ids["A"], "A within 2 kb downstream padding")
	assert.True(t, ids["B"], "B within 2 kb upstream padding")

	// Negative distance is clamped to 0.
	assert.Empty(t, tree.FindOverlapsWithinDistance(999, -5))
	assert.Len(t, tree.FindOverlapsWithinDistance(1500, -100), 1)
}

func TestIntervalTree_MatchesLinearScan(t *testing.T) {
	// Verify interval tree produces same results as linear scan
	transcripts := []*Transcript{
		{ID: "A", Start: 1000, End: 5000},
		{ID: "B", Start: 2000, End: 3000},
		{ID: "C", Start: 4000, End: 8000},
		{ID: "D", Start: 6000, End: 7000},
		{ID: "E", Start: 9000, End: 10000},
	}
	tree := BuildIntervalTree(transcripts)

	for pos := int64(0); pos <= 11000; pos += 500 {
		// Linear scan
		var linear []*Transcript
		for _, tx := range transcripts {
			if tx.Contains(pos) {
				linear = append(linear, tx)
			}
		}
		// Tree query
		treeResult := tree.FindOverlaps(pos)

		linearIDs := map[string]bool{}
		for _, tx := range linear {
			linearIDs[tx.ID] = true
		}
		treeIDs := map[string]bool{}
		for _, tx := range treeResult {
			treeIDs[tx.ID] = true
		}

		assert.Equal(t, linearIDs, treeIDs, "pos=%d", pos)
	}
}
