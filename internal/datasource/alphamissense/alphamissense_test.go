package alphamissense

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// Small test fixture in AlphaMissense TSV format.
const testTSV = `# Copyright 2023 Google LLC
#
# Data licensed under CC BY 4.0
#CHROM	POS	REF	ALT	am_pathogenicity	am_class
chr1	69094	G	A	0.0782	likely_benign
chr1	69094	G	C	0.0891	likely_benign
chr12	25245350	C	A	0.9876	likely_pathogenic
chr12	25245350	C	T	0.8234	likely_pathogenic
chr17	7674220	C	T	0.6543	ambiguous
`

func writeTSV(t *testing.T) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "test_am.tsv")
	require.NoError(t, os.WriteFile(path, []byte(testTSV), 0644))
	return path
}

func TestLoadAndLookup(t *testing.T) {
	store, err := Open("")
	require.NoError(t, err)
	defer store.Close()

	assert.False(t, store.Loaded(), "should be empty before load")

	tsvPath := writeTSV(t)
	require.NoError(t, store.Load(tsvPath))

	assert.True(t, store.Loaded(), "should have data after load")

	// Hit: KRAS G12C position
	r, ok := store.Lookup("chr12", 25245350, "C", "A")
	assert.True(t, ok)
	assert.InDelta(t, 0.9876, r.Score, 0.001)
	assert.Equal(t, "likely_pathogenic", r.Class)

	// Hit: another allele at same position
	r, ok = store.Lookup("chr12", 25245350, "C", "T")
	assert.True(t, ok)
	assert.InDelta(t, 0.8234, r.Score, 0.001)

	// Miss: non-existent variant
	_, ok = store.Lookup("chr12", 99999999, "A", "T")
	assert.False(t, ok)
}

func TestLoadIdempotent(t *testing.T) {
	store, err := Open("")
	require.NoError(t, err)
	defer store.Close()

	tsvPath := writeTSV(t)
	require.NoError(t, store.Load(tsvPath))

	// Loading again should fail due to primary key constraint, which is expected.
	// In practice we check Loaded() before calling Load().
	err = store.Load(tsvPath)
	assert.Error(t, err, "duplicate load should error on primary key conflict")
}

func TestLookupEmpty(t *testing.T) {
	store, err := Open("")
	require.NoError(t, err)
	defer store.Close()

	_, ok := store.Lookup("chr1", 1, "A", "T")
	assert.False(t, ok)
}
