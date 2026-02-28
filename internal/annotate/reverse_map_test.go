package annotate

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// createKRASCache returns a cache with the KRAS transcript for reverse mapping tests.
func createKRASCache() *cache.Cache {
	c := cache.New()
	c.AddTranscript(createKRASTranscript())
	return c
}

func TestCDSToGenomic_ReverseStrand(t *testing.T) {
	tr := createKRASTranscript()

	tests := []struct {
		name    string
		cdsPos  int64
		wantPos int64
	}{
		{"CDS pos 1 (start codon first base)", 1, 25245384},
		{"CDS pos 2", 2, 25245383},
		{"CDS pos 3 (end of start codon)", 3, 25245382},
		{"CDS pos 34 (G12C codon pos 1)", 34, 25245351},
		{"CDS pos 35 (G12V codon pos 2)", 35, 25245350},
		{"CDS pos 36 (codon 12 pos 3)", 36, 25245349},
		{"CDS pos 0 (invalid)", 0, 0},
		{"CDS pos -1 (invalid)", -1, 0},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := CDSToGenomic(tt.cdsPos, tr)
			assert.Equal(t, tt.wantPos, got)
		})
	}
}

func TestCDSToGenomic_RoundTrip(t *testing.T) {
	// Verify CDSToGenomic is the inverse of GenomicToCDS
	tr := createKRASTranscript()

	// Test several known CDS positions
	for _, cdsPos := range []int64{1, 2, 3, 34, 35, 36, 100, 111} {
		genomicPos := CDSToGenomic(cdsPos, tr)
		require.NotZero(t, genomicPos, "CDSToGenomic(%d) returned 0", cdsPos)

		roundTrip := GenomicToCDS(genomicPos, tr)
		assert.Equal(t, cdsPos, roundTrip, "round trip failed: CDS %d → genomic %d → CDS %d", cdsPos, genomicPos, roundTrip)
	}
}

func TestCDSToGenomic_ForwardStrand(t *testing.T) {
	// Simple forward strand transcript: two coding exons
	tr := &cache.Transcript{
		ID:       "ENST_FWD",
		Chrom:    "1",
		Start:    1000,
		End:      2000,
		Strand:   1,
		Biotype:  "protein_coding",
		CDSStart: 1010,
		CDSEnd:   1180,
		Exons: []cache.Exon{
			{Number: 1, Start: 1000, End: 1050, CDSStart: 1010, CDSEnd: 1050, Frame: 0}, // 41 bp CDS
			{Number: 2, Start: 1100, End: 1200, CDSStart: 1100, CDSEnd: 1180, Frame: 2}, // 81 bp CDS
		},
	}

	tests := []struct {
		cdsPos  int64
		wantPos int64
	}{
		{1, 1010},   // First CDS base
		{41, 1050},  // Last base of exon 1 CDS
		{42, 1100},  // First base of exon 2 CDS
		{122, 1180}, // Last CDS base
	}

	for _, tt := range tests {
		got := CDSToGenomic(tt.cdsPos, tr)
		assert.Equal(t, tt.wantPos, got, "CDSToGenomic(%d)", tt.cdsPos)

		// Round trip
		roundTrip := GenomicToCDS(got, tr)
		assert.Equal(t, tt.cdsPos, roundTrip, "round trip CDS %d → genomic %d → CDS %d", tt.cdsPos, got, roundTrip)
	}
}

func TestReverseMapProteinChange_KRASG12C(t *testing.T) {
	c := createKRASCache()

	variants, err := ReverseMapProteinChange(c, "KRAS", 'G', 12, 'C')
	require.NoError(t, err)
	require.Len(t, variants, 1, "expected exactly 1 genomic variant for G12C")

	v := variants[0]
	assert.Equal(t, "12", v.Chrom)
	assert.Equal(t, int64(25245351), v.Pos)
	assert.Equal(t, "C", v.Ref)
	assert.Equal(t, "A", v.Alt)
}

func TestReverseMapProteinChange_KRASG12V(t *testing.T) {
	c := createKRASCache()

	variants, err := ReverseMapProteinChange(c, "KRAS", 'G', 12, 'V')
	require.NoError(t, err)
	require.Len(t, variants, 1, "expected exactly 1 genomic variant for G12V")

	v := variants[0]
	assert.Equal(t, "12", v.Chrom)
	assert.Equal(t, int64(25245350), v.Pos)
	assert.Equal(t, "C", v.Ref)
	assert.Equal(t, "A", v.Alt)
}

func TestReverseMapProteinChange_KRASG12D(t *testing.T) {
	c := createKRASCache()

	variants, err := ReverseMapProteinChange(c, "KRAS", 'G', 12, 'D')
	require.NoError(t, err)
	require.NotEmpty(t, variants)

	// Verify all variants produce the correct annotation
	for _, v := range variants {
		assert.Equal(t, "12", v.Chrom)
	}
}

func TestReverseMapProteinChange_GeneNotFound(t *testing.T) {
	c := createKRASCache()

	_, err := ReverseMapProteinChange(c, "NONEXISTENT", 'G', 12, 'C')
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "not found")
}

func TestReverseMapProteinChange_RefMismatch(t *testing.T) {
	c := createKRASCache()

	// Position 12 is Glycine (G), not Alanine (A)
	_, err := ReverseMapProteinChange(c, "KRAS", 'A', 12, 'C')
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "mismatch")
}

func TestReverseMapHGVSc_KRAS_c35GT(t *testing.T) {
	c := createKRASCache()

	variants, err := ReverseMapHGVSc(c, "KRAS", "35G>T")
	require.NoError(t, err)
	require.Len(t, variants, 1)

	v := variants[0]
	assert.Equal(t, "12", v.Chrom)
	assert.Equal(t, int64(25245350), v.Pos)
	assert.Equal(t, "C", v.Ref) // Complement of G on reverse strand
	assert.Equal(t, "A", v.Alt) // Complement of T on reverse strand
}

func TestReverseMapHGVSc_KRAS_c34GT(t *testing.T) {
	c := createKRASCache()

	variants, err := ReverseMapHGVSc(c, "KRAS", "34G>T")
	require.NoError(t, err)
	require.Len(t, variants, 1)

	v := variants[0]
	assert.Equal(t, "12", v.Chrom)
	assert.Equal(t, int64(25245351), v.Pos)
	assert.Equal(t, "C", v.Ref)
	assert.Equal(t, "A", v.Alt)
}

func TestReverseMapHGVSc_GeneNotFound(t *testing.T) {
	c := createKRASCache()

	_, err := ReverseMapHGVSc(c, "NONEXISTENT", "35G>T")
	assert.Error(t, err)
}

func TestReverseMapHGVSc_UnsupportedNotation(t *testing.T) {
	c := createKRASCache()

	_, err := ReverseMapHGVSc(c, "KRAS", "35delG")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "unsupported")
}

func TestReverseMapHGVSc_ByTranscriptID(t *testing.T) {
	c := createKRASCache()

	variants, err := ReverseMapHGVSc(c, "ENST00000311936", "35G>T")
	require.NoError(t, err)
	require.Len(t, variants, 1)

	v := variants[0]
	assert.Equal(t, int64(25245350), v.Pos)
}

func TestFindTranscriptsByGene(t *testing.T) {
	c := createKRASCache()

	transcripts := c.FindTranscriptsByGene("KRAS")
	require.Len(t, transcripts, 1)
	assert.Equal(t, "ENST00000311936", transcripts[0].ID)

	transcripts = c.FindTranscriptsByGene("NONEXISTENT")
	assert.Empty(t, transcripts)
}
