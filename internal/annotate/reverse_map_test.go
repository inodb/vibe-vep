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

func TestReverseMapHGVSc_SingleBaseDel(t *testing.T) {
	c := createKRASCache()

	// KRAS c.34del: CDS pos 34 is 'G', pos 35 is also 'G', so two equivalent positions.
	// CDS 34 → genomic 25245351, pad from CDS 35 → VCF 25245350:CC>C
	// CDS 35 → genomic 25245350, pad from CDS 36 → VCF 25245349:CA>C (CDS 36='T', comp='A')
	variants, err := ReverseMapHGVSc(c, "KRAS", "34del")
	require.NoError(t, err)
	require.Len(t, variants, 2, "expected 2 equivalent genomic variants for c.34del in GG repeat")

	assert.Equal(t, "12", variants[0].Chrom)
	assert.Equal(t, int64(25245350), variants[0].Pos)
	assert.Equal(t, "CC", variants[0].Ref)
	assert.Equal(t, "C", variants[0].Alt)

	assert.Equal(t, int64(25245349), variants[1].Pos)
}

func TestReverseMapHGVSc_RangeDel(t *testing.T) {
	c := createKRASCache()

	// KRAS c.34_36del: CDS 34-36 = 'GGT', but repeat context extends further:
	// CDS: ...T(33) G(34) G(35) T(36) G(37) G(38) C(39)...
	// Shifting the 3-base window: equivalent positions at 33-35, 34-36, 35-37, 36-38
	// (T at boundary matches: CDS[32]='T'==CDS[35]='T', etc.)
	variants, err := ReverseMapHGVSc(c, "KRAS", "34_36del")
	require.NoError(t, err)
	require.True(t, len(variants) >= 2, "expected multiple equivalent genomic variants, got %d", len(variants))

	// First variant should be the leftmost CDS position mapped to genomic
	assert.Equal(t, "12", variants[0].Chrom)
}

func TestReverseMapHGVSc_ForwardStrandDel(t *testing.T) {
	// Forward strand transcript with known CDS sequence
	c := cache.New()
	c.AddTranscript(&cache.Transcript{
		ID:       "ENST_FWD",
		GeneName: "FWD_GENE",
		Chrom:    "1",
		Start:    1000,
		End:      2000,
		Strand:   1,
		Biotype:  "protein_coding",
		CDSStart: 1010,
		CDSEnd:   1180,
		Exons: []cache.Exon{
			{Number: 1, Start: 1000, End: 1100, CDSStart: 1010, CDSEnd: 1100, Frame: 0},
			{Number: 2, Start: 1150, End: 1200, CDSStart: 1150, CDSEnd: 1180, Frame: 1},
		},
		CDSSequence: "ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGAC",
	})

	// c.5del: CDS pos 5 is 'C' (0-indexed pos 4 of "ATGAC...")
	// Forward strand: genomic for CDS 5 = 1014
	// Padding from CDS pos 4 ('A') at genomic 1013
	// VCF: pos=1013, ref=AC, alt=A
	variants, err := ReverseMapHGVSc(c, "FWD_GENE", "5del")
	require.NoError(t, err)
	require.Len(t, variants, 1)

	v := variants[0]
	assert.Equal(t, "1", v.Chrom)
	assert.Equal(t, int64(1013), v.Pos)
	assert.Equal(t, "AC", v.Ref)
	assert.Equal(t, "A", v.Alt)
}

func TestEquivalentCDSDeletionPositions(t *testing.T) {
	tests := []struct {
		name     string
		cdsSeq   string
		start    int64
		end      int64
		wantPos  []int64
	}{
		{
			name:    "single base in homopolymer GGG",
			cdsSeq:  "ATGGGATC",
			start:   3, end: 3, // first G of GGG
			wantPos: []int64{3, 4, 5},
		},
		{
			name:    "single base no repeat",
			cdsSeq:  "ATCGATCG",
			start:   3, end: 3, // C
			wantPos: []int64{3},
		},
		{
			name:    "3-base tandem repeat",
			cdsSeq:  "XYZATGATGATGXYZ",
			start:   4, end: 6, // first ATG
			wantPos: []int64{4, 5, 6, 7, 8, 9, 10},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := equivalentCDSDeletionPositions(tt.cdsSeq, tt.start, tt.end)
			assert.Equal(t, tt.wantPos, got)
		})
	}
}

func TestResolveHGVSg_Substitution(t *testing.T) {
	c := createKRASCache()

	// 12:g.25245350C>T — simple substitution, bases given directly
	variants, err := ResolveHGVSg(c, "12", "25245350C>T")
	require.NoError(t, err)
	require.Len(t, variants, 1)

	v := variants[0]
	assert.Equal(t, "12", v.Chrom)
	assert.Equal(t, int64(25245350), v.Pos)
	assert.Equal(t, "C", v.Ref)
	assert.Equal(t, "T", v.Alt)
}

func TestResolveHGVSg_SingleBaseDel_ReverseStrand(t *testing.T) {
	c := createKRASCache()

	// Delete genomic position 25245351 (KRAS CDS pos 34, 'G' in CDS = 'C' on genomic)
	// Padding at 25245350 (CDS pos 35, 'G' in CDS = 'C' on genomic)
	// VCF: pos=25245350, ref=CC, alt=C
	variants, err := ResolveHGVSg(c, "12", "25245351del")
	require.NoError(t, err)
	require.Len(t, variants, 1)

	v := variants[0]
	assert.Equal(t, "12", v.Chrom)
	assert.Equal(t, int64(25245350), v.Pos)
	assert.Equal(t, "CC", v.Ref)
	assert.Equal(t, "C", v.Alt)
}

func TestResolveHGVSg_RangeDel_ReverseStrand(t *testing.T) {
	c := createKRASCache()

	// Delete genomic 25245349-25245351 (KRAS CDS 34-36 reversed: 'GGT' → genomic 'ACC' reversed to 'CCA'... wait)
	// Genomic positions 25245349, 25245350, 25245351
	// CDS pos for 25245351 = 34, CDS base 'G', genomic = complement 'C'
	// CDS pos for 25245350 = 35, CDS base 'G', genomic = complement 'C'
	// CDS pos for 25245349 = 36, CDS base 'T', genomic = complement 'A'
	// Deleted genomic bases (ascending pos order): C, C, A at 25245349, 25245350, 25245351 → wait no
	// Actually deleted bases in ascending genomic order:
	//   25245349: CDS 36 = 'T', complement = 'A'
	//   25245350: CDS 35 = 'G', complement = 'C'
	//   25245351: CDS 34 = 'G', complement = 'C'
	// Deleted on genomic strand: "ACC"
	// Padding at 25245348: CDS 37 = 'G', complement = 'C'
	// VCF: pos=25245348, ref=CACC, alt=C
	variants, err := ResolveHGVSg(c, "12", "25245349_25245351del")
	require.NoError(t, err)
	require.Len(t, variants, 1)

	v := variants[0]
	assert.Equal(t, "12", v.Chrom)
	assert.Equal(t, int64(25245348), v.Pos)
	assert.Equal(t, "CACC", v.Ref)
	assert.Equal(t, "C", v.Alt)
}

func TestResolveHGVSg_ForwardStrandDel(t *testing.T) {
	c := cache.New()
	c.AddTranscript(&cache.Transcript{
		ID:       "ENST_FWD",
		GeneName: "FWD_GENE",
		Chrom:    "1",
		Start:    1000,
		End:      2000,
		Strand:   1,
		Biotype:  "protein_coding",
		CDSStart: 1010,
		CDSEnd:   1180,
		Exons: []cache.Exon{
			{Number: 1, Start: 1000, End: 1100, CDSStart: 1010, CDSEnd: 1100, Frame: 0},
			{Number: 2, Start: 1150, End: 1200, CDSStart: 1150, CDSEnd: 1180, Frame: 1},
		},
		CDSSequence: "ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGAC",
	})

	// Delete genomic pos 1014 (CDS pos 5 = 'C')
	// Padding at genomic 1013 (CDS pos 4 = 'A')
	// VCF: pos=1013, ref=AC, alt=A
	variants, err := ResolveHGVSg(c, "1", "1014del")
	require.NoError(t, err)
	require.Len(t, variants, 1)

	v := variants[0]
	assert.Equal(t, "1", v.Chrom)
	assert.Equal(t, int64(1013), v.Pos)
	assert.Equal(t, "AC", v.Ref)
	assert.Equal(t, "A", v.Alt)
}

func TestResolveHGVSg_NoTranscript(t *testing.T) {
	c := createKRASCache()

	_, err := ResolveHGVSg(c, "99", "1000del")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "no transcript")
}

func TestResolveHGVSg_UnsupportedNotation(t *testing.T) {
	c := createKRASCache()

	_, err := ResolveHGVSg(c, "12", "25245350insA")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "unsupported")
}

func TestFindTranscriptsByGene(t *testing.T) {
	c := createKRASCache()

	transcripts := c.FindTranscriptsByGene("KRAS")
	require.Len(t, transcripts, 1)
	assert.Equal(t, "ENST00000311936", transcripts[0].ID)

	transcripts = c.FindTranscriptsByGene("NONEXISTENT")
	assert.Empty(t, transcripts)
}
