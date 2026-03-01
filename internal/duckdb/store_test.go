package duckdb

import (
	"testing"
	"time"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func openInMemory(t *testing.T) *Store {
	t.Helper()
	s, err := Open("")
	require.NoError(t, err)
	t.Cleanup(func() { s.Close() })
	return s
}

// --- Variant cache tests (DuckDB) ---

func TestOpenClose(t *testing.T) {
	s := openInMemory(t)
	assert.NotNil(t, s.DB())
}

func TestWriteAndLookupVariants(t *testing.T) {
	s := openInMemory(t)

	results := []VariantResult{
		{
			Chrom: "12", Pos: 25245350, Ref: "C", Alt: "A",
			Ann: &annotate.Annotation{
				TranscriptID: "ENST00000311936.8",
				GeneName:     "KRAS", GeneID: "ENSG00000133703",
				Consequence: "missense_variant", Impact: "MODERATE",
				CDSPosition: 35, ProteinPosition: 12,
				AminoAcidChange: "G/V", CodonChange: "gGt/gTt",
				IsCanonical: true, Allele: "A",
				Biotype: "protein_coding", ExonNumber: "2/6",
				HGVSp: "p.Gly12Val", HGVSc: "c.35G>T",
			},
		},
		{
			Chrom: "12", Pos: 25245350, Ref: "C", Alt: "A",
			Ann: &annotate.Annotation{
				TranscriptID: "ENST00000256078.10",
				GeneName:     "KRAS", GeneID: "ENSG00000133703",
				Consequence: "missense_variant", Impact: "MODERATE",
				IsCanonical: false, Allele: "A",
				Biotype: "protein_coding",
			},
		},
	}

	err := s.WriteVariantResults(results)
	require.NoError(t, err)

	anns, err := s.LookupVariant("12", 25245350, "C", "A")
	require.NoError(t, err)
	require.Len(t, anns, 2)
	assert.Equal(t, "ENST00000311936.8", anns[0].TranscriptID)
	assert.Equal(t, "p.Gly12Val", anns[0].HGVSp)
	assert.Equal(t, "ENST00000256078.10", anns[1].TranscriptID)

	anns, err = s.LookupVariant("12", 99999, "C", "A")
	require.NoError(t, err)
	assert.Empty(t, anns)
}

func TestLookupVariantEmpty(t *testing.T) {
	s := openInMemory(t)

	anns, err := s.LookupVariant("1", 100, "A", "T")
	require.NoError(t, err)
	assert.Empty(t, anns)
}

func TestClearVariantResults(t *testing.T) {
	s := openInMemory(t)

	results := []VariantResult{
		{
			Chrom: "1", Pos: 100, Ref: "A", Alt: "T",
			Ann: &annotate.Annotation{
				TranscriptID: "ENST00000001.1",
				GeneName:     "TEST",
				Consequence:  "missense_variant", Impact: "MODERATE",
				Allele: "T", Biotype: "protein_coding",
			},
		},
	}
	require.NoError(t, s.WriteVariantResults(results))

	anns, err := s.LookupVariant("1", 100, "A", "T")
	require.NoError(t, err)
	require.Len(t, anns, 1)

	require.NoError(t, s.ClearVariantResults())

	anns, err = s.LookupVariant("1", 100, "A", "T")
	require.NoError(t, err)
	assert.Empty(t, anns)
}

func TestSearchByGene(t *testing.T) {
	s := openInMemory(t)

	results := []VariantResult{
		{
			Chrom: "12", Pos: 25245350, Ref: "C", Alt: "A",
			Ann: &annotate.Annotation{
				TranscriptID: "ENST00000311936.8",
				GeneName:     "KRAS", Consequence: "missense_variant", Impact: "MODERATE",
				AminoAcidChange: "G/V", Allele: "A", Biotype: "protein_coding",
			},
		},
		{
			Chrom: "7", Pos: 140753336, Ref: "A", Alt: "T",
			Ann: &annotate.Annotation{
				TranscriptID: "ENST00000288602.11",
				GeneName:     "BRAF", Consequence: "missense_variant", Impact: "MODERATE",
				Allele: "T", Biotype: "protein_coding",
			},
		},
	}
	require.NoError(t, s.WriteVariantResults(results))

	krasResults, err := s.SearchByGene("KRAS")
	require.NoError(t, err)
	require.Len(t, krasResults, 1)
	assert.Equal(t, "KRAS", krasResults[0].Ann.GeneName)

	brafResults, err := s.SearchByGene("BRAF")
	require.NoError(t, err)
	require.Len(t, brafResults, 1)

	noneResults, err := s.SearchByGene("NOTEXIST")
	require.NoError(t, err)
	assert.Empty(t, noneResults)
}

func TestSearchByProteinChange(t *testing.T) {
	s := openInMemory(t)

	results := []VariantResult{
		{
			Chrom: "12", Pos: 25245350, Ref: "C", Alt: "A",
			Ann: &annotate.Annotation{
				TranscriptID: "ENST00000311936.8",
				GeneName:     "KRAS", Consequence: "missense_variant", Impact: "MODERATE",
				AminoAcidChange: "G/V", Allele: "A", Biotype: "protein_coding",
			},
		},
		{
			Chrom: "12", Pos: 25245351, Ref: "G", Alt: "T",
			Ann: &annotate.Annotation{
				TranscriptID: "ENST00000311936.8",
				GeneName:     "KRAS", Consequence: "missense_variant", Impact: "MODERATE",
				AminoAcidChange: "G/C", Allele: "T", Biotype: "protein_coding",
			},
		},
	}
	require.NoError(t, s.WriteVariantResults(results))

	found, err := s.SearchByProteinChange("KRAS", "G/V")
	require.NoError(t, err)
	require.Len(t, found, 1)
	assert.Equal(t, int64(25245350), found[0].Pos)

	found, err = s.SearchByProteinChange("KRAS", "G/C")
	require.NoError(t, err)
	require.Len(t, found, 1)
	assert.Equal(t, int64(25245351), found[0].Pos)

	found, err = s.SearchByProteinChange("KRAS", "X/Y")
	require.NoError(t, err)
	assert.Empty(t, found)
}

// --- Transcript cache tests (gob) ---

func TestTranscriptCacheWriteAndLoad(t *testing.T) {
	dir := t.TempDir()
	tc := NewTranscriptCache(dir)

	c := cache.New()
	c.AddTranscript(&cache.Transcript{
		ID: "ENST00000311936.8", GeneID: "ENSG00000133703", GeneName: "KRAS",
		Chrom: "12", Start: 25205246, End: 25250936, Strand: -1,
		Biotype: "protein_coding", IsCanonical: true, IsMANESelect: true,
		CDSStart: 25209431, CDSEnd: 25245384,
		CDSSequence: "ATGACTGAATATAAACTTGTGG", ProteinSequence: "MTEYK",
		Exons: []cache.Exon{
			{Number: 1, Start: 25205246, End: 25209911, CDSStart: 25209431, CDSEnd: 25209911, Frame: 0},
			{Number: 2, Start: 25215441, End: 25215560, CDSStart: 25215441, CDSEnd: 25215560, Frame: 2},
		},
	})
	c.AddTranscript(&cache.Transcript{
		ID: "ENST00000256078.10", GeneID: "ENSG00000133703", GeneName: "KRAS",
		Chrom: "12", Start: 25205246, End: 25250936, Strand: -1,
		Biotype: "protein_coding",
		Exons: []cache.Exon{
			{Number: 1, Start: 25205246, End: 25209911, CDSStart: 25209431, CDSEnd: 25209911, Frame: 0},
		},
	})

	now := time.Now()
	fp := FileFingerprint{Size: 1000, ModTime: now}
	require.NoError(t, tc.Write(c, fp, fp, fp))

	// Load back
	c2 := cache.New()
	require.NoError(t, tc.Load(c2))

	assert.Equal(t, 2, c2.TranscriptCount())
	assert.Equal(t, []string{"12"}, c2.Chromosomes())

	kras := c2.GetTranscript("ENST00000311936.8")
	require.NotNil(t, kras)
	assert.Equal(t, "KRAS", kras.GeneName)
	assert.Equal(t, int8(-1), kras.Strand)
	assert.True(t, kras.IsCanonical)
	assert.True(t, kras.IsMANESelect)
	assert.Equal(t, "ATGACTGAATATAAACTTGTGG", kras.CDSSequence)
	assert.Equal(t, "MTEYK", kras.ProteinSequence)

	require.Len(t, kras.Exons, 2)
	assert.Equal(t, 1, kras.Exons[0].Number)
	assert.Equal(t, int64(25205246), kras.Exons[0].Start)
	assert.Equal(t, 0, kras.Exons[0].Frame)
}

func TestTranscriptCacheValidation(t *testing.T) {
	dir := t.TempDir()
	tc := NewTranscriptCache(dir)

	now := time.Now()
	gtf := FileFingerprint{Size: 1000, ModTime: now}
	fasta := FileFingerprint{Size: 2000, ModTime: now}
	canonical := FileFingerprint{Size: 500, ModTime: now}

	// No cache yet → invalid
	assert.False(t, tc.Valid(gtf, fasta, canonical))

	// Write cache
	c := cache.New()
	c.AddTranscript(&cache.Transcript{
		ID: "ENST00000001.1", Chrom: "1", Start: 100, End: 200, Strand: 1,
	})
	require.NoError(t, tc.Write(c, gtf, fasta, canonical))

	// Same fingerprints → valid
	assert.True(t, tc.Valid(gtf, fasta, canonical))

	// Different size → stale
	gtfChanged := gtf
	gtfChanged.Size = 9999
	assert.False(t, tc.Valid(gtfChanged, fasta, canonical))

	// Different modtime → stale
	fastaChanged := fasta
	fastaChanged.ModTime = now.Add(time.Hour)
	assert.False(t, tc.Valid(gtf, fastaChanged, canonical))
}

func TestTranscriptCacheClear(t *testing.T) {
	dir := t.TempDir()
	tc := NewTranscriptCache(dir)

	now := time.Now()
	fp := FileFingerprint{Size: 100, ModTime: now}

	c := cache.New()
	c.AddTranscript(&cache.Transcript{
		ID: "ENST00000001.1", Chrom: "1", Start: 100, End: 200, Strand: 1,
	})
	require.NoError(t, tc.Write(c, fp, fp, fp))
	assert.True(t, tc.Valid(fp, fp, fp))

	tc.Clear()
	assert.False(t, tc.Valid(fp, fp, fp))
}
