package output

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestDiffWriter_Header(t *testing.T) {
	var buf bytes.Buffer
	dw := NewDiffWriter(&buf,
		[]string{"Hugo_Symbol", "Variant_Classification"},
		[]string{"Hugo_Symbol", "Variant_Classification"},
		[]string{"Hugo_Symbol", "Variant_Classification"},
		false, 0)

	require.NoError(t, dw.WriteHeader())
	require.NoError(t, dw.Flush())

	line := strings.TrimSpace(buf.String())
	assert.Equal(t, "Variant\tL_Hugo_Symbol\tR_Hugo_Symbol\tL_Variant_Classification\tR_Variant_Classification", line)
}

func TestDiffWriter_MatchHidden(t *testing.T) {
	var buf bytes.Buffer
	dw := NewDiffWriter(&buf,
		[]string{"Hugo_Symbol"},
		[]string{"Hugo_Symbol"},
		[]string{"Hugo_Symbol"},
		false, 0) // showAll=false

	dw.WriteHeader()

	// Same values → match → should be hidden
	err := dw.WriteDiff("12:100 A>T",
		map[string]string{"Hugo_Symbol": "KRAS"},
		map[string]string{"Hugo_Symbol": "KRAS"})
	require.NoError(t, err)
	require.NoError(t, dw.Flush())

	lines := strings.Split(strings.TrimSpace(buf.String()), "\n")
	assert.Len(t, lines, 1, "only header should be present")

	total, _, _, diffsShown := dw.Stats()
	assert.Equal(t, 1, total)
	assert.Equal(t, 0, diffsShown)
}

func TestDiffWriter_DiffShown(t *testing.T) {
	var buf bytes.Buffer
	dw := NewDiffWriter(&buf,
		[]string{"Hugo_Symbol"},
		[]string{"Hugo_Symbol"},
		[]string{"Hugo_Symbol"},
		false, 0)

	dw.WriteHeader()

	err := dw.WriteDiff("12:100 A>T",
		map[string]string{"Hugo_Symbol": "KRAS"},
		map[string]string{"Hugo_Symbol": "BRAF"})
	require.NoError(t, err)
	require.NoError(t, dw.Flush())

	lines := strings.Split(strings.TrimSpace(buf.String()), "\n")
	assert.Len(t, lines, 2, "header + 1 diff row")
	assert.Contains(t, lines[1], "KRAS")
	assert.Contains(t, lines[1], "BRAF")
}

func TestDiffWriter_ShowAll(t *testing.T) {
	var buf bytes.Buffer
	dw := NewDiffWriter(&buf,
		[]string{"Hugo_Symbol"},
		[]string{"Hugo_Symbol"},
		[]string{"Hugo_Symbol"},
		true, 0) // showAll=true

	dw.WriteHeader()

	dw.WriteDiff("12:100 A>T",
		map[string]string{"Hugo_Symbol": "KRAS"},
		map[string]string{"Hugo_Symbol": "KRAS"})
	require.NoError(t, dw.Flush())

	lines := strings.Split(strings.TrimSpace(buf.String()), "\n")
	assert.Len(t, lines, 2, "match row should be shown with --all")
}

func TestDiffWriter_MaxDiffs(t *testing.T) {
	var buf bytes.Buffer
	dw := NewDiffWriter(&buf,
		[]string{"Hugo_Symbol"},
		[]string{"Hugo_Symbol"},
		[]string{"Hugo_Symbol"},
		false, 2) // maxDiffs=2

	dw.WriteHeader()

	for i := 0; i < 5; i++ {
		dw.WriteDiff("12:100 A>T",
			map[string]string{"Hugo_Symbol": "KRAS"},
			map[string]string{"Hugo_Symbol": "BRAF"})
	}
	require.NoError(t, dw.Flush())

	lines := strings.Split(strings.TrimSpace(buf.String()), "\n")
	assert.Len(t, lines, 3, "header + 2 diff rows (maxDiffs=2)")

	_, _, _, diffsShown := dw.Stats()
	assert.Equal(t, 2, diffsShown)
}

func TestDiffWriter_Summary(t *testing.T) {
	var buf, summary bytes.Buffer
	dw := NewDiffWriter(&buf,
		[]string{"Hugo_Symbol", "Variant_Classification"},
		[]string{"Hugo_Symbol", "Variant_Classification"},
		[]string{"Hugo_Symbol", "Variant_Classification"},
		false, 0)

	// 2 matches, 1 diff
	dw.WriteDiff("12:100 A>T",
		map[string]string{"Hugo_Symbol": "KRAS", "Variant_Classification": "Missense_Mutation"},
		map[string]string{"Hugo_Symbol": "KRAS", "Variant_Classification": "Missense_Mutation"})
	dw.WriteDiff("12:200 G>C",
		map[string]string{"Hugo_Symbol": "KRAS", "Variant_Classification": "Missense_Mutation"},
		map[string]string{"Hugo_Symbol": "BRAF", "Variant_Classification": "Silent"})
	dw.WriteLeftOnly("12:300 C>T")
	dw.WriteRightOnly("12:400 T>G")

	dw.WriteSummary(&summary, "left.maf", "right.maf", 3, 3)

	out := summary.String()
	assert.Contains(t, out, "left.maf")
	assert.Contains(t, out, "right.maf")
	assert.Contains(t, out, "Shared:")
	assert.Contains(t, out, "Left only:")
	assert.Contains(t, out, "Right only:")
	assert.Contains(t, out, "Hugo_Symbol")
	assert.Contains(t, out, "50.0%")
	// Transition should show KRAS → BRAF
	assert.Contains(t, out, "KRAS")
	assert.Contains(t, out, "BRAF")
}

func TestDiffWriter_LeftRightOnly(t *testing.T) {
	var buf bytes.Buffer
	dw := NewDiffWriter(&buf,
		[]string{"Hugo_Symbol"},
		[]string{"Hugo_Symbol"},
		[]string{"Hugo_Symbol"},
		false, 0)

	dw.WriteLeftOnly("12:100 A>T")
	dw.WriteLeftOnly("12:200 G>C")
	dw.WriteRightOnly("12:300 C>A")

	_, leftOnly, rightOnly, _ := dw.Stats()
	assert.Equal(t, 2, leftOnly)
	assert.Equal(t, 1, rightOnly)
}

func TestDiffWriter_MappedColumns(t *testing.T) {
	var buf bytes.Buffer
	// Simulate comparing Protein_Change (left) vs HGVSp_Short (right)
	dw := NewDiffWriter(&buf,
		[]string{"Protein_Change"},           // display name
		[]string{"Protein_Change"},           // left column name
		[]string{"HGVSp_Short"},              // right column name
		false, 0)

	dw.WriteHeader()
	err := dw.WriteDiff("12:100 A>T",
		map[string]string{"Protein_Change": "p.G12C"},
		map[string]string{"HGVSp_Short": "p.G12C"})
	require.NoError(t, err)
	require.NoError(t, dw.Flush())

	// Should match — no diff shown (showAll=false)
	lines := strings.Split(strings.TrimSpace(buf.String()), "\n")
	assert.Len(t, lines, 1, "match should be hidden")

	total, _, _, _ := dw.Stats()
	assert.Equal(t, 1, total)
}

func TestResolveColumns(t *testing.T) {
	leftHeader := []string{"Hugo_Symbol", "Variant_Classification", "HGVSp_Short", "Custom_L"}
	rightHeader := []string{"Hugo_Symbol", "Variant_Classification", "HGVSp_Short", "Custom_R"}

	t.Run("auto_intersect", func(t *testing.T) {
		display, left, right, warns := ResolveColumns(leftHeader, rightHeader, nil, nil)
		assert.Equal(t, []string{"Hugo_Symbol", "Variant_Classification", "HGVSp_Short"}, display)
		assert.Equal(t, display, left)
		assert.Equal(t, display, right)
		assert.Empty(t, warns)
	})

	t.Run("filter_columns", func(t *testing.T) {
		display, _, _, _ := ResolveColumns(leftHeader, rightHeader,
			[]string{"Hugo_Symbol", "HGVSp_Short"}, nil)
		assert.Equal(t, []string{"Hugo_Symbol", "HGVSp_Short"}, display)
	})

	t.Run("column_map", func(t *testing.T) {
		display, left, right, warns := ResolveColumns(leftHeader, rightHeader,
			nil, map[string]string{"Custom_L": "Custom_R"})
		assert.Contains(t, display, "Custom_L")
		assert.Contains(t, left, "Custom_L")
		assert.Contains(t, right, "Custom_R")
		assert.Empty(t, warns)
	})

	t.Run("missing_column_warning", func(t *testing.T) {
		_, _, _, warns := ResolveColumns(leftHeader, rightHeader,
			[]string{"NonExistent"}, nil)
		assert.Len(t, warns, 1)
		assert.Contains(t, warns[0], "NonExistent")
	})
}

func TestParseColumnMap(t *testing.T) {
	t.Run("valid", func(t *testing.T) {
		m, err := ParseColumnMap("Protein_Change=HGVSp_Short,Effect=Variant_Classification")
		require.NoError(t, err)
		assert.Equal(t, "HGVSp_Short", m["Protein_Change"])
		assert.Equal(t, "Variant_Classification", m["Effect"])
	})

	t.Run("empty", func(t *testing.T) {
		m, err := ParseColumnMap("")
		require.NoError(t, err)
		assert.Nil(t, m)
	})

	t.Run("invalid", func(t *testing.T) {
		_, err := ParseColumnMap("bad_format")
		assert.Error(t, err)
	})
}

func TestNormalizeVariantKey(t *testing.T) {
	assert.Equal(t, "12:100 A>T", NormalizeVariantKey("chr12", "100", "A", "T"))
	assert.Equal(t, "12:100 A>T", NormalizeVariantKey("12", "100", "A", "T"))
	assert.Equal(t, "12:100 >T", NormalizeVariantKey("12", "100", "-", "T"))
}

func TestFormatCount(t *testing.T) {
	assert.Equal(t, "0", formatCount(0))
	assert.Equal(t, "100", formatCount(100))
	assert.Equal(t, "1,000", formatCount(1000))
	assert.Equal(t, "1,000,000", formatCount(1000000))
	assert.Equal(t, "1,052,366", formatCount(1052366))
}

func TestReadMAFFile(t *testing.T) {
	// Create a temporary MAF file
	dir := t.TempDir()
	path := filepath.Join(dir, "test.maf")

	content := strings.Join([]string{
		"#comment line",
		"Hugo_Symbol\tChromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2\tVariant_Classification",
		"KRAS\t12\t25245350\tC\tA\tMissense_Mutation",
		"TP53\t17\t7577539\tG\tA\tMissense_Mutation",
		"KRAS\t12\t25245350\tC\tA\tSilent",
	}, "\n")

	require.NoError(t, os.WriteFile(path, []byte(content), 0644))

	header, variants, keys, err := ReadMAFFile(path)
	require.NoError(t, err)

	assert.Contains(t, header, "Hugo_Symbol")
	assert.Len(t, keys, 2, "two unique variant keys")
	assert.Equal(t, "12:25245350 C>A", keys[0])

	// Duplicate key should have 2 entries
	krasRows := variants["12:25245350 C>A"]
	assert.Len(t, krasRows, 2)
	assert.Equal(t, "KRAS", krasRows[0]["Hugo_Symbol"])
	assert.Equal(t, "Missense_Mutation", krasRows[0]["Variant_Classification"])
	assert.Equal(t, "Silent", krasRows[1]["Variant_Classification"])
}

func TestCompareFiles_Identical(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "test.maf")

	content := strings.Join([]string{
		"Hugo_Symbol\tChromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2\tVariant_Classification",
		"KRAS\t12\t25245350\tC\tA\tMissense_Mutation",
		"TP53\t17\t7577539\tG\tA\tMissense_Mutation",
	}, "\n")
	require.NoError(t, os.WriteFile(path, []byte(content), 0644))

	header, variants, keys, err := ReadMAFFile(path)
	require.NoError(t, err)

	var out, summary bytes.Buffer
	err = CompareFiles(header, header, variants, variants, keys, keys,
		"left.maf", "right.maf", nil, nil, false, 0, nil, &out, &summary)
	require.NoError(t, err)

	// No diffs should be output (only header)
	lines := strings.Split(strings.TrimSpace(out.String()), "\n")
	assert.Len(t, lines, 1, "only header, no diff rows")

	// Summary should show 100% match
	assert.Contains(t, summary.String(), "100.0%")
}

func TestCompareFiles_WithDiffs(t *testing.T) {
	dir := t.TempDir()

	leftContent := strings.Join([]string{
		"Hugo_Symbol\tChromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2\tVariant_Classification",
		"KRAS\t12\t25245350\tC\tA\tMissense_Mutation",
		"TP53\t17\t7577539\tG\tA\tMissense_Mutation",
		"BRAF\t7\t140753336\tA\tT\tMissense_Mutation",
	}, "\n")

	rightContent := strings.Join([]string{
		"Hugo_Symbol\tChromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2\tVariant_Classification",
		"KRAS\t12\t25245350\tC\tA\tSilent",
		"TP53\t17\t7577539\tG\tA\tMissense_Mutation",
		"EGFR\t7\t55249071\tC\tT\tMissense_Mutation",
	}, "\n")

	leftPath := filepath.Join(dir, "left.maf")
	rightPath := filepath.Join(dir, "right.maf")
	require.NoError(t, os.WriteFile(leftPath, []byte(leftContent), 0644))
	require.NoError(t, os.WriteFile(rightPath, []byte(rightContent), 0644))

	lHeader, lVariants, lKeys, err := ReadMAFFile(leftPath)
	require.NoError(t, err)
	rHeader, rVariants, rKeys, err := ReadMAFFile(rightPath)
	require.NoError(t, err)

	var out, summary bytes.Buffer
	err = CompareFiles(lHeader, rHeader, lVariants, rVariants, lKeys, rKeys,
		"left.maf", "right.maf", nil, nil, false, 0, nil, &out, &summary)
	require.NoError(t, err)

	// Should have 1 diff row (KRAS classification changed)
	lines := strings.Split(strings.TrimSpace(out.String()), "\n")
	assert.Len(t, lines, 2, "header + 1 diff")
	assert.Contains(t, lines[1], "Missense_Mutation")
	assert.Contains(t, lines[1], "Silent")

	// Summary should show left-only and right-only
	s := summary.String()
	assert.Contains(t, s, "Left only:")
	assert.Contains(t, s, "Right only:")
}

func TestCompareFiles_ColumnFilter(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "test.maf")

	content := strings.Join([]string{
		"Hugo_Symbol\tChromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2\tVariant_Classification\tHGVSp_Short",
		"KRAS\t12\t25245350\tC\tA\tMissense_Mutation\tp.G12C",
	}, "\n")
	require.NoError(t, os.WriteFile(path, []byte(content), 0644))

	header, variants, keys, err := ReadMAFFile(path)
	require.NoError(t, err)

	var out, summary bytes.Buffer
	// Only compare Hugo_Symbol
	err = CompareFiles(header, header, variants, variants, keys, keys,
		"left.maf", "right.maf", []string{"Hugo_Symbol"}, nil, false, 0, nil, &out, &summary)
	require.NoError(t, err)

	// Header should only have L_Hugo_Symbol and R_Hugo_Symbol
	lines := strings.Split(strings.TrimSpace(out.String()), "\n")
	assert.Contains(t, lines[0], "L_Hugo_Symbol")
	assert.NotContains(t, lines[0], "Variant_Classification")
}
