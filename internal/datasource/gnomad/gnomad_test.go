package gnomad

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestParseVCFLine(t *testing.T) {
	line := "chr1\t10042\t.\tT\tA\t.\tPASS\tAC=4;AN=1526;AF=0.002621;nhomalt=0"
	entry, chrom, ok := ParseVCFLine(line)
	assert.True(t, ok)
	assert.Equal(t, "1", chrom)
	assert.Equal(t, int64(10042), entry.Pos)
	assert.Equal(t, "T", entry.Ref)
	assert.Equal(t, "A", entry.Alt)
	assert.InDelta(t, 0.002621, entry.AF, 1e-7)
	assert.Equal(t, 4, entry.AC)
	assert.Equal(t, 1526, entry.AN)
	assert.Equal(t, 0, entry.Nhomalt)
}

func TestParseVCFLineNoPASS(t *testing.T) {
	line := "1\t100\t.\tA\tG\t.\tAC0\tAC=0;AN=100;AF=0"
	_, _, ok := ParseVCFLine(line)
	assert.False(t, ok, "should skip non-PASS variants")
}

func TestParseVCFLineNoAF(t *testing.T) {
	line := "1\t100\t.\tA\tG\t.\tPASS\tAC=0;AN=100"
	_, _, ok := ParseVCFLine(line)
	assert.False(t, ok, "should skip variants without AF")
}

func TestParseVCFLineMultiAllelic(t *testing.T) {
	line := "1\t100\t.\tA\tG,T\t.\tPASS\tAC=5,2;AN=1000;AF=0.005,0.002;nhomalt=0,0"
	entry, _, ok := ParseVCFLine(line)
	assert.True(t, ok)
	assert.Equal(t, "G", entry.Alt, "should take first ALT only")
	assert.InDelta(t, 0.005, entry.AF, 1e-7)
	assert.Equal(t, 5, entry.AC)
	assert.Equal(t, 0, entry.Nhomalt)
}

func TestParseVCFLineWithExomeAF(t *testing.T) {
	line := "1\t100\t.\tA\tG\t.\tPASS\tAC=10;AN=2000;AF=0.005;AF_exome=0.004;nhomalt=1"
	entry, _, ok := ParseVCFLine(line)
	assert.True(t, ok)
	assert.InDelta(t, 0.005, entry.AF, 1e-7)
	assert.InDelta(t, 0.004, entry.AFExome, 1e-7)
	assert.Equal(t, 1, entry.Nhomalt)
}

func TestParseVCFLineDotFilter(t *testing.T) {
	line := "1\t100\t.\tA\tG\t.\t.\tAC=5;AN=1000;AF=0.005"
	_, _, ok := ParseVCFLine(line)
	assert.True(t, ok, "should accept '.' filter")
}

func TestNormalizeChrom(t *testing.T) {
	assert.Equal(t, "1", NormalizeChrom("chr1"))
	assert.Equal(t, "1", NormalizeChrom("1"))
	assert.Equal(t, "X", NormalizeChrom("chrX"))
}

func TestFormatAF(t *testing.T) {
	assert.Equal(t, "", FormatAF(0))
	assert.Equal(t, "0.005", FormatAF(0.005))
	assert.Equal(t, "1.23e-05", FormatAF(0.0000123))
}

func TestParseInfoFloat(t *testing.T) {
	info := "AC=4;AN=1526;AF=0.002621;nhomalt=0"
	assert.InDelta(t, 0.002621, parseInfoFloat(info, "AF="), 1e-7)
	assert.Equal(t, float64(0), parseInfoFloat(info, "MISSING="))
}

func TestParseInfoInt(t *testing.T) {
	info := "AC=4;AN=1526;AF=0.002621;nhomalt=3"
	assert.Equal(t, 4, parseInfoInt(info, "AC="))
	assert.Equal(t, 1526, parseInfoInt(info, "AN="))
	assert.Equal(t, 3, parseInfoInt(info, "nhomalt="))
	assert.Equal(t, 0, parseInfoInt(info, "MISSING="))
}
