package vcf

import (
	"strings"
	"testing"
)

// FuzzVCFParser fuzzes the VCF parser with arbitrary data lines.
// The parser should never panic on malformed input.
func FuzzVCFParser(f *testing.F) {
	// Seed with valid VCF data.
	seeds := []string{
		"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n12\t25245350\trs121913529\tC\tA\t100\tPASS\tDP=50",
		"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tA\tG\t.\t.\t.",
		"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n7\t55259515\t.\tT\tG,C\t50\tPASS\tDP=100",
		"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr17\t7577539\t.\tG\tA\t.\tPASS\t.",
		// Indels.
		"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tAG\tA\t.\t.\t.",
		"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tA\tACG\t.\t.\t.",
		// Edge cases.
		"",
		"##fileformat=VCFv4.2\n",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
	}

	for _, s := range seeds {
		f.Add(s)
	}

	f.Fuzz(func(t *testing.T, data string) {
		p, err := NewParserFromReader(strings.NewReader(data))
		if err != nil {
			return // invalid header is fine
		}
		// Read all variants until EOF or error — must not panic.
		for {
			v, err := p.Next()
			if err != nil || v == nil {
				break
			}
			// Exercise NormalizeChrom and SplitMultiAllelic.
			_ = v.NormalizeChrom()
			_ = SplitMultiAllelic(v)
		}
	})
}
