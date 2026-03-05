package maf

import (
	"strings"
	"testing"
)

// FuzzMAFParser fuzzes the MAF parser with arbitrary data lines.
// The parser should never panic on malformed input.
func FuzzMAFParser(f *testing.F) {
	// Minimal valid MAF: header + one data line.
	header := "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tConsequence\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2"
	seeds := []string{
		header + "\nKRAS\t3845\t.\tGRCh38\t12\t25245350\t25245350\t+\tmissense_variant\tMissense_Mutation\tSNP\tC\tC\tA",
		header + "\nTP53\t7157\t.\tGRCh37\t17\t7577539\t7577539\t+\tstop_gained\tNonsense_Mutation\tSNP\tG\tG\tA",
		// Deletion.
		header + "\nBRCA1\t672\t.\tGRCh38\t17\t41244936\t41244937\t+\tframeshift_variant\tFrame_Shift_Del\tDEL\tAG\tAG\t-",
		// Insertion.
		header + "\nEGFR\t1956\t.\tGRCh38\t7\t55259515\t55259516\t+\tinframe_insertion\tIn_Frame_Ins\tINS\t-\t-\tGGT",
		// Empty.
		"",
		header + "\n",
		// Just header.
		header,
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
			v, ann, err := p.NextWithAnnotation()
			if err != nil || v == nil {
				break
			}
			// Exercise methods on the returned objects.
			_ = v.NormalizeChrom()
			if ann != nil {
				_ = ann.HugoSymbol
				_ = ann.RawFields
			}
		}
	})
}
