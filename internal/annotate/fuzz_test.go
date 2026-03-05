package annotate

import "testing"

// FuzzParseVariantSpec fuzzes the variant specification parser.
// It should never panic on arbitrary input strings.
func FuzzParseVariantSpec(f *testing.F) {
	seeds := []string{
		// Genomic formats.
		"12:25245350:C:A",
		"chr12:25245350:C>A",
		"12-25245350-C-A",
		"X:100:A:T",
		"chrM:1000:G:C",
		// Protein formats.
		"KRAS G12C",
		"KRAS p.G12C",
		"KRAS p.Gly12Cys",
		"TP53 R248W",
		"BRAF V600E",
		"EGFR L858R",
		"PIK3CA H1047R",
		// Stop codon.
		"TP53 R196*",
		"TP53 p.R196*",
		"TP53 p.R196X",
		// HGVSc formats.
		"KRAS c.35G>T",
		"ENST00000311936:c.35G>T",
		"BRCA1 c.5266dupC",
		// Edge cases.
		"",
		" ",
		":",
		"12:",
		"::",
		"KRAS",
		"p.G12C",
		"c.35G>T",
		"12:abc:C:A",
		"chr99:100:A:T",
	}

	for _, s := range seeds {
		f.Add(s)
	}

	f.Fuzz(func(t *testing.T, input string) {
		// ParseVariantSpec should return an error for invalid input, never panic.
		spec, err := ParseVariantSpec(input)
		if err != nil {
			return
		}
		// If parsing succeeded, exercise the fields.
		_ = spec.Type
		_ = spec.Chrom
		_ = spec.Pos
		_ = spec.Ref
		_ = spec.Alt
		_ = spec.GeneName
		_ = spec.RefAA
		_ = spec.Position
		_ = spec.AltAA
		_ = spec.TranscriptID
		_ = spec.CDSChange
	})
}
