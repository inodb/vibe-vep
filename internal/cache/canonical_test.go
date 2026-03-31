package cache

import (
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestParseBiomartCanonicals(t *testing.T) {
	// Biomart format: col 0 = hgnc, col 2 = ensembl canonical, col 11 = mskcc canonical
	header := "hgnc_symbol\tcanon_gene\tcanon_tx\texplanation\tgn_tx\tgn_explanation\tuniprot_tx\tuniprot_explanation\toncokb_tx\toncokb_explanation\tother_tx\tmskcc_tx\n"
	row1 := "KRAS\tENSG00000133703\tENST00000311936\tensembl\tENST00000256078\tgn\tENST00000256078\tuniprot\tENST00000256078\toncokb\tother\tENST00000311936\n"
	row2 := "TP53\tENSG00000141510\tENST00000269305\tensembl\tENST00000269305\tgn\tENST00000269305\tuniprot\tENST00000269305\toncokb\tother\tENST00000413465\n"
	row3 := "EMPTY\tENSG00000000001\tENST00000000001\tensembl\tENST00000000001\tgn\tENST00000000001\tuniprot\tENST00000000001\toncokb\tother\tnan\n"
	input := header + row1 + row2 + row3

	mskcc, ensembl, _, err := parseBiomartCanonicals(strings.NewReader(input))
	require.NoError(t, err)

	// MSK canonicals (col 11)
	assert.Equal(t, "ENST00000311936", mskcc["KRAS"])
	assert.Equal(t, "ENST00000413465", mskcc["TP53"])
	assert.NotContains(t, mskcc, "EMPTY", "nan values should be skipped")

	// Ensembl canonicals (col 2)
	assert.Equal(t, "ENST00000311936", ensembl["KRAS"])
	assert.Equal(t, "ENST00000269305", ensembl["TP53"])
	assert.Equal(t, "ENST00000000001", ensembl["EMPTY"])
}

func TestBiomartCanonicals_CancerGenesDiffer(t *testing.T) {
	// Spot-check 11 cancer genes where MSK and Ensembl pick different transcripts.
	// Biomart cols: 0=hgnc, 2=ensembl_canonical, 11=mskcc_canonical
	// Columns 1,3-10 are filler.
	filler := "\tfiller\t%s\tfiller\tfiller\tfiller\tfiller\tfiller\tfiller\tfiller\tfiller\t%s\n"
	header := "hgnc_symbol\tcol1\tcanon_tx\tcol3\tcol4\tcol5\tcol6\tcol7\tcol8\tcol9\tcol10\tmskcc_tx\n"

	genes := []struct {
		gene       string
		ensemblTx  string
		mskTx      string
	}{
		{"ABL1", "ENST00000372348", "ENST00000318560"},
		{"AKT1", "ENST00000554581", "ENST00000349310"},
		{"APC", "ENST00000457016", "ENST00000257430"},
		{"BRCA1", "ENST00000471181", "ENST00000357654"},
		{"BRCA2", "ENST00000544455", "ENST00000380152"},
		{"CHEK2", "ENST00000382580", "ENST00000328354"},
		{"FGFR2", "ENST00000457416", "ENST00000358487"},
		{"FGFR3", "ENST00000340107", "ENST00000260795"},
		{"KRAS", "ENST00000256078", "ENST00000311936"},
		{"MET", "ENST00000318493", "ENST00000397752"},
		{"NF1", "ENST00000358273", "ENST00000356175"},
	}

	var sb strings.Builder
	sb.WriteString(header)
	for _, g := range genes {
		sb.WriteString(g.gene)
		sb.WriteString(strings.Replace(filler, "%s", g.ensemblTx, 1))
		// second %s is mskcc
		line := sb.String()
		_ = line
		// Actually let's build properly
	}

	// Rebuild properly
	var input strings.Builder
	input.WriteString(header)
	for _, g := range genes {
		input.WriteString(g.gene)
		for i := 1; i <= 10; i++ {
			if i == 2 {
				input.WriteByte('\t')
				input.WriteString(g.ensemblTx)
			} else {
				input.WriteString("\tfiller")
			}
		}
		input.WriteByte('\t')
		input.WriteString(g.mskTx)
		input.WriteByte('\n')
	}

	mskcc, ensembl, _, err := parseBiomartCanonicals(strings.NewReader(input.String()))
	require.NoError(t, err)

	for _, g := range genes {
		t.Run(g.gene, func(t *testing.T) {
			assert.Equal(t, g.mskTx, mskcc[g.gene], "MSK canonical for %s", g.gene)
			assert.Equal(t, g.ensemblTx, ensembl[g.gene], "Ensembl canonical for %s", g.gene)
			assert.NotEqual(t, mskcc[g.gene], ensembl[g.gene],
				"%s: MSK and Ensembl should pick different transcripts", g.gene)
		})
	}
}

func TestApplyCanonicalOverrides_Independent(t *testing.T) {
	// Verify IsCanonicalMSK and IsCanonicalEnsembl are set independently.
	// Use KRAS as example: Ensembl picks ENST00000256078, MSK picks ENST00000311936.
	c := New()

	tx1 := &Transcript{
		ID:                 "ENST00000256078.10",
		GeneName:           "KRAS",
		Chrom:              "12",
		Start:              25205246,
		End:                25250929,
		Strand:             -1,
		Biotype:            "protein_coding",
		IsCanonicalMSK:     true, // GTF default
		IsCanonicalEnsembl: true, // GTF tag
	}
	tx2 := &Transcript{
		ID:       "ENST00000311936.8",
		GeneName: "KRAS",
		Chrom:    "12",
		Start:    25205246,
		End:      25250929,
		Strand:   -1,
		Biotype:  "protein_coding",
	}
	c.AddTranscript(tx1)
	c.AddTranscript(tx2)

	loader := &GENCODELoader{
		mskCanonicalOverrides: CanonicalOverrides{
			"KRAS": "ENST00000311936",
		},
		ensCanonicalOverrides: CanonicalOverrides{
			"KRAS": "ENST00000256078",
		},
	}
	loader.applyCanonicalOverrides(c)

	// tx1 (ENST00000256078) should be Ensembl canonical but NOT MSK canonical
	assert.False(t, tx1.IsCanonicalMSK, "ENST00000256078 should NOT be MSK canonical")
	assert.True(t, tx1.IsCanonicalEnsembl, "ENST00000256078 SHOULD be Ensembl canonical")

	// tx2 (ENST00000311936) should be MSK canonical but NOT Ensembl canonical
	assert.True(t, tx2.IsCanonicalMSK, "ENST00000311936 SHOULD be MSK canonical")
	assert.False(t, tx2.IsCanonicalEnsembl, "ENST00000311936 should NOT be Ensembl canonical")
}

func TestApplyCanonicalOverrides_BothAgree(t *testing.T) {
	// For genes like BRAF/TP53/EGFR where both sources agree,
	// the same transcript should have both flags true.
	c := New()

	brafCanonical := &Transcript{
		ID:                 "ENST00000288602.11",
		GeneName:           "BRAF",
		Chrom:              "7",
		Start:              140719327,
		End:                140924929,
		Strand:             -1,
		Biotype:            "protein_coding",
		IsCanonicalMSK:     true,
		IsCanonicalEnsembl: true,
	}
	brafOther := &Transcript{
		ID:       "ENST00000496384.7",
		GeneName: "BRAF",
		Chrom:    "7",
		Start:    140719327,
		End:      140924929,
		Strand:   -1,
		Biotype:  "protein_coding",
	}
	c.AddTranscript(brafCanonical)
	c.AddTranscript(brafOther)

	loader := &GENCODELoader{
		mskCanonicalOverrides: CanonicalOverrides{
			"BRAF": "ENST00000288602",
		},
		ensCanonicalOverrides: CanonicalOverrides{
			"BRAF": "ENST00000288602",
		},
	}
	loader.applyCanonicalOverrides(c)

	assert.True(t, brafCanonical.IsCanonicalMSK, "BRAF canonical should be MSK canonical")
	assert.True(t, brafCanonical.IsCanonicalEnsembl, "BRAF canonical should be Ensembl canonical")
	assert.False(t, brafOther.IsCanonicalMSK, "BRAF other should NOT be MSK canonical")
	assert.False(t, brafOther.IsCanonicalEnsembl, "BRAF other should NOT be Ensembl canonical")
}

func TestParseMSKCCOverrides(t *testing.T) {
	// MSKCC isoform overrides format: gene_name, refseq_id, enst_id, note
	input := "gene_name\trefseq_id\tenst_id\tnote\n" +
		"IKZF1\tNM_006060.4\tENST00000331340.3\t\n" +
		"APC\tNM_000038.5\tENST00000257430.4\t\n" +
		"TP53\tNM_000546.5\tENST00000269305.4\t\n"

	overrides, err := ParseMSKCCOverrides(strings.NewReader(input))
	require.NoError(t, err)

	assert.Len(t, overrides, 3)
	assert.Equal(t, "ENST00000331340", overrides["IKZF1"])
	assert.Equal(t, "ENST00000257430", overrides["APC"])
	assert.Equal(t, "ENST00000269305", overrides["TP53"])
}
