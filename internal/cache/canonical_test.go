package cache

import (
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestParseCanonicalOverrides_MSKCC(t *testing.T) {
	// Biomart format with MSKCC at col 8
	header := "hgnc_symbol\tcanon_gene\tcanon_tx\texplanation\tgn_tx\tgn_explanation\tuniprot_tx\tuniprot_explanation\tmskcc_tx\tmskcc_explanation\toncokb_tx\toncokb_explanation\n"
	row1 := "KRAS\tENSG00000133703\tENST00000311936\tensembl\tENST00000256078\tgn\tENST00000256078\tuniprot\tENST00000311936\tmskcc\tENST00000256078\toncokb\n"
	row2 := "TP53\tENSG00000141510\tENST00000269305\tensembl\tENST00000269305\tgn\tENST00000269305\tuniprot\tENST00000413465\tmskcc\tENST00000269305\toncokb\n"
	row3 := "EMPTY\tENSG00000000001\tENST00000000001\tensembl\tENST00000000001\tgn\tENST00000000001\tuniprot\tnan\tmskcc\tENST00000000001\toncokb\n"
	input := header + row1 + row2 + row3

	overrides, err := parseCanonicalOverrides(strings.NewReader(input))
	require.NoError(t, err)

	assert.Equal(t, "ENST00000311936", overrides["KRAS"])
	assert.Equal(t, "ENST00000413465", overrides["TP53"])
	assert.NotContains(t, overrides, "EMPTY", "nan values should be skipped")
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
