package oncokb

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestLoadCancerGeneList(t *testing.T) {
	// Find the cancerGeneList.tsv in the repo root
	path := findCancerGeneList(t)

	cgl, err := LoadCancerGeneList(path)
	require.NoError(t, err)
	require.NotEmpty(t, cgl)

	// Check known cancer genes
	tests := []struct {
		gene     string
		geneType string
	}{
		{"ABL1", "ONCOGENE"},
		{"AKT1", "ONCOGENE"},
		{"TP53", "TSG"},
		{"BRCA1", "TSG"},
		{"KRAS", "ONCOGENE"},
	}
	for _, tt := range tests {
		t.Run(tt.gene, func(t *testing.T) {
			ann, ok := cgl[tt.gene]
			require.True(t, ok, "gene %s should be in cancer gene list", tt.gene)
			assert.Equal(t, tt.gene, ann.HugoSymbol)
			assert.Equal(t, tt.geneType, ann.GeneType)
		})
	}
}

func TestLoadCancerGeneList_NotFound(t *testing.T) {
	_, err := LoadCancerGeneList("/nonexistent/path.tsv")
	assert.Error(t, err)
}

func TestCancerGeneList_IsCancerGene(t *testing.T) {
	cgl := CancerGeneList{
		"TP53": &Annotation{HugoSymbol: "TP53", GeneType: "TSG"},
	}
	assert.True(t, cgl.IsCancerGene("TP53"))
	assert.False(t, cgl.IsCancerGene("UNKNOWN"))
}

func findCancerGeneList(t *testing.T) string {
	t.Helper()
	// Walk up to find repo root
	for _, rel := range []string{
		filepath.Join("..", "..", "..", "cancerGeneList.tsv"),
		"cancerGeneList.tsv",
	} {
		if _, err := os.Stat(rel); err == nil {
			return rel
		}
	}
	t.Skip("cancerGeneList.tsv not found")
	return ""
}
