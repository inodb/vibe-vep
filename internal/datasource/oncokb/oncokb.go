// Package oncokb provides OncoKB cancer gene list loading and annotation.
package oncokb

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

// Annotation holds OncoKB gene-level annotations.
type Annotation struct {
	HugoSymbol string
	GeneType   string // "ONCOGENE", "TSG", or "ONCOGENE,TSG"
}

// CancerGeneList maps Hugo Symbol to Annotation.
type CancerGeneList map[string]*Annotation

// IsCancerGene returns true if the gene is in the cancer gene list.
func (c CancerGeneList) IsCancerGene(gene string) bool {
	_, ok := c[gene]
	return ok
}

// LoadCancerGeneList loads an OncoKB cancerGeneList.tsv file.
// The TSV must have columns "Hugo Symbol" and "Gene Type" in the header.
func LoadCancerGeneList(path string) (CancerGeneList, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open cancer gene list: %w", err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)

	// Read header to find column indices
	if !scanner.Scan() {
		return nil, fmt.Errorf("cancer gene list: empty file")
	}
	header := strings.Split(scanner.Text(), "\t")

	hugoIdx := -1
	geneTypeIdx := -1
	for i, col := range header {
		switch col {
		case "Hugo Symbol":
			hugoIdx = i
		case "Gene Type":
			geneTypeIdx = i
		}
	}
	if hugoIdx < 0 {
		return nil, fmt.Errorf("cancer gene list: missing 'Hugo Symbol' column")
	}
	if geneTypeIdx < 0 {
		return nil, fmt.Errorf("cancer gene list: missing 'Gene Type' column")
	}

	cgl := make(CancerGeneList)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) <= hugoIdx || len(fields) <= geneTypeIdx {
			continue
		}
		hugo := strings.TrimSpace(fields[hugoIdx])
		geneType := strings.TrimSpace(fields[geneTypeIdx])
		if hugo == "" {
			continue
		}
		cgl[hugo] = &Annotation{
			HugoSymbol: hugo,
			GeneType:   geneType,
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("reading cancer gene list: %w", err)
	}

	return cgl, nil
}
