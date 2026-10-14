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

	// OncoKBAnnotated marks the genes OncoKB curates, a subset of the list.
	// genome-nexus prefers a curated gene over a merely listed one when two
	// genes overlapping a variant each offer a canonical transcript.
	OncoKBAnnotated bool
}

// CuratedGenes returns the symbols OncoKB curates.
func (c CancerGeneList) CuratedGenes() map[string]bool {
	genes := make(map[string]bool)
	for symbol, a := range c {
		if a.OncoKBAnnotated {
			genes[symbol] = true
		}
	}
	return genes
}

// Genes returns every symbol in the list, curated or not.
func (c CancerGeneList) Genes() map[string]bool {
	genes := make(map[string]bool, len(c))
	for symbol := range c {
		genes[symbol] = true
	}
	return genes
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
	annotatedIdx := -1
	for i, col := range header {
		switch col {
		case "Hugo Symbol":
			hugoIdx = i
		case "Gene Type":
			geneTypeIdx = i
		case "OncoKB Annotated":
			annotatedIdx = i
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
		annotated := false
		if annotatedIdx >= 0 && annotatedIdx < len(fields) {
			annotated = strings.EqualFold(strings.TrimSpace(fields[annotatedIdx]), "Yes")
		}
		cgl[hugo] = &Annotation{
			HugoSymbol:      hugo,
			GeneType:        geneType,
			OncoKBAnnotated: annotated,
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("reading cancer gene list: %w", err)
	}

	return cgl, nil
}
