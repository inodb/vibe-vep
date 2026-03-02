// Package cache provides VEP cache loading functionality.
package cache

import (
	"sort"
)

// Cache provides access to VEP transcript data for variant annotation.
type Cache struct {
	// transcripts stores transcripts indexed by chromosome
	transcripts map[string][]*Transcript
	// geneIndex maps gene name to transcripts (built lazily)
	geneIndex map[string][]*Transcript
	// trees stores interval trees indexed by chromosome (built by BuildIndex)
	trees map[string]*IntervalTree
}

// New creates a new empty cache.
func New() *Cache {
	return &Cache{
		transcripts: make(map[string][]*Transcript),
	}
}

// AddTranscript adds a transcript to the cache.
func (c *Cache) AddTranscript(t *Transcript) {
	chrom := t.Chrom
	c.transcripts[chrom] = append(c.transcripts[chrom], t)
	// Invalidate gene index so it gets rebuilt on next use
	c.geneIndex = nil
}

// BuildIndex builds interval trees for all chromosomes for O(log n + k) lookup,
// and pre-computes per-transcript CDS offsets and exonic base counts.
// Should be called after all transcripts are loaded.
func (c *Cache) BuildIndex() {
	c.trees = make(map[string]*IntervalTree, len(c.transcripts))
	for chrom, txs := range c.transcripts {
		for _, t := range txs {
			t.BuildCDSIndex()
		}
		c.trees[chrom] = BuildIntervalTree(txs)
	}
}

// FindTranscripts returns all transcripts that overlap a given genomic position.
func (c *Cache) FindTranscripts(chrom string, pos int64) []*Transcript {
	// Use interval tree if built
	if c.trees != nil {
		if tree, ok := c.trees[chrom]; ok {
			return tree.FindOverlaps(pos)
		}
		return nil
	}

	// Fallback to linear scan
	transcripts, ok := c.transcripts[chrom]
	if !ok {
		return nil
	}

	var result []*Transcript
	for _, t := range transcripts {
		if t.Contains(pos) {
			result = append(result, t)
		}
	}
	return result
}

// GetTranscript returns a specific transcript by ID, or nil if not found.
func (c *Cache) GetTranscript(id string) *Transcript {
	for _, transcripts := range c.transcripts {
		for _, t := range transcripts {
			if t.ID == id {
				return t
			}
		}
	}
	return nil
}

// TranscriptCount returns the total number of transcripts in the cache.
func (c *Cache) TranscriptCount() int {
	count := 0
	for _, transcripts := range c.transcripts {
		count += len(transcripts)
	}
	return count
}

// Chromosomes returns a sorted list of chromosomes in the cache.
func (c *Cache) Chromosomes() []string {
	chroms := make([]string, 0, len(c.transcripts))
	for chrom := range c.transcripts {
		chroms = append(chroms, chrom)
	}
	sort.Strings(chroms)
	return chroms
}

// FindTranscriptsByChrom returns all transcripts for a chromosome.
func (c *Cache) FindTranscriptsByChrom(chrom string) []*Transcript {
	return c.transcripts[chrom]
}

// FindTranscriptsByGene returns all transcripts for a given gene name.
func (c *Cache) FindTranscriptsByGene(geneName string) []*Transcript {
	if c.geneIndex == nil {
		c.buildGeneIndex()
	}
	return c.geneIndex[geneName]
}

// buildGeneIndex constructs the gene name → transcripts index.
func (c *Cache) buildGeneIndex() {
	c.geneIndex = make(map[string][]*Transcript)
	for _, transcripts := range c.transcripts {
		for _, t := range transcripts {
			if t.GeneName != "" {
				c.geneIndex[t.GeneName] = append(c.geneIndex[t.GeneName], t)
			}
		}
	}
}
