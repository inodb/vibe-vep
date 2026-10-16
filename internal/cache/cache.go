// Package cache provides VEP cache loading functionality.
package cache

import (
	"sort"
	"strings"
)

// Cache provides access to VEP transcript data for variant annotation.
type Cache struct {
	// transcripts stores transcripts indexed by chromosome
	transcripts map[string][]*Transcript
	// geneIndex maps gene name to transcripts (built lazily)
	geneIndex map[string][]*Transcript
	// geneIDIndex maps Ensembl gene ID (ENSG...) to transcripts
	geneIDIndex map[string][]*Transcript
	// proteinIndex maps Ensembl protein ID (ENSP...) to transcripts
	proteinIndex map[string][]*Transcript
	// entrezIndex maps Entrez gene ID to transcripts
	entrezIndex map[string][]*Transcript
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
	// Invalidate the lookup indexes so they get rebuilt on next use
	c.geneIndex = nil
	c.geneIDIndex = nil
	c.proteinIndex = nil
	c.entrezIndex = nil
}

// BuildIndex builds interval trees for all chromosomes for O(log n + k) lookup,
// and pre-computes per-transcript CDS offsets and exonic base counts.
// Should be called after all transcripts are loaded.
func (c *Cache) BuildIndex() {
	// Build the lookup indexes here rather than leaving them to the first
	// reader: the server answers concurrently, and a lazy build under two
	// simultaneous requests is a data race.
	c.buildLookupIndexes()

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
	return c.FindTranscriptsWithinDistance(chrom, pos, 0)
}

// FindTranscriptsWithinDistance returns all transcripts whose body extended by
// `distance` bp on either side contains pos. Passing distance=0 is equivalent
// to FindTranscripts. Used to emit upstream_gene_variant / downstream_gene_variant
// consequences for variants outside but near a transcript body, matching
// Ensembl VEP's --distance behaviour (default 5000 bp).
func (c *Cache) FindTranscriptsWithinDistance(chrom string, pos, distance int64) []*Transcript {
	if distance < 0 {
		distance = 0
	}

	// Use interval tree if built
	if c.trees != nil {
		if tree, ok := c.trees[chrom]; ok {
			return tree.FindOverlapsWithinDistance(pos, distance)
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
		if pos >= t.Start-distance && pos <= t.End+distance {
			result = append(result, t)
		}
	}
	return result
}

// FindTranscriptsInRange returns all transcripts whose body extended by distance
// overlaps the genomic range [start, end]. Used for multi-position variants
// (deletions, MNPs) so transcripts near either end of the variant are found.
func (c *Cache) FindTranscriptsInRange(chrom string, start, end, distance int64) []*Transcript {
	if distance < 0 {
		distance = 0
	}

	if c.trees != nil {
		if tree, ok := c.trees[chrom]; ok {
			return tree.FindOverlapsInRange(start, end, distance)
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
		if t.Start-distance <= end && t.End+distance >= start {
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

// GetTranscriptByPrefix returns the first transcript whose ID starts with the given prefix
// followed by a dot or end of string. This supports unversioned ID lookups
// (e.g. "ENST00000357654" matches "ENST00000357654.9").
func (c *Cache) GetTranscriptByPrefix(prefix string) *Transcript {
	for _, transcripts := range c.transcripts {
		for _, t := range transcripts {
			if t.ID == prefix || (len(t.ID) > len(prefix) && t.ID[:len(prefix)] == prefix && t.ID[len(prefix)] == '.') {
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
		c.buildLookupIndexes()
	}
	return c.geneIndex[geneName]
}

// FindTranscriptsByGeneID returns all transcripts for an Ensembl gene ID,
// with or without its version suffix.
func (c *Cache) FindTranscriptsByGeneID(geneID string) []*Transcript {
	if c.geneIDIndex == nil {
		c.buildLookupIndexes()
	}
	return c.geneIDIndex[stripIDVersion(geneID)]
}

// FindTranscriptsByProteinID returns all transcripts for an Ensembl protein
// ID, with or without its version suffix.
func (c *Cache) FindTranscriptsByProteinID(proteinID string) []*Transcript {
	if c.proteinIndex == nil {
		c.buildLookupIndexes()
	}
	return c.proteinIndex[stripIDVersion(proteinID)]
}

// FindTranscriptsByEntrezID returns all transcripts for an Entrez gene ID.
func (c *Cache) FindTranscriptsByEntrezID(entrezID string) []*Transcript {
	if c.entrezIndex == nil {
		c.buildLookupIndexes()
	}
	return c.entrezIndex[entrezID]
}

// buildLookupIndexes constructs the gene name, gene ID, protein ID and Entrez
// gene ID indexes in one pass over the transcripts.
func (c *Cache) buildLookupIndexes() {
	c.geneIndex = make(map[string][]*Transcript)
	c.geneIDIndex = make(map[string][]*Transcript)
	c.proteinIndex = make(map[string][]*Transcript)
	c.entrezIndex = make(map[string][]*Transcript)
	for _, transcripts := range c.transcripts {
		for _, t := range transcripts {
			if t.GeneName != "" {
				c.geneIndex[t.GeneName] = append(c.geneIndex[t.GeneName], t)
			}
			if id := stripIDVersion(t.GeneID); id != "" {
				c.geneIDIndex[id] = append(c.geneIDIndex[id], t)
			}
			if id := stripIDVersion(t.ProteinID); id != "" {
				c.proteinIndex[id] = append(c.proteinIndex[id], t)
			}
			if t.EntrezGeneID != "" {
				c.entrezIndex[t.EntrezGeneID] = append(c.entrezIndex[t.EntrezGeneID], t)
			}
		}
	}
}

// stripIDVersion drops the ".N" suffix Ensembl identifiers carry, so a lookup
// matches whether or not the caller supplied one.
func stripIDVersion(id string) string {
	if i := strings.IndexByte(id, '.'); i >= 0 {
		return id[:i]
	}
	return id
}
