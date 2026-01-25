// Package cache provides VEP cache loading functionality.
package cache

import (
	"sort"
)

// Cache provides access to VEP transcript data for variant annotation.
type Cache struct {
	// transcripts stores transcripts indexed by chromosome
	transcripts map[string][]*Transcript
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
}

// FindTranscripts returns all transcripts that overlap a given genomic position.
func (c *Cache) FindTranscripts(chrom string, pos int64) []*Transcript {
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
