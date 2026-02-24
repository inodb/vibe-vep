// Package cache provides transcript cache loading functionality.
package cache

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

// Loader loads transcript data from JSON cache files (used for testing).
type Loader struct {
	cacheDir string
	species  string
	assembly string
}

// NewLoader creates a new cache loader.
func NewLoader(cacheDir, species, assembly string) *Loader {
	return &Loader{
		cacheDir: cacheDir,
		species:  species,
		assembly: assembly,
	}
}

// Load loads all transcripts for a given chromosome into the cache.
func (l *Loader) Load(c *Cache, chrom string) error {
	chromDir := l.chromPath(chrom)

	// Check if chromosome directory exists
	if _, err := os.Stat(chromDir); os.IsNotExist(err) {
		return nil // No data for this chromosome
	}

	// Find JSON files
	jsonFiles, err := filepath.Glob(filepath.Join(chromDir, "*.json"))
	if err != nil {
		return fmt.Errorf("glob json files: %w", err)
	}

	for _, f := range jsonFiles {
		if err := l.loadJSONFile(c, f); err != nil {
			return fmt.Errorf("load json file %s: %w", f, err)
		}
	}

	return nil
}

// LoadAll loads all chromosomes into the cache.
func (l *Loader) LoadAll(c *Cache) error {
	speciesDir := l.speciesPath()

	entries, err := os.ReadDir(speciesDir)
	if err != nil {
		if os.IsNotExist(err) {
			return fmt.Errorf("cache directory not found: %s", speciesDir)
		}
		return fmt.Errorf("read cache directory: %w", err)
	}

	for _, entry := range entries {
		if entry.IsDir() {
			chrom := entry.Name()
			if err := l.Load(c, chrom); err != nil {
				return fmt.Errorf("load chromosome %s: %w", chrom, err)
			}
		}
	}

	return nil
}

// speciesPath returns the path to the species/assembly directory.
func (l *Loader) speciesPath() string {
	speciesDir := filepath.Join(l.cacheDir, l.species)

	// Try to find the correct version directory by scanning for *_{assembly} pattern
	entries, err := os.ReadDir(speciesDir)
	if err == nil {
		suffix := "_" + l.assembly
		for _, entry := range entries {
			if entry.IsDir() && strings.HasSuffix(entry.Name(), suffix) {
				return filepath.Join(speciesDir, entry.Name())
			}
		}
	}

	// Fallback to legacy format
	assemblyDir := fmt.Sprintf("%s_%s", l.species, l.assembly)
	return filepath.Join(l.cacheDir, l.species, assemblyDir)
}

// chromPath returns the path to a chromosome directory.
func (l *Loader) chromPath(chrom string) string {
	return filepath.Join(l.speciesPath(), chrom)
}

// loadJSONFile loads transcripts from a JSON file.
func (l *Loader) loadJSONFile(c *Cache, path string) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()

	var transcripts []*Transcript
	if err := json.NewDecoder(f).Decode(&transcripts); err != nil {
		return fmt.Errorf("decode json: %w", err)
	}

	for _, t := range transcripts {
		c.AddTranscript(t)
	}

	return nil
}
