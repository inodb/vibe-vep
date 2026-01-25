// Package cache provides VEP cache loading functionality.
package cache

import (
	"bufio"
	"compress/gzip"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

// Loader loads transcript data from VEP cache files.
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

	// Find all region files
	files, err := filepath.Glob(filepath.Join(chromDir, "*.gz"))
	if err != nil {
		return fmt.Errorf("glob region files: %w", err)
	}

	// Also check for JSON files (test data)
	jsonFiles, err := filepath.Glob(filepath.Join(chromDir, "*.json"))
	if err != nil {
		return fmt.Errorf("glob json files: %w", err)
	}

	// Load gzipped region files
	for _, f := range files {
		if err := l.loadRegionFile(c, f); err != nil {
			return fmt.Errorf("load region file %s: %w", f, err)
		}
	}

	// Load JSON files (test data)
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
			return fmt.Errorf("cache directory not found: %s\nHint: Download VEP cache with: vep_install -a cf -s %s -y %s",
				speciesDir, l.species, l.assembly)
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

// LoadRegion loads transcripts overlapping a specific genomic region.
func (l *Loader) LoadRegion(c *Cache, chrom string, start, end int64) error {
	chromDir := l.chromPath(chrom)

	if _, err := os.Stat(chromDir); os.IsNotExist(err) {
		return nil
	}

	// Find region files that might contain this position
	files, err := l.findRegionFiles(chromDir, start, end)
	if err != nil {
		return err
	}

	for _, f := range files {
		if err := l.loadRegionFile(c, f); err != nil {
			return fmt.Errorf("load region file %s: %w", f, err)
		}
	}

	return nil
}

// speciesPath returns the path to the species/assembly directory.
func (l *Loader) speciesPath() string {
	// VEP cache structure: ~/.vep/homo_sapiens/115_GRCh38/
	// We simplify to: {cacheDir}/{species}/{species}_{assembly}/
	assemblyDir := fmt.Sprintf("%s_%s", l.species, l.assembly)
	return filepath.Join(l.cacheDir, l.species, assemblyDir)
}

// chromPath returns the path to a chromosome directory.
func (l *Loader) chromPath(chrom string) string {
	return filepath.Join(l.speciesPath(), chrom)
}

// findRegionFiles finds VEP cache files that might contain data for the given region.
// VEP cache files are named like "25000001-26000000.gz" (1Mb regions).
func (l *Loader) findRegionFiles(chromDir string, start, end int64) ([]string, error) {
	entries, err := os.ReadDir(chromDir)
	if err != nil {
		return nil, fmt.Errorf("read chromosome directory: %w", err)
	}

	var files []string
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}

		name := entry.Name()

		// Check for JSON files (test data)
		if strings.HasSuffix(name, ".json") {
			files = append(files, filepath.Join(chromDir, name))
			continue
		}

		// Skip non-gz files and variation files
		if !strings.HasSuffix(name, ".gz") || strings.Contains(name, "all_vars") {
			continue
		}

		// Parse region from filename (e.g., "25000001-26000000.gz")
		baseName := strings.TrimSuffix(name, ".gz")
		parts := strings.Split(baseName, "-")
		if len(parts) != 2 {
			continue
		}

		regionStart, err1 := strconv.ParseInt(parts[0], 10, 64)
		regionEnd, err2 := strconv.ParseInt(parts[1], 10, 64)
		if err1 != nil || err2 != nil {
			continue
		}

		// Check if regions overlap
		if start <= regionEnd && end >= regionStart {
			files = append(files, filepath.Join(chromDir, name))
		}
	}

	return files, nil
}

// loadRegionFile loads transcripts from a gzipped VEP cache region file.
func (l *Loader) loadRegionFile(c *Cache, path string) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()

	gz, err := gzip.NewReader(f)
	if err != nil {
		return fmt.Errorf("open gzip reader: %w", err)
	}
	defer gz.Close()

	// Try to parse as JSON first (for test data)
	// In production, this would use Sereal decoder
	var transcripts []*Transcript
	if err := json.NewDecoder(gz).Decode(&transcripts); err != nil {
		// If JSON fails, try Sereal format
		// For now, just return the error since we don't have Sereal support yet
		return fmt.Errorf("decode transcripts (JSON): %w", err)
	}

	for _, t := range transcripts {
		c.AddTranscript(t)
	}

	return nil
}

// loadJSONFile loads transcripts from a JSON file (for testing).
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

// CacheInfo contains metadata about a VEP cache.
type CacheInfo struct {
	Species       string
	Assembly      string
	Version       string
	VersionNumber int
	Chromosomes   []string
}

// ReadInfo reads cache metadata from info.txt.
func (l *Loader) ReadInfo() (*CacheInfo, error) {
	infoPath := filepath.Join(l.speciesPath(), "info.txt")

	f, err := os.Open(infoPath)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, fmt.Errorf("cache info not found: %s", infoPath)
		}
		return nil, err
	}
	defer f.Close()

	info := &CacheInfo{
		Species:  l.species,
		Assembly: l.assembly,
	}

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := scanner.Text()
		parts := strings.SplitN(line, "\t", 2)
		if len(parts) != 2 {
			continue
		}

		key, value := parts[0], parts[1]
		switch key {
		case "version":
			info.Version = value
			if v, err := strconv.Atoi(value); err == nil {
				info.VersionNumber = v
			}
		case "assembly":
			info.Assembly = value
		case "species":
			info.Species = value
		}
	}

	// List available chromosomes
	entries, err := os.ReadDir(l.speciesPath())
	if err == nil {
		for _, entry := range entries {
			if entry.IsDir() {
				info.Chromosomes = append(info.Chromosomes, entry.Name())
			}
		}
		sort.Strings(info.Chromosomes)
	}

	return info, scanner.Err()
}
