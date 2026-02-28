package duckdb

import (
	"encoding/gob"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"github.com/inodb/vibe-vep/internal/cache"
)

// TranscriptCache manages gob-serialized transcript data on disk.
// Files are stored alongside the GENCODE source files:
//
//	~/.vibe-vep/{assembly}/transcripts.gob       (serialized transcripts)
//	~/.vibe-vep/{assembly}/transcripts.gob.meta  (source file fingerprints)
type TranscriptCache struct {
	dir string // cache directory (e.g. ~/.vibe-vep/grch38)
}

// NewTranscriptCache creates a transcript cache for the given directory.
func NewTranscriptCache(dir string) *TranscriptCache {
	return &TranscriptCache{dir: dir}
}

func (tc *TranscriptCache) gobPath() string {
	return filepath.Join(tc.dir, "transcripts.gob")
}

func (tc *TranscriptCache) metaPath() string {
	return filepath.Join(tc.dir, "transcripts.gob.meta")
}

// Valid checks whether the cached transcripts match the current source files.
func (tc *TranscriptCache) Valid(gtf, fasta, canonical FileFingerprint) bool {
	meta, err := tc.readMeta()
	if err != nil {
		return false
	}

	checks := []struct{ key, val string }{
		{"gtf_size", strconv.FormatInt(gtf.Size, 10)},
		{"gtf_modtime", gtf.ModTime.UTC().Format(time.RFC3339Nano)},
		{"fasta_size", strconv.FormatInt(fasta.Size, 10)},
		{"fasta_modtime", fasta.ModTime.UTC().Format(time.RFC3339Nano)},
		{"canonical_size", strconv.FormatInt(canonical.Size, 10)},
		{"canonical_modtime", canonical.ModTime.UTC().Format(time.RFC3339Nano)},
	}

	for _, c := range checks {
		if meta[c.key] != c.val {
			return false
		}
	}

	// Verify gob file exists
	if _, err := os.Stat(tc.gobPath()); err != nil {
		return false
	}
	return true
}

// Load reads serialized transcripts from disk into the cache.
func (tc *TranscriptCache) Load(c *cache.Cache) error {
	f, err := os.Open(tc.gobPath())
	if err != nil {
		return fmt.Errorf("open transcript cache: %w", err)
	}
	defer f.Close()

	var data map[string][]*cache.Transcript
	if err := gob.NewDecoder(f).Decode(&data); err != nil {
		return fmt.Errorf("decode transcript cache: %w", err)
	}

	for _, transcripts := range data {
		for _, t := range transcripts {
			c.AddTranscript(t)
		}
	}
	return nil
}

// Write serializes all transcripts from the cache to disk.
func (tc *TranscriptCache) Write(c *cache.Cache, gtf, fasta, canonical FileFingerprint) error {
	// Serialize transcripts
	data := make(map[string][]*cache.Transcript)
	for _, chrom := range c.Chromosomes() {
		data[chrom] = c.FindTranscriptsByChrom(chrom)
	}

	f, err := os.Create(tc.gobPath())
	if err != nil {
		return fmt.Errorf("create transcript cache: %w", err)
	}

	if err := gob.NewEncoder(f).Encode(data); err != nil {
		f.Close()
		os.Remove(tc.gobPath())
		return fmt.Errorf("encode transcript cache: %w", err)
	}
	if err := f.Close(); err != nil {
		return fmt.Errorf("close transcript cache: %w", err)
	}

	// Write metadata
	return tc.writeMeta(gtf, fasta, canonical)
}

// Clear removes the cached transcript files.
func (tc *TranscriptCache) Clear() {
	os.Remove(tc.gobPath())
	os.Remove(tc.metaPath())
}

func (tc *TranscriptCache) writeMeta(gtf, fasta, canonical FileFingerprint) error {
	lines := []string{
		"gtf_size=" + strconv.FormatInt(gtf.Size, 10),
		"gtf_modtime=" + gtf.ModTime.UTC().Format(time.RFC3339Nano),
		"fasta_size=" + strconv.FormatInt(fasta.Size, 10),
		"fasta_modtime=" + fasta.ModTime.UTC().Format(time.RFC3339Nano),
		"canonical_size=" + strconv.FormatInt(canonical.Size, 10),
		"canonical_modtime=" + canonical.ModTime.UTC().Format(time.RFC3339Nano),
		"created_at=" + time.Now().UTC().Format(time.RFC3339),
		"",
	}
	return os.WriteFile(tc.metaPath(), []byte(strings.Join(lines, "\n")), 0644)
}

func (tc *TranscriptCache) readMeta() (map[string]string, error) {
	data, err := os.ReadFile(tc.metaPath())
	if err != nil {
		return nil, err
	}

	meta := make(map[string]string)
	for _, line := range strings.Split(string(data), "\n") {
		if k, v, ok := strings.Cut(line, "="); ok {
			meta[k] = v
		}
	}
	return meta, nil
}
