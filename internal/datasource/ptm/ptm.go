// Package ptm provides post-translational modification lookups.
// Data comes from the genome-nexus-importer PTM export.
package ptm

import (
	"bufio"
	"compress/gzip"
	"encoding/json"
	"fmt"
	"os"
	"strings"
)

// PTM represents a single post-translational modification entry.
type PTM struct {
	UniprotEntry         string   `json:"uniprotEntry"`
	UniprotAccession     string   `json:"uniprotAccession"`
	EnsemblTranscriptIds []string `json:"ensemblTranscriptIds"`
	Position             int      `json:"position"`
	Type                 string   `json:"type"`
	PubmedIds            []string `json:"pubmedIds"`
	Sequence             string   `json:"sequence"`
}

// rawPTM mirrors the JSON structure in the ptm.json.gz file (snake_case keys).
type rawPTM struct {
	UniprotEntry         string   `json:"uniprot_entry"`
	UniprotAccession     string   `json:"uniprot_accession"`
	EnsemblTranscriptIds []string `json:"ensembl_transcript_ids"`
	Position             int      `json:"position"`
	Type                 string   `json:"type"`
	PubmedIds            []string `json:"pubmed_ids"`
	Sequence             string   `json:"sequence"`
}

// Store holds PTM data indexed by unversioned transcript ID.
type Store struct {
	byTranscript map[string][]PTM // unversioned transcript ID -> PTMs
}

// Load reads a gzipped JSONL file of PTM entries and returns a populated Store.
func Load(path string) (*Store, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	gz, err := gzip.NewReader(f)
	if err != nil {
		return nil, fmt.Errorf("gzip open: %w", err)
	}
	defer gz.Close()

	s := &Store{byTranscript: make(map[string][]PTM)}

	scanner := bufio.NewScanner(gz)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	for scanner.Scan() {
		var raw rawPTM
		if err := json.Unmarshal(scanner.Bytes(), &raw); err != nil {
			continue // skip malformed lines
		}
		if len(raw.EnsemblTranscriptIds) == 0 {
			continue // skip entries without transcript mappings
		}

		entry := PTM{
			UniprotEntry:         raw.UniprotEntry,
			UniprotAccession:     raw.UniprotAccession,
			EnsemblTranscriptIds: raw.EnsemblTranscriptIds,
			Position:             raw.Position,
			Type:                 raw.Type,
			PubmedIds:            raw.PubmedIds,
			Sequence:             raw.Sequence,
		}

		for _, txID := range raw.EnsemblTranscriptIds {
			key := stripVersion(txID)
			s.byTranscript[key] = append(s.byTranscript[key], entry)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scan: %w", err)
	}

	return s, nil
}

// LookupByTranscript returns PTMs for the given transcript ID (version is stripped).
func (s *Store) LookupByTranscript(transcriptID string) []PTM {
	if s == nil {
		return nil
	}
	return s.byTranscript[stripVersion(transcriptID)]
}

// TranscriptCount returns the number of transcripts with PTM entries.
func (s *Store) TranscriptCount() int {
	return len(s.byTranscript)
}

// stripVersion removes the version suffix from an Ensembl ID.
func stripVersion(id string) string {
	if idx := strings.IndexByte(id, '.'); idx >= 0 {
		return id[:idx]
	}
	return id
}
