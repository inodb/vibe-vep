// Package uniprot provides Ensembl transcript to UniProt accession mappings.
package uniprot

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

// Store holds transcript-to-UniProt accession mappings.
type Store struct {
	byTranscript map[string]string // unversioned transcript ID -> UniProt accession
}

// Load reads a TSV mapping file (enst_id\tfinal_uniprot_id) and returns a populated Store.
func Load(path string) (*Store, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	s := &Store{byTranscript: make(map[string]string)}

	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	// Skip header line.
	if !scanner.Scan() {
		return nil, fmt.Errorf("empty mapping file")
	}

	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.SplitN(line, "\t", 3)
		if len(fields) < 2 {
			continue
		}
		txID := strings.TrimSpace(fields[0])
		uniprotID := strings.TrimSpace(fields[1])
		if txID == "" || uniprotID == "" {
			continue
		}
		s.byTranscript[stripVersion(txID)] = uniprotID
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scan: %w", err)
	}

	return s, nil
}

// LookupByTranscript returns the UniProt accession for the given transcript ID (version is stripped).
func (s *Store) LookupByTranscript(transcriptID string) string {
	if s == nil {
		return ""
	}
	return s.byTranscript[stripVersion(transcriptID)]
}

// Count returns the number of transcript mappings loaded.
func (s *Store) Count() int {
	return len(s.byTranscript)
}

// stripVersion removes the version suffix from an Ensembl ID.
func stripVersion(id string) string {
	if idx := strings.IndexByte(id, '.'); idx >= 0 {
		return id[:idx]
	}
	return id
}
