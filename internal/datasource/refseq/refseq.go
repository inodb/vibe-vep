// Package refseq provides Ensembl transcript to RefSeq accession mappings,
// parsed from GENCODE's per-release metadata.RefSeq file.
//
// The file is a headerless TSV shipped alongside the GENCODE GTF:
//
//	ENST00000269305.4	NM_000546.5	NP_000537.3
//	ENST00000269305.4	NM_001126112.2	NP_001119584.1
//
// Column 0 is the versioned Ensembl transcript ID, column 1 the versioned
// RefSeq mRNA accession, and column 2 (optional) the RefSeq protein
// accession, which is not used here. A transcript may appear on several
// lines when it maps to more than one RefSeq mRNA; VEP reports these as a
// list in refseq_transcript_ids, so file order is preserved.
package refseq

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strings"
)

// Store holds transcript-to-RefSeq mRNA accession mappings.
type Store struct {
	byVersioned map[string][]string // ENST00000269305.4 -> ["NM_000546.5", ...]
	byStripped  map[string][]string // ENST00000269305   -> ["NM_000546.5", ...]
}

// Load reads a GENCODE metadata.RefSeq file (optionally gzipped) and returns
// a populated Store.
func Load(path string) (*Store, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var r io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return nil, fmt.Errorf("gzip: %w", err)
		}
		defer gz.Close()
		r = gz
	}

	return parse(r)
}

func parse(r io.Reader) (*Store, error) {
	s := &Store{
		byVersioned: make(map[string][]string),
		byStripped:  make(map[string][]string),
	}

	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	for scanner.Scan() {
		fields := strings.SplitN(scanner.Text(), "\t", 3)
		if len(fields) < 2 {
			continue
		}
		txID := strings.TrimSpace(fields[0])
		refSeqID := strings.TrimSpace(fields[1])
		// Skip blank rows and any header line that may appear in future releases.
		if refSeqID == "" || !strings.HasPrefix(txID, "ENST") {
			continue
		}

		if !contains(s.byVersioned[txID], refSeqID) {
			s.byVersioned[txID] = append(s.byVersioned[txID], refSeqID)
		}
		stripped := stripVersion(txID)
		if !contains(s.byStripped[stripped], refSeqID) {
			s.byStripped[stripped] = append(s.byStripped[stripped], refSeqID)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scan: %w", err)
	}

	return s, nil
}

// LookupByTranscript returns the RefSeq mRNA accessions for the given
// transcript ID, in GENCODE file order. The versioned ID is tried first so a
// transcript keeps the accessions of its own version; the version-stripped ID
// is the fallback for callers whose GENCODE release differs from the mapping.
func (s *Store) LookupByTranscript(transcriptID string) []string {
	if s == nil {
		return nil
	}
	if ids, ok := s.byVersioned[transcriptID]; ok {
		return ids
	}
	return s.byStripped[stripVersion(transcriptID)]
}

// Map returns the version-stripped transcript ID -> RefSeq accessions
// mapping, for callers that join on unversioned IDs.
//
// Stripping is deliberate: the GTF and this metadata file are located by
// independent globs, so a GENCODE release bump can leave a stale GTF beside a
// new mapping file. Versioned keys would silently drop such a pairing to zero
// matches. Within one release an ENST appears with exactly one version, so
// stripping never merges distinct version rows or perturbs element order.
func (s *Store) Map() map[string][]string {
	if s == nil {
		return nil
	}
	return s.byStripped
}

// Count returns the number of transcripts with at least one RefSeq accession.
func (s *Store) Count() int {
	if s == nil {
		return 0
	}
	return len(s.byVersioned)
}

func contains(ids []string, id string) bool {
	for _, x := range ids {
		if x == id {
			return true
		}
	}
	return false
}

// stripVersion removes the version suffix from an Ensembl ID.
func stripVersion(id string) string {
	if idx := strings.IndexByte(id, '.'); idx >= 0 {
		return id[:idx]
	}
	return id
}
