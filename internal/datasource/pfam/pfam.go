// Package pfam provides PFAM protein domain lookups.
// Data comes from Ensembl BioMart exports and the PFAM database.
package pfam

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// Domain holds PFAM domain metadata.
type Domain struct {
	Accession   string // e.g. "PF07714"
	Name        string // e.g. "Pkinase_Tyr"
	Description string // e.g. "Protein tyrosine kinase"
}

// DomainRange represents a PFAM domain mapping on a transcript.
type DomainRange struct {
	PfamDomainID string // PFAM accession (e.g. "PF07714")
	Start        int    // protein position start
	End          int    // protein position end
}

// Store holds PFAM domain data for lookups.
type Store struct {
	domains     map[string]Domain        // accession -> domain info
	transcripts map[string][]DomainRange // unversioned transcript ID -> domain ranges
}

// NewStore creates an empty Store.
func NewStore() *Store {
	return &Store{
		domains:     make(map[string]Domain),
		transcripts: make(map[string][]DomainRange),
	}
}

// Load reads both pfamA.txt and ensembl_biomart_pfam.txt files and returns a populated Store.
func Load(pfamAPath, biomartPath string) (*Store, error) {
	s := NewStore()

	if err := s.loadPfamA(pfamAPath); err != nil {
		return nil, fmt.Errorf("load pfamA: %w", err)
	}

	if err := s.loadBiomart(biomartPath); err != nil {
		return nil, fmt.Errorf("load biomart pfam: %w", err)
	}

	return s, nil
}

// loadPfamA parses pfamA.txt (TSV: pfamA_acc, pfamA_id, description).
func (s *Store) loadPfamA(path string) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	// Read header
	if !scanner.Scan() {
		return fmt.Errorf("empty pfamA file")
	}

	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) < 3 {
			continue
		}
		acc := strings.TrimSpace(fields[0])
		if acc == "" {
			continue
		}
		s.domains[acc] = Domain{
			Accession:   acc,
			Name:        strings.TrimSpace(fields[1]),
			Description: strings.TrimSpace(fields[2]),
		}
	}

	return scanner.Err()
}

// loadBiomart parses ensembl_biomart_pfam.txt.
// Columns: Gene stable ID, Transcript stable ID, Gene name, Pfam domain ID, Pfam domain start, Pfam domain end
func (s *Store) loadBiomart(path string) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	// Read header
	if !scanner.Scan() {
		return fmt.Errorf("empty biomart pfam file")
	}

	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) < 6 {
			continue
		}

		txID := strings.TrimSpace(fields[1])
		pfamID := strings.TrimSpace(fields[3])
		startStr := strings.TrimSpace(fields[4])
		endStr := strings.TrimSpace(fields[5])

		if txID == "" || pfamID == "" || startStr == "" || endStr == "" {
			continue
		}

		start, err := strconv.Atoi(startStr)
		if err != nil {
			continue
		}
		end, err := strconv.Atoi(endStr)
		if err != nil {
			continue
		}

		// Strip version suffix from transcript ID (e.g. ENST00000288602.11 -> ENST00000288602)
		txID = stripVersion(txID)

		s.transcripts[txID] = append(s.transcripts[txID], DomainRange{
			PfamDomainID: pfamID,
			Start:        start,
			End:          end,
		})
	}

	return scanner.Err()
}

// LookupDomain returns the Domain for a given accession, if found.
func (s *Store) LookupDomain(accession string) (Domain, bool) {
	d, ok := s.domains[accession]
	return d, ok
}

// LookupDomains returns Domain info for multiple accessions.
func (s *Store) LookupDomains(accessions []string) []Domain {
	var result []Domain
	for _, acc := range accessions {
		if d, ok := s.domains[acc]; ok {
			result = append(result, d)
		}
	}
	return result
}

// LookupTranscript returns the PFAM domain ranges for a transcript (unversioned ID).
func (s *Store) LookupTranscript(transcriptID string) []DomainRange {
	return s.transcripts[stripVersion(transcriptID)]
}

// DomainCount returns the number of PFAM domains loaded.
func (s *Store) DomainCount() int {
	return len(s.domains)
}

// TranscriptCount returns the number of transcripts with PFAM domain mappings.
func (s *Store) TranscriptCount() int {
	return len(s.transcripts)
}

// stripVersion removes the version suffix from an Ensembl ID.
// e.g. "ENST00000288602.11" -> "ENST00000288602"
func stripVersion(id string) string {
	if idx := strings.IndexByte(id, '.'); idx >= 0 {
		return id[:idx]
	}
	return id
}
