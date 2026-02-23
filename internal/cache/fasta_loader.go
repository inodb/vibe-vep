// Package cache provides VEP cache loading functionality.
package cache

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

// FASTALoader loads CDS sequences from GENCODE FASTA files.
type FASTALoader struct {
	path      string
	sequences map[string]string    // transcript_id -> full sequence
	cdsRanges map[string][2]int    // transcript_id -> [cdsStart, cdsEnd] (1-based from header)
}

// NewFASTALoader creates a new FASTA loader.
func NewFASTALoader(path string) *FASTALoader {
	return &FASTALoader{
		path:      path,
		sequences: make(map[string]string),
		cdsRanges: make(map[string][2]int),
	}
}

// Load parses the FASTA file and stores sequences indexed by transcript ID.
func (l *FASTALoader) Load() error {
	f, err := os.Open(l.path)
	if err != nil {
		return fmt.Errorf("open FASTA file: %w", err)
	}
	defer f.Close()

	var reader io.Reader = f

	// Handle gzipped files
	if strings.HasSuffix(l.path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return fmt.Errorf("open gzip reader: %w", err)
		}
		defer gz.Close()
		reader = gz
	}

	return l.parseFASTA(reader)
}

// parseFASTA parses FASTA content.
// GENCODE CDS FASTA headers look like:
// >ENST00000456328.2|ENSG00000290825.1|OTTHUMG00000002860.3|OTTHUMT00000007999.2|DDX11L2-202|DDX11L2|459|UTR5:1-200|CDS:201-459|UTR3:460-1657|
func (l *FASTALoader) parseFASTA(reader io.Reader) error {
	scanner := bufio.NewScanner(reader)
	// Increase buffer size for long sequences
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024) // 10MB max line

	var currentID string
	var currentSeq strings.Builder

	for scanner.Scan() {
		line := scanner.Text()

		if strings.HasPrefix(line, ">") {
			// Save previous sequence
			if currentID != "" && currentSeq.Len() > 0 {
				l.sequences[currentID] = currentSeq.String()
			}

			// Parse new header
			currentID = l.parseHeader(line)
			// Parse and store CDS range from header
			if cdsStart, cdsEnd, ok := parseCDSRange(line); ok {
				l.cdsRanges[currentID] = [2]int{cdsStart, cdsEnd}
			}
			currentSeq.Reset()
		} else {
			// Accumulate sequence
			currentSeq.WriteString(strings.TrimSpace(line))
		}
	}

	// Save last sequence
	if currentID != "" && currentSeq.Len() > 0 {
		l.sequences[currentID] = currentSeq.String()
	}

	if err := scanner.Err(); err != nil {
		return fmt.Errorf("scan FASTA: %w", err)
	}

	return nil
}

// parseHeader extracts transcript ID from FASTA header.
// Handles both GENCODE format and simple Ensembl format.
func (l *FASTALoader) parseHeader(header string) string {
	// Remove leading >
	header = strings.TrimPrefix(header, ">")

	// GENCODE format: ENST00000456328.2|ENSG00000290825.1|...
	// Simple format: ENST00000456328.2 description
	// Also handle: >ENST00000456328

	// First, split on pipe for GENCODE format
	if idx := strings.Index(header, "|"); idx != -1 {
		id := header[:idx]
		return stripVersion(id)
	}

	// Split on space for simple format
	if idx := strings.Index(header, " "); idx != -1 {
		id := header[:idx]
		return stripVersion(id)
	}

	// No delimiter, use whole header (minus version)
	return stripVersion(header)
}

// parseCDSRange extracts CDS start and end positions from a GENCODE FASTA header.
// Header format: >ENST...|...|CDS:90-920|...
// Returns 1-based start and end positions.
func parseCDSRange(header string) (start, end int, ok bool) {
	// Look for CDS:start-end pattern in pipe-delimited fields
	for _, field := range strings.Split(header, "|") {
		field = strings.TrimSpace(field)
		if !strings.HasPrefix(field, "CDS:") {
			continue
		}
		rangeStr := field[4:] // strip "CDS:"
		parts := strings.SplitN(rangeStr, "-", 2)
		if len(parts) != 2 {
			return 0, 0, false
		}
		s, err1 := strconv.Atoi(parts[0])
		e, err2 := strconv.Atoi(parts[1])
		if err1 != nil || err2 != nil {
			return 0, 0, false
		}
		return s, e, true
	}
	return 0, 0, false
}

// GetSequence returns the CDS sequence for a transcript ID.
// If CDS boundaries were parsed from the FASTA header, only the CDS portion is returned.
// The ID should be without version suffix (e.g., "ENST00000311936").
func (l *FASTALoader) GetSequence(transcriptID string) string {
	// Try without version first
	id := stripVersion(transcriptID)
	seq, ok := l.sequences[id]
	if !ok {
		seq, ok = l.sequences[transcriptID]
		if !ok {
			return ""
		}
		id = transcriptID
	}

	// Extract CDS portion if boundaries are known
	if cdsRange, hasCDS := l.cdsRanges[id]; hasCDS {
		start := cdsRange[0] - 1 // Convert 1-based to 0-based
		end := cdsRange[1]
		if start >= 0 && end <= len(seq) && start < end {
			return seq[start:end]
		}
	}

	return seq
}

// SequenceCount returns the number of loaded sequences.
func (l *FASTALoader) SequenceCount() int {
	return len(l.sequences)
}

// HasSequence checks if a sequence exists for the given transcript ID.
func (l *FASTALoader) HasSequence(transcriptID string) bool {
	id := stripVersion(transcriptID)
	_, ok := l.sequences[id]
	if ok {
		return true
	}
	_, ok = l.sequences[transcriptID]
	return ok
}
