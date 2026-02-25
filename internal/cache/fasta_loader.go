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
	path       string
	sequences  map[string]string // versioned transcript_id -> full sequence
	cdsRanges  map[string][2]int // versioned transcript_id -> [cdsStart, cdsEnd] (1-based from header)
	baseToFull map[string]string // unversioned ID -> versioned ID (for fallback lookup)
}

// NewFASTALoader creates a new FASTA loader.
func NewFASTALoader(path string) *FASTALoader {
	return &FASTALoader{
		path:       path,
		sequences:  make(map[string]string),
		cdsRanges:  make(map[string][2]int),
		baseToFull: make(map[string]string),
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
			// Build base-to-full mapping for unversioned lookups
			if base := stripVersion(currentID); base != currentID {
				l.baseToFull[base] = currentID
			}
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
		return header[:idx]
	}

	// Split on space for simple format
	if idx := strings.Index(header, " "); idx != -1 {
		return header[:idx]
	}

	return header
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
func (l *FASTALoader) GetSequence(transcriptID string) string {
	id, seq, ok := l.lookupSequence(transcriptID)
	if !ok {
		return ""
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

// GetCDSPlusDownstream returns the CDS sequence plus up to maxExtra bases of
// 3'UTR (downstream of CDS). This is used for scanning past the stop codon
// in frameshift and stop-lost predictions.
func (l *FASTALoader) GetCDSPlusDownstream(transcriptID string, maxExtra int) string {
	id, seq, ok := l.lookupSequence(transcriptID)
	if !ok {
		return ""
	}

	cdsRange, hasCDS := l.cdsRanges[id]
	if !hasCDS {
		return seq
	}

	start := cdsRange[0] - 1 // 0-based
	end := cdsRange[1]       // exclusive end of CDS

	// Extend into 3'UTR by up to maxExtra bases
	extEnd := end + maxExtra
	if extEnd > len(seq) {
		extEnd = len(seq)
	}

	if start >= 0 && start < extEnd {
		return seq[start:extEnd]
	}
	return ""
}

// HasSequence checks if a sequence exists for the given transcript ID.
func (l *FASTALoader) HasSequence(transcriptID string) bool {
	_, _, ok := l.lookupSequence(transcriptID)
	return ok
}

// lookupSequence finds a sequence by transcript ID, trying exact match first,
// then falling back to base ID → versioned ID mapping.
func (l *FASTALoader) lookupSequence(transcriptID string) (id, seq string, ok bool) {
	// Try exact match first
	if seq, ok := l.sequences[transcriptID]; ok {
		return transcriptID, seq, true
	}
	// Try base ID → versioned ID mapping (e.g., "ENST00000311936" → "ENST00000311936.8")
	base := stripVersion(transcriptID)
	if full, ok := l.baseToFull[base]; ok {
		if seq, ok := l.sequences[full]; ok {
			return full, seq, true
		}
	}
	return "", "", false
}
