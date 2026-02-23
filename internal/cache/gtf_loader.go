// Package cache provides VEP cache loading functionality.
package cache

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
)

// GTFLoader loads transcript data from GENCODE GTF files.
type GTFLoader struct {
	path string
}

// NewGTFLoader creates a new GTF loader.
func NewGTFLoader(path string) *GTFLoader {
	return &GTFLoader{path: path}
}

// Load loads all transcripts from the GTF file into the cache.
func (l *GTFLoader) Load(c *Cache) error {
	return l.loadGTF(c, "")
}

// LoadChromosome loads transcripts for a specific chromosome.
func (l *GTFLoader) LoadChromosome(c *Cache, chrom string) error {
	return l.loadGTF(c, chrom)
}

// loadGTF parses the GTF file and populates the cache.
// If filterChrom is non-empty, only loads that chromosome.
func (l *GTFLoader) loadGTF(c *Cache, filterChrom string) error {
	f, err := os.Open(l.path)
	if err != nil {
		return fmt.Errorf("open GTF file: %w", err)
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

	// Parse GTF and build transcript map
	transcripts, err := l.parseGTF(reader, filterChrom)
	if err != nil {
		return err
	}

	// Add transcripts to cache
	for _, t := range transcripts {
		c.AddTranscript(t)
	}

	return nil
}

// gtfFeature represents a parsed GTF line.
type gtfFeature struct {
	chrom      string
	source     string
	featureType string
	start      int64
	end        int64
	score      string
	strand     string
	phase      string
	attributes map[string]string
}

// parseGTF parses GTF content and returns transcripts.
func (l *GTFLoader) parseGTF(reader io.Reader, filterChrom string) (map[string]*Transcript, error) {
	scanner := bufio.NewScanner(reader)
	// Increase buffer size for long lines
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 1024*1024)

	transcripts := make(map[string]*Transcript)
	exonsByTranscript := make(map[string][]Exon)
	cdsByTranscript := make(map[string][][2]int64) // start, end pairs

	lineNum := 0
	for scanner.Scan() {
		lineNum++
		line := scanner.Text()

		// Skip comments and empty lines
		if strings.HasPrefix(line, "#") || line == "" {
			continue
		}

		feat, err := l.parseLine(line)
		if err != nil {
			continue // Skip malformed lines
		}

		// Filter by chromosome if specified
		if filterChrom != "" && feat.chrom != normalizeChrom(filterChrom) {
			continue
		}

		transcriptID := feat.attributes["transcript_id"]
		if transcriptID == "" {
			continue
		}

		// Strip version suffix for consistent lookup
		transcriptID = stripVersion(transcriptID)

		switch feat.featureType {
		case "transcript":
			t := &Transcript{
				ID:          transcriptID,
				GeneID:      stripVersion(feat.attributes["gene_id"]),
				GeneName:    feat.attributes["gene_name"],
				Chrom:       feat.chrom,
				Start:       feat.start,
				End:         feat.end,
				Strand:      parseStrand(feat.strand),
				Biotype:     feat.attributes["transcript_type"],
				IsCanonical: feat.attributes["tag"] == "Ensembl_canonical" || strings.Contains(feat.attributes["tag"], "Ensembl_canonical"),
			}
			// Check for MANE Select tag
			if strings.Contains(feat.attributes["tag"], "MANE_Select") {
				t.IsMANESelect = true
			}
			transcripts[transcriptID] = t

		case "exon":
			exonNum, _ := strconv.Atoi(feat.attributes["exon_number"])
			exon := Exon{
				Number: exonNum,
				Start:  feat.start,
				End:    feat.end,
				Frame:  -1, // Will be set from CDS
			}
			exonsByTranscript[transcriptID] = append(exonsByTranscript[transcriptID], exon)

		case "CDS":
			// Track CDS regions for this transcript
			cdsByTranscript[transcriptID] = append(cdsByTranscript[transcriptID], [2]int64{feat.start, feat.end})

		case "start_codon":
			if t, ok := transcripts[transcriptID]; ok {
				if t.Strand == 1 {
					// Forward strand: start codon is beginning of CDS
					if t.CDSStart == 0 || feat.start < t.CDSStart {
						t.CDSStart = feat.start
					}
				} else {
					// Reverse strand: start codon is end of CDS (genomically)
					if feat.end > t.CDSEnd {
						t.CDSEnd = feat.end
					}
				}
			}

		case "stop_codon":
			if t, ok := transcripts[transcriptID]; ok {
				if t.Strand == 1 {
					// Forward strand: stop codon is end of CDS
					if feat.end > t.CDSEnd {
						t.CDSEnd = feat.end
					}
				} else {
					// Reverse strand: stop codon is beginning of CDS (genomically)
					if t.CDSStart == 0 || feat.start < t.CDSStart {
						t.CDSStart = feat.start
					}
				}
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scan GTF: %w", err)
	}

	// Assemble transcripts with exons and CDS info
	for id, t := range transcripts {
		exons := exonsByTranscript[id]
		if len(exons) == 0 {
			continue
		}

		// Sort exons by genomic position
		sort.Slice(exons, func(i, j int) bool {
			return exons[i].Start < exons[j].Start
		})

		// Calculate CDS boundaries from CDS features if not set by start/stop codons
		if cdsRegions := cdsByTranscript[id]; len(cdsRegions) > 0 {
			minStart := cdsRegions[0][0]
			maxEnd := cdsRegions[0][1]
			for _, region := range cdsRegions[1:] {
				if region[0] < minStart {
					minStart = region[0]
				}
				if region[1] > maxEnd {
					maxEnd = region[1]
				}
			}
			if t.CDSStart == 0 {
				t.CDSStart = minStart
			}
			if t.CDSEnd == 0 {
				t.CDSEnd = maxEnd
			}
		}

		// Calculate CDS portions of exons and reading frames
		if t.CDSStart > 0 && t.CDSEnd > 0 {
			cdsPosition := int64(0)
			for i := range exons {
				e := &exons[i]
				// Check if exon overlaps CDS
				if e.End >= t.CDSStart && e.Start <= t.CDSEnd {
					e.CDSStart = max(e.Start, t.CDSStart)
					e.CDSEnd = min(e.End, t.CDSEnd)

					// Calculate frame based on CDS position
					if t.Strand == 1 {
						e.Frame = int(cdsPosition % 3)
						cdsPosition += e.CDSEnd - e.CDSStart + 1
					}
				}
			}

			// For reverse strand, calculate frames in reverse order
			if t.Strand == -1 {
				cdsPosition = 0
				for i := len(exons) - 1; i >= 0; i-- {
					e := &exons[i]
					if e.CDSStart > 0 && e.CDSEnd > 0 {
						e.Frame = int(cdsPosition % 3)
						cdsPosition += e.CDSEnd - e.CDSStart + 1
					}
				}
			}
		}

		t.Exons = exons
	}

	return transcripts, nil
}

// parseLine parses a single GTF line.
func (l *GTFLoader) parseLine(line string) (*gtfFeature, error) {
	fields := strings.Split(line, "\t")
	if len(fields) < 9 {
		return nil, fmt.Errorf("invalid GTF line: expected 9 fields, got %d", len(fields))
	}

	start, err := strconv.ParseInt(fields[3], 10, 64)
	if err != nil {
		return nil, fmt.Errorf("parse start: %w", err)
	}

	end, err := strconv.ParseInt(fields[4], 10, 64)
	if err != nil {
		return nil, fmt.Errorf("parse end: %w", err)
	}

	feat := &gtfFeature{
		chrom:       normalizeChrom(fields[0]),
		source:      fields[1],
		featureType: fields[2],
		start:       start,
		end:         end,
		score:       fields[5],
		strand:      fields[6],
		phase:       fields[7],
		attributes:  parseAttributes(fields[8]),
	}

	return feat, nil
}

// parseAttributes parses GTF attribute column.
// Format: key "value"; key "value"; ...
func parseAttributes(attrStr string) map[string]string {
	attrs := make(map[string]string)

	// Split by semicolon
	parts := strings.Split(attrStr, ";")
	for _, part := range parts {
		part = strings.TrimSpace(part)
		if part == "" {
			continue
		}

		// Find the first space to separate key from value
		idx := strings.Index(part, " ")
		if idx == -1 {
			continue
		}

		key := part[:idx]
		value := strings.TrimSpace(part[idx+1:])

		// Remove quotes
		value = strings.Trim(value, "\"")

		attrs[key] = value
	}

	return attrs
}

// parseStrand converts strand string to int8.
func parseStrand(s string) int8 {
	if s == "-" {
		return -1
	}
	return 1
}

// stripVersion removes the version suffix from an Ensembl ID.
// e.g., "ENST00000456328.2" -> "ENST00000456328"
func stripVersion(id string) string {
	if idx := strings.LastIndex(id, "."); idx != -1 {
		return id[:idx]
	}
	return id
}

// normalizeChrom normalizes chromosome names by removing "chr" prefix.
// This ensures consistency between different data sources (GENCODE uses "chr1", VCF/MAF often use "1").
func normalizeChrom(chrom string) string {
	if strings.HasPrefix(chrom, "chr") {
		return chrom[3:]
	}
	return chrom
}

// GENCODELoader combines GTF and FASTA loaders for complete annotation data.
type GENCODELoader struct {
	gtfPath            string
	fastaPath          string
	gtf                *GTFLoader
	fasta              *FASTALoader
	canonicalOverrides CanonicalOverrides
}

// NewGENCODELoader creates a loader for GENCODE GTF + FASTA files.
func NewGENCODELoader(gtfPath, fastaPath string) *GENCODELoader {
	return &GENCODELoader{
		gtfPath:   gtfPath,
		fastaPath: fastaPath,
		gtf:       NewGTFLoader(gtfPath),
	}
}

// SetCanonicalOverrides sets Genome Nexus canonical transcript overrides.
// When applied, for each gene with an override, the matching transcript is marked
// as canonical and other transcripts for that gene are unmarked.
func (l *GENCODELoader) SetCanonicalOverrides(overrides CanonicalOverrides) {
	l.canonicalOverrides = overrides
}

// Load loads all transcripts and sequences into the cache.
func (l *GENCODELoader) Load(c *Cache) error {
	// Load GTF annotations
	if err := l.gtf.Load(c); err != nil {
		return fmt.Errorf("load GTF: %w", err)
	}

	// Apply canonical overrides if set
	if len(l.canonicalOverrides) > 0 {
		l.applyCanonicalOverrides(c)
	}

	// Load FASTA sequences if provided
	if l.fastaPath != "" {
		l.fasta = NewFASTALoader(l.fastaPath)
		if err := l.fasta.Load(); err != nil {
			return fmt.Errorf("load FASTA: %w", err)
		}

		// Attach sequences to transcripts
		for _, chrom := range c.Chromosomes() {
			for _, t := range c.FindTranscriptsByChrom(chrom) {
				if seq := l.fasta.GetSequence(t.ID); seq != "" {
					t.CDSSequence = seq
				}
			}
		}
	}

	return nil
}

// applyCanonicalOverrides applies Genome Nexus canonical transcript overrides.
func (l *GENCODELoader) applyCanonicalOverrides(c *Cache) {
	// Group transcripts by gene name
	geneTranscripts := make(map[string][]*Transcript)
	for _, chrom := range c.Chromosomes() {
		for _, t := range c.FindTranscriptsByChrom(chrom) {
			if t.GeneName != "" {
				geneTranscripts[t.GeneName] = append(geneTranscripts[t.GeneName], t)
			}
		}
	}

	// Apply overrides
	for gene, canonicalID := range l.canonicalOverrides {
		transcripts, ok := geneTranscripts[gene]
		if !ok {
			continue
		}

		// Check if the canonical transcript exists for this gene
		found := false
		for _, t := range transcripts {
			if t.ID == canonicalID {
				found = true
				break
			}
		}

		if !found {
			continue
		}

		// Set the override: mark only the canonical transcript
		for _, t := range transcripts {
			t.IsCanonical = (t.ID == canonicalID)
		}
	}
}

// LoadAll implements TranscriptLoader interface.
func (l *GENCODELoader) LoadAll(c *Cache) error {
	return l.Load(c)
}

// GetSequence returns the CDS sequence for a transcript.
func (l *GENCODELoader) GetSequence(transcriptID string) string {
	if l.fasta == nil {
		return ""
	}
	return l.fasta.GetSequence(transcriptID)
}
