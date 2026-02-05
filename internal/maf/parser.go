// Package maf provides MAF (Mutation Annotation Format) file parsing functionality.
package maf

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/vcf"
)

// Standard MAF column names
const (
	ColChromosome       = "Chromosome"
	ColStartPosition    = "Start_Position"
	ColEndPosition      = "End_Position"
	ColReferenceAllele  = "Reference_Allele"
	ColTumorSeqAllele2  = "Tumor_Seq_Allele2"
	ColHugoSymbol       = "Hugo_Symbol"
	ColConsequence      = "Consequence"
	ColHGVSpShort       = "HGVSp_Short"
	ColTranscriptID     = "Transcript_ID"
	ColVariantType      = "Variant_Type"
	ColNCBIBuild        = "NCBI_Build"
)

// ColumnIndices holds the indices of important MAF columns.
type ColumnIndices struct {
	Chromosome      int
	StartPosition   int
	EndPosition     int
	ReferenceAllele int
	TumorSeqAllele2 int
	HugoSymbol      int
	Consequence     int
	HGVSpShort      int
	TranscriptID    int
	VariantType     int
	NCBIBuild       int
}

// MAFAnnotation holds the original MAF annotation data for validation.
type MAFAnnotation struct {
	HugoSymbol   string
	Consequence  string
	HGVSpShort   string
	TranscriptID string
	VariantType  string
	NCBIBuild    string
}

// Parser reads variants from a MAF file.
type Parser struct {
	reader     *bufio.Reader
	file       *os.File
	gzipReader *gzip.Reader
	lineNumber int
	columns    ColumnIndices
	headerLine string
}

// NewParser creates a new MAF parser for the given file.
// Supports both plain MAF and gzipped MAF (.maf.gz) files.
func NewParser(path string) (*Parser, error) {
	if path == "-" {
		return NewParserFromReader(os.Stdin)
	}

	file, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open maf file: %w", err)
	}

	p := &Parser{file: file}

	// Check for gzip magic bytes
	buf := make([]byte, 2)
	_, err = file.Read(buf)
	if err != nil {
		file.Close()
		return nil, fmt.Errorf("read maf header: %w", err)
	}

	// Seek back to beginning
	_, err = file.Seek(0, 0)
	if err != nil {
		file.Close()
		return nil, fmt.Errorf("seek maf file: %w", err)
	}

	// Check for gzip magic number (0x1f, 0x8b)
	if buf[0] == 0x1f && buf[1] == 0x8b {
		p.gzipReader, err = gzip.NewReader(file)
		if err != nil {
			file.Close()
			return nil, fmt.Errorf("create gzip reader: %w", err)
		}
		p.reader = bufio.NewReader(p.gzipReader)
	} else {
		p.reader = bufio.NewReader(file)
	}

	// Parse header
	if err := p.parseHeader(); err != nil {
		p.Close()
		return nil, err
	}

	return p, nil
}

// NewParserFromReader creates a parser from an io.Reader (e.g., stdin).
func NewParserFromReader(r io.Reader) (*Parser, error) {
	p := &Parser{
		reader: bufio.NewReader(r),
	}

	if err := p.parseHeader(); err != nil {
		return nil, err
	}

	return p, nil
}

// parseHeader reads and parses the MAF header line to find column indices.
func (p *Parser) parseHeader() error {
	for {
		line, err := p.reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				return &ParseError{
					Line:    p.lineNumber,
					Message: "no header line found",
				}
			}
			return fmt.Errorf("read header: %w", err)
		}
		p.lineNumber++

		line = strings.TrimRight(line, "\r\n")

		// Skip comment lines (start with #)
		if strings.HasPrefix(line, "#") {
			continue
		}

		// Skip empty lines
		if line == "" {
			continue
		}

		// This should be the header line
		p.headerLine = line
		if err := p.parseColumnIndices(line); err != nil {
			return err
		}
		return nil
	}
}

// parseColumnIndices parses the header line to find column indices.
func (p *Parser) parseColumnIndices(headerLine string) error {
	columns := strings.Split(headerLine, "\t")

	// Initialize all indices to -1 (not found)
	p.columns = ColumnIndices{
		Chromosome:      -1,
		StartPosition:   -1,
		EndPosition:     -1,
		ReferenceAllele: -1,
		TumorSeqAllele2: -1,
		HugoSymbol:      -1,
		Consequence:     -1,
		HGVSpShort:      -1,
		TranscriptID:    -1,
		VariantType:     -1,
		NCBIBuild:       -1,
	}

	for i, col := range columns {
		switch col {
		case ColChromosome:
			p.columns.Chromosome = i
		case ColStartPosition:
			p.columns.StartPosition = i
		case ColEndPosition:
			p.columns.EndPosition = i
		case ColReferenceAllele:
			p.columns.ReferenceAllele = i
		case ColTumorSeqAllele2:
			p.columns.TumorSeqAllele2 = i
		case ColHugoSymbol:
			p.columns.HugoSymbol = i
		case ColConsequence:
			p.columns.Consequence = i
		case ColHGVSpShort:
			p.columns.HGVSpShort = i
		case ColTranscriptID:
			p.columns.TranscriptID = i
		case ColVariantType:
			p.columns.VariantType = i
		case ColNCBIBuild:
			p.columns.NCBIBuild = i
		}
	}

	// Validate required columns
	if p.columns.Chromosome == -1 {
		return &ParseError{
			Line:    p.lineNumber,
			Message: "required column 'Chromosome' not found in header",
		}
	}
	if p.columns.StartPosition == -1 {
		return &ParseError{
			Line:    p.lineNumber,
			Message: "required column 'Start_Position' not found in header",
		}
	}
	if p.columns.ReferenceAllele == -1 {
		return &ParseError{
			Line:    p.lineNumber,
			Message: "required column 'Reference_Allele' not found in header",
		}
	}
	if p.columns.TumorSeqAllele2 == -1 {
		return &ParseError{
			Line:    p.lineNumber,
			Message: "required column 'Tumor_Seq_Allele2' not found in header",
		}
	}

	return nil
}

// Next reads the next variant from the MAF file.
// Returns nil, nil when there are no more variants.
func (p *Parser) Next() (*vcf.Variant, error) {
	line, err := p.reader.ReadString('\n')
	if err != nil {
		if err == io.EOF {
			return nil, nil
		}
		return nil, fmt.Errorf("read variant line: %w", err)
	}
	p.lineNumber++

	line = strings.TrimRight(line, "\r\n")
	if line == "" {
		return p.Next() // Skip empty lines
	}

	// Skip comment lines
	if strings.HasPrefix(line, "#") {
		return p.Next()
	}

	return p.parseLine(line)
}

// NextWithAnnotation reads the next variant along with its MAF annotation data.
// This is useful for validation against existing annotations.
func (p *Parser) NextWithAnnotation() (*vcf.Variant, *MAFAnnotation, error) {
	line, err := p.reader.ReadString('\n')
	if err != nil {
		if err == io.EOF {
			return nil, nil, nil
		}
		return nil, nil, fmt.Errorf("read variant line: %w", err)
	}
	p.lineNumber++

	line = strings.TrimRight(line, "\r\n")
	if line == "" {
		return p.NextWithAnnotation() // Skip empty lines
	}

	// Skip comment lines
	if strings.HasPrefix(line, "#") {
		return p.NextWithAnnotation()
	}

	return p.parseLineWithAnnotation(line)
}

// parseLine parses a single MAF data line into a Variant.
func (p *Parser) parseLine(line string) (*vcf.Variant, error) {
	v, _, err := p.parseLineWithAnnotation(line)
	return v, err
}

// parseLineWithAnnotation parses a single MAF data line into a Variant and MAFAnnotation.
func (p *Parser) parseLineWithAnnotation(line string) (*vcf.Variant, *MAFAnnotation, error) {
	fields := strings.Split(line, "\t")

	// Ensure we have enough columns
	minCols := max(p.columns.Chromosome, p.columns.StartPosition, p.columns.ReferenceAllele, p.columns.TumorSeqAllele2)
	if len(fields) <= minCols {
		return nil, nil, &ParseError{
			Line:    p.lineNumber,
			Message: fmt.Sprintf("expected at least %d columns, found %d", minCols+1, len(fields)),
		}
	}

	pos, err := strconv.ParseInt(fields[p.columns.StartPosition], 10, 64)
	if err != nil {
		return nil, nil, &ParseError{
			Line:    p.lineNumber,
			Message: fmt.Sprintf("invalid position: %s", fields[p.columns.StartPosition]),
		}
	}

	ref := fields[p.columns.ReferenceAllele]
	alt := fields[p.columns.TumorSeqAllele2]

	// Handle MAF deletion convention where ref might be longer and alt is "-"
	if alt == "-" {
		alt = ""
	}
	if ref == "-" {
		ref = ""
	}

	v := &vcf.Variant{
		Chrom:  fields[p.columns.Chromosome],
		Pos:    pos,
		ID:     ".",
		Ref:    ref,
		Alt:    alt,
		Qual:   0,
		Filter: ".",
		Info:   make(map[string]interface{}),
	}

	// Build annotation from available columns
	ann := &MAFAnnotation{}

	if p.columns.HugoSymbol >= 0 && p.columns.HugoSymbol < len(fields) {
		ann.HugoSymbol = fields[p.columns.HugoSymbol]
	}
	if p.columns.Consequence >= 0 && p.columns.Consequence < len(fields) {
		ann.Consequence = fields[p.columns.Consequence]
	}
	if p.columns.HGVSpShort >= 0 && p.columns.HGVSpShort < len(fields) {
		ann.HGVSpShort = fields[p.columns.HGVSpShort]
	}
	if p.columns.TranscriptID >= 0 && p.columns.TranscriptID < len(fields) {
		ann.TranscriptID = fields[p.columns.TranscriptID]
	}
	if p.columns.VariantType >= 0 && p.columns.VariantType < len(fields) {
		ann.VariantType = fields[p.columns.VariantType]
	}
	if p.columns.NCBIBuild >= 0 && p.columns.NCBIBuild < len(fields) {
		ann.NCBIBuild = fields[p.columns.NCBIBuild]
	}

	return v, ann, nil
}

// Header returns the MAF header line.
func (p *Parser) Header() string {
	return p.headerLine
}

// Columns returns the parsed column indices.
func (p *Parser) Columns() ColumnIndices {
	return p.columns
}

// LineNumber returns the current line number being processed.
func (p *Parser) LineNumber() int {
	return p.lineNumber
}

// Close closes the parser and underlying file.
func (p *Parser) Close() error {
	if p.gzipReader != nil {
		p.gzipReader.Close()
	}
	if p.file != nil {
		return p.file.Close()
	}
	return nil
}

// ParseError represents an error during MAF parsing with line context.
type ParseError struct {
	Line    int
	Message string
}

func (e *ParseError) Error() string {
	return fmt.Sprintf("maf parse error at line %d: %s", e.Line, e.Message)
}

// max returns the maximum of the provided integers.
func max(values ...int) int {
	m := values[0]
	for _, v := range values[1:] {
		if v > m {
			m = v
		}
	}
	return m
}
