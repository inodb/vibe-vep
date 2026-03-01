// Package vcf provides VCF file parsing functionality.
package vcf

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

// Parser reads variants from a VCF file.
type Parser struct {
	reader      *bufio.Reader
	file        *os.File
	gzipReader  *gzip.Reader
	lineNumber  int
	header      []string
	sampleNames []string // sample names from #CHROM header line
}

// NewParser creates a new VCF parser for the given file.
// Supports both plain VCF and gzipped VCF (.vcf.gz) files.
func NewParser(path string) (*Parser, error) {
	if path == "-" {
		return NewParserFromReader(os.Stdin)
	}

	file, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open vcf file: %w", err)
	}

	p := &Parser{file: file}

	// Check for gzip magic bytes
	buf := make([]byte, 2)
	_, err = file.Read(buf)
	if err != nil {
		file.Close()
		return nil, fmt.Errorf("read vcf header: %w", err)
	}

	// Seek back to beginning
	_, err = file.Seek(0, 0)
	if err != nil {
		file.Close()
		return nil, fmt.Errorf("seek vcf file: %w", err)
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

// parseHeader reads and stores VCF header lines.
func (p *Parser) parseHeader() error {
	for {
		line, err := p.reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
			return fmt.Errorf("read header: %w", err)
		}
		p.lineNumber++

		line = strings.TrimRight(line, "\r\n")

		if strings.HasPrefix(line, "##") {
			p.header = append(p.header, line)
			continue
		}

		if strings.HasPrefix(line, "#CHROM") {
			p.header = append(p.header, line)
			// Extract sample names from columns after FORMAT (index 9+)
			fields := strings.Split(line, "\t")
			if len(fields) > 9 {
				p.sampleNames = fields[9:]
			}
			return nil
		}

		// Non-header line encountered without #CHROM
		return &ParseError{
			Line:    p.lineNumber,
			Message: "expected #CHROM header line",
		}
	}

	return &ParseError{
		Line:    p.lineNumber,
		Message: "no #CHROM header line found",
	}
}

// Next reads the next variant from the VCF file.
// Returns nil, nil when there are no more variants.
func (p *Parser) Next() (*Variant, error) {
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

	return p.parseLine(line)
}

// parseLine parses a single VCF data line into a Variant.
func (p *Parser) parseLine(line string) (*Variant, error) {
	fields := strings.Split(line, "\t")
	if len(fields) < 8 {
		return nil, &ParseError{
			Line:    p.lineNumber,
			Message: fmt.Sprintf("expected at least 8 columns, found %d", len(fields)),
		}
	}

	pos, err := strconv.ParseInt(fields[1], 10, 64)
	if err != nil {
		return nil, &ParseError{
			Line:    p.lineNumber,
			Message: fmt.Sprintf("invalid position: %s", fields[1]),
		}
	}

	qual := 0.0
	if fields[5] != "." {
		qual, _ = strconv.ParseFloat(fields[5], 64)
	}

	v := &Variant{
		Chrom:  fields[0],
		Pos:    pos,
		ID:     fields[2],
		Ref:    fields[3],
		Alt:    fields[4],
		Qual:   qual,
		Filter: fields[6],
		Info:   parseInfo(fields[7]),
	}

	// Capture FORMAT + sample columns if present
	if len(fields) > 8 {
		v.SampleColumns = strings.Join(fields[8:], "\t")
	}

	return v, nil
}

// parseInfo parses the INFO field into a map.
func parseInfo(info string) map[string]interface{} {
	result := make(map[string]interface{})
	if info == "." {
		return result
	}

	for _, kv := range strings.Split(info, ";") {
		parts := strings.SplitN(kv, "=", 2)
		if len(parts) == 2 {
			result[parts[0]] = parts[1]
		} else {
			// Flag-type INFO field
			result[parts[0]] = true
		}
	}

	return result
}

// SplitMultiAllelic splits a multi-allelic variant into separate variants.
func SplitMultiAllelic(v *Variant) []*Variant {
	alts := strings.Split(v.Alt, ",")
	if len(alts) == 1 {
		return []*Variant{v}
	}

	variants := make([]*Variant, len(alts))
	for i, alt := range alts {
		variants[i] = &Variant{
			Chrom:         v.Chrom,
			Pos:           v.Pos,
			ID:            v.ID,
			Ref:           v.Ref,
			Alt:           alt,
			Qual:          v.Qual,
			Filter:        v.Filter,
			Info:          v.Info, // Note: INFO is shared, may need deep copy for some use cases
			SampleColumns: v.SampleColumns,
		}
	}

	return variants
}

// Header returns the VCF header lines.
func (p *Parser) Header() []string {
	return p.header
}

// SampleNames returns sample names from the #CHROM header line.
// Returns nil if no sample columns are present.
func (p *Parser) SampleNames() []string {
	return p.sampleNames
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

// ParseError represents an error during VCF parsing with line context.
type ParseError struct {
	Line    int
	Message string
}

func (e *ParseError) Error() string {
	return fmt.Sprintf("vcf parse error at line %d: %s", e.Line, e.Message)
}
