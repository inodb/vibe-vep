// Package gnomad provides gnomAD VCF parsing helpers.
// The actual store is in the genomicindex package (unified SQLite).
package gnomad

import (
	"strconv"
	"strings"
)

// Version constants for gnomAD releases per assembly.
const (
	VersionGRCh38 = "4.1"
	VersionGRCh37 = "2.1.1"
)

// Entry represents a single gnomAD variant annotation.
type Entry struct {
	Pos     int64
	Ref     string
	Alt     string
	AF      float64 // Overall allele frequency
	AFExome float64 // Exome allele frequency (v2 only, 0 if unavailable)
	AC      int     // Allele count
	AN      int     // Allele number (total alleles)
	Nhomalt int     // Number of homozygous alternate individuals
}

// ParseVCFLine parses a single gnomAD VCF data line into an Entry.
// Returns the entry, normalized chromosome, and whether parsing succeeded.
func ParseVCFLine(line string) (Entry, string, bool) {
	// VCF: CHROM POS ID REF ALT QUAL FILTER INFO ...
	fields := strings.SplitN(line, "\t", 9)
	if len(fields) < 8 {
		return Entry{}, "", false
	}

	chrom := NormalizeChrom(fields[0])
	pos, err := strconv.ParseInt(fields[1], 10, 64)
	if err != nil {
		return Entry{}, "", false
	}

	ref := fields[3]
	alt := fields[4]
	info := fields[7]

	// Handle multi-allelic: take first ALT only.
	if idx := strings.IndexByte(alt, ','); idx >= 0 {
		alt = alt[:idx]
	}

	// Only include PASS variants.
	filter := fields[6]
	if filter != "PASS" && filter != "." {
		return Entry{}, "", false
	}

	af := parseInfoFloat(info, "AF=")
	if af == 0 {
		return Entry{}, "", false
	}

	entry := Entry{
		Pos:     pos,
		Ref:     ref,
		Alt:     alt,
		AF:      af,
		AFExome: parseInfoFloat(info, "AF_exome="),
		AC:      parseInfoInt(info, "AC="),
		AN:      parseInfoInt(info, "AN="),
		Nhomalt: parseInfoInt(info, "nhomalt="),
	}

	return entry, chrom, true
}

// parseInfoFloat extracts a float value from a VCF INFO field.
func parseInfoFloat(info, key string) float64 {
	idx := strings.Index(info, key)
	if idx < 0 {
		return 0
	}
	val := info[idx+len(key):]
	if end := strings.IndexByte(val, ';'); end >= 0 {
		val = val[:end]
	}
	// Handle multi-allelic values: take first.
	if end := strings.IndexByte(val, ','); end >= 0 {
		val = val[:end]
	}
	f, _ := strconv.ParseFloat(val, 64)
	return f
}

// parseInfoInt extracts an integer value from a VCF INFO field.
func parseInfoInt(info, key string) int {
	idx := strings.Index(info, key)
	if idx < 0 {
		return 0
	}
	val := info[idx+len(key):]
	if end := strings.IndexByte(val, ';'); end >= 0 {
		val = val[:end]
	}
	// Handle multi-allelic values: take first.
	if end := strings.IndexByte(val, ','); end >= 0 {
		val = val[:end]
	}
	n, _ := strconv.Atoi(val)
	return n
}

// NormalizeChrom removes "chr" prefix for consistent lookups.
func NormalizeChrom(chrom string) string {
	return strings.TrimPrefix(chrom, "chr")
}

// FormatAF formats an allele frequency value for output.
func FormatAF(v float64) string {
	if v == 0 {
		return ""
	}
	return strconv.FormatFloat(v, 'g', 6, 64)
}
