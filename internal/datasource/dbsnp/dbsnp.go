// Package dbsnp provides dbSNP VCF parsing helpers.
// The actual store is in the genomicindex package (unified SQLite).
package dbsnp

import (
	"strconv"
	"strings"
)

// Entry represents a single dbSNP variant.
type Entry struct {
	Pos int64
	Ref string
	Alt string
	ID  string // RS ID (e.g., "rs12345")
}

// ParseVCFLine parses a single dbSNP VCF data line into an Entry.
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

	id := fields[2]
	ref := fields[3]
	alt := fields[4]

	// Skip entries without RS ID.
	if id == "." || id == "" {
		return Entry{}, "", false
	}

	// Handle multi-allelic: take first ALT only.
	if idx := strings.IndexByte(alt, ','); idx >= 0 {
		alt = alt[:idx]
	}

	entry := Entry{
		Pos: pos,
		Ref: ref,
		Alt: alt,
		ID:  id,
	}

	return entry, chrom, true
}

// NormalizeChrom removes "chr" prefix and converts NCBI RefSeq accessions
// (e.g., "NC_000001.11") to plain chromosome numbers.
func NormalizeChrom(chrom string) string {
	chrom = strings.TrimPrefix(chrom, "chr")

	// dbSNP VCFs use RefSeq accessions like NC_000001.11 for chr1.
	if strings.HasPrefix(chrom, "NC_") {
		return refseqToChrom(chrom)
	}

	return chrom
}

// refseqToChrom converts NCBI RefSeq accession to chromosome name.
// NC_000001.* → "1", NC_000022.* → "22", NC_000023.* → "X", NC_000024.* → "Y".
func refseqToChrom(acc string) string {
	// Strip version suffix.
	base := acc
	if idx := strings.IndexByte(acc, '.'); idx >= 0 {
		base = acc[:idx]
	}

	// NC_000001 through NC_000024.
	if len(base) != 9 || !strings.HasPrefix(base, "NC_") {
		return acc // return as-is if unexpected format
	}

	numStr := base[3:] // "000001"
	n, err := strconv.Atoi(numStr)
	if err != nil || n < 1 || n > 24 {
		return acc
	}

	switch n {
	case 23:
		return "X"
	case 24:
		return "Y"
	default:
		return strconv.Itoa(n)
	}
}
