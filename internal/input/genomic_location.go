// Package input provides parsers for streaming variant input formats.
package input

import (
	"encoding/json"
	"fmt"
	"strings"

	"github.com/inodb/vibe-vep/internal/vcf"
)

// GenomicLocation represents a genome-nexus style variant location.
type GenomicLocation struct {
	Chromosome     string `json:"chromosome"`
	Start          int64  `json:"start"`
	End            int64  `json:"end"`
	ReferenceAllele string `json:"referenceAllele"`
	VariantAllele  string `json:"variantAllele"`
}

// ParseGenomicLocation parses a JSON line into a GenomicLocation.
func ParseGenomicLocation(line []byte) (GenomicLocation, error) {
	var gl GenomicLocation
	if err := json.Unmarshal(line, &gl); err != nil {
		return GenomicLocation{}, fmt.Errorf("parse genomic location: %w", err)
	}
	if gl.Chromosome == "" {
		return GenomicLocation{}, fmt.Errorf("missing chromosome field")
	}
	if gl.Start == 0 {
		return GenomicLocation{}, fmt.Errorf("missing start field")
	}
	return gl, nil
}

// ToVariant converts a GenomicLocation to a vcf.Variant.
// Handles MAF-style conventions: "-" alleles become empty strings.
func (gl GenomicLocation) ToVariant() *vcf.Variant {
	ref := gl.ReferenceAllele
	alt := gl.VariantAllele
	if ref == "-" {
		ref = ""
	}
	if alt == "-" {
		alt = ""
	}

	return &vcf.Variant{
		Chrom: strings.TrimPrefix(gl.Chromosome, "chr"),
		Pos:   gl.Start,
		Ref:   ref,
		Alt:   alt,
	}
}

// FormatInput returns a comma-separated string representation for the "input" field.
func (gl GenomicLocation) FormatInput() string {
	return fmt.Sprintf("%s,%d,%d,%s,%s",
		gl.Chromosome, gl.Start, gl.End,
		gl.ReferenceAllele, gl.VariantAllele)
}
