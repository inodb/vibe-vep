// Package vcf provides VCF file parsing functionality.
package vcf

// Variant represents a single genomic variant from a VCF file.
type Variant struct {
	Chrom  string                 // Chromosome name (e.g., "12", "chr12")
	Pos    int64                  // 1-based genomic position
	ID     string                 // Variant identifier (e.g., rs ID)
	Ref    string                 // Reference allele
	Alt    string                 // Alternate allele (single allele after splitting)
	Qual   float64                // Quality score
	Filter string                 // Filter status (PASS or filter name)
	Info   map[string]interface{} // INFO field key-value pairs
}

// IsSNV returns true if the variant is a single nucleotide variant.
func (v *Variant) IsSNV() bool {
	return len(v.Ref) == 1 && len(v.Alt) == 1
}

// IsIndel returns true if the variant is an insertion or deletion.
func (v *Variant) IsIndel() bool {
	return len(v.Ref) != len(v.Alt)
}

// IsInsertion returns true if the variant is an insertion.
func (v *Variant) IsInsertion() bool {
	return len(v.Alt) > len(v.Ref)
}

// IsDeletion returns true if the variant is a deletion.
func (v *Variant) IsDeletion() bool {
	return len(v.Ref) > len(v.Alt)
}

// NormalizeChrom returns the chromosome name without "chr" prefix.
func (v *Variant) NormalizeChrom() string {
	if len(v.Chrom) > 3 && v.Chrom[:3] == "chr" {
		return v.Chrom[3:]
	}
	return v.Chrom
}
