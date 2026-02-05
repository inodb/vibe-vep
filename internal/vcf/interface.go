// Package vcf provides VCF file parsing functionality.
package vcf

// VariantParser is the interface for parsers that read variants.
// Both VCF and MAF parsers implement this interface.
type VariantParser interface {
	// Next reads the next variant.
	// Returns nil, nil when there are no more variants.
	Next() (*Variant, error)

	// Close closes the parser and releases resources.
	Close() error

	// LineNumber returns the current line number being processed.
	LineNumber() int
}
