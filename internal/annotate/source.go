package annotate

import "github.com/inodb/vibe-vep/internal/vcf"

// AnnotationSource adds external data to variant annotations.
type AnnotationSource interface {
	Name() string         // e.g. "alphamissense"
	Version() string      // e.g. "2023"
	Columns() []ColumnDef // columns this source provides
	Annotate(v *vcf.Variant, anns []*Annotation)
}

// ColumnDef describes a column provided by an annotation source.
type ColumnDef struct {
	Name        string // short name, e.g. "score"
	Description string // human-readable description
}

// CoreColumns defines the columns produced by vibe-vep's core prediction.
var CoreColumns = []ColumnDef{
	{Name: "hugo_symbol", Description: "Gene symbol"},
	{Name: "consequence", Description: "SO consequence term"},
	{Name: "variant_classification", Description: "MAF variant classification"},
	{Name: "transcript_id", Description: "Ensembl transcript ID"},
	{Name: "hgvsc", Description: "HGVS coding DNA notation"},
	{Name: "hgvsp", Description: "HGVS protein notation (3-letter)"},
	{Name: "hgvsp_short", Description: "HGVS protein notation (1-letter)"},
}

// SetExtra sets a value in the annotation's Extra map.
func (a *Annotation) SetExtra(source, field, value string) {
	if a.Extra == nil {
		a.Extra = make(map[string]string)
	}
	a.Extra[source+"."+field] = value
}

// GetExtra returns a value from the annotation's Extra map.
func (a *Annotation) GetExtra(source, field string) string {
	if a.Extra == nil {
		return ""
	}
	return a.Extra[source+"."+field]
}
