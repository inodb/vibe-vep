package oncokb

import (
	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Source wraps a CancerGeneList as an annotate.AnnotationSource.
type Source struct {
	cgl CancerGeneList
}

// NewSource creates an AnnotationSource backed by the given CancerGeneList.
func NewSource(cgl CancerGeneList) *Source {
	return &Source{cgl: cgl}
}

func (s *Source) Name() string    { return "oncokb" }
func (s *Source) Version() string { return "cancerGeneList.tsv" }

func (s *Source) Columns() []annotate.ColumnDef {
	return []annotate.ColumnDef{
		{Name: "gene_type", Description: "Gene classification (ONCOGENE/TSG)"},
	}
}

// Annotate adds OncoKB gene type to annotations whose gene is in the cancer gene list.
func (s *Source) Annotate(v *vcf.Variant, anns []*annotate.Annotation) {
	for _, ann := range anns {
		if ga, ok := s.cgl[ann.GeneName]; ok {
			ann.SetExtra("oncokb", "gene_type", ga.GeneType)
		}
	}
}
