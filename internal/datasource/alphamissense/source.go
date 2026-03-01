package alphamissense

import (
	"fmt"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Source wraps an AlphaMissense Store as an annotate.AnnotationSource.
type Source struct {
	store *Store
}

// NewSource creates an AnnotationSource backed by the given Store.
func NewSource(store *Store) *Source {
	return &Source{store: store}
}

func (s *Source) Name() string    { return "alphamissense" }
func (s *Source) Version() string { return "2023" }

func (s *Source) Columns() []annotate.ColumnDef {
	return []annotate.ColumnDef{
		{Name: "score", Description: "Pathogenicity score (0-1)"},
		{Name: "class", Description: "likely_benign/ambiguous/likely_pathogenic"},
	}
}

// Annotate adds AlphaMissense scores to missense annotations.
func (s *Source) Annotate(v *vcf.Variant, anns []*annotate.Annotation) {
	chrom := v.NormalizeChrom()
	if len(chrom) > 0 && chrom[0] != 'c' {
		chrom = "chr" + chrom
	}
	for _, ann := range anns {
		if isMissense(ann.Consequence) {
			if r, ok := s.store.Lookup(chrom, v.Pos, v.Ref, v.Alt); ok {
				ann.SetExtra("alphamissense", "score", fmt.Sprintf("%.4f", r.Score))
				ann.SetExtra("alphamissense", "class", r.Class)
			}
		}
	}
}

// Store returns the underlying AlphaMissense store.
func (s *Source) Store() *Store {
	return s.store
}

// isMissense returns true if the consequence includes missense_variant.
func isMissense(consequence string) bool {
	for rest := consequence; rest != ""; {
		term := rest
		if i := strings.IndexByte(rest, ','); i >= 0 {
			term = rest[:i]
			rest = rest[i+1:]
		} else {
			rest = ""
		}
		if term == annotate.ConsequenceMissenseVariant {
			return true
		}
	}
	return false
}
