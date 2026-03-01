package alphamissense

import (
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Pre-built keys for Extra map to avoid per-variant string concatenation.
const (
	extraKeyScore = "alphamissense.score"
	extraKeyClass = "alphamissense.class"
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
				ann.SetExtraKey(extraKeyScore, formatScore(r.Score))
				ann.SetExtraKey(extraKeyClass, r.Class)
			}
		}
	}
}

// Store returns the underlying AlphaMissense store.
func (s *Source) Store() *Store {
	return s.store
}

// formatScore formats a float64 score as a 4-decimal string without fmt.Sprintf.
func formatScore(score float64) string {
	return strconv.FormatFloat(score, 'f', 4, 64)
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
