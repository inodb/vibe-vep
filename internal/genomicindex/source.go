package genomicindex

import (
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Pre-built keys for Extra map to avoid per-variant string concatenation.
const (
	extraKeyAMScore  = "alphamissense.score"
	extraKeyAMClass  = "alphamissense.class"
	extraKeyCVClnSig = "clinvar.clnsig"
	extraKeyCVRevStat = "clinvar.clnrevstat"
	extraKeyCVClnDN  = "clinvar.clndn"
	extraKeySigMut   = "signal.mutation_status"
	extraKeySigCount = "signal.count_carriers"
	extraKeySigFreq  = "signal.frequency"
)

// GenomicSource is a unified AnnotationSource that combines AlphaMissense,
// ClinVar, and SIGNAL lookups into a single SQLite point query.
type GenomicSource struct {
	store   *Store
	version string
}

// NewSource creates a GenomicSource backed by the given Store.
func NewSource(store *Store, version string) *GenomicSource {
	return &GenomicSource{store: store, version: version}
}

func (s *GenomicSource) Name() string                   { return "" }
func (s *GenomicSource) Version() string                 { return s.version }
func (s *GenomicSource) MatchLevel() annotate.MatchLevel { return annotate.MatchGenomic }

func (s *GenomicSource) Columns() []annotate.ColumnDef {
	return []annotate.ColumnDef{
		// AlphaMissense
		{Name: "alphamissense.score", Description: "Pathogenicity score (0-1)"},
		{Name: "alphamissense.class", Description: "likely_benign/ambiguous/likely_pathogenic"},
		// ClinVar
		{Name: "clinvar.clnsig", Description: "Clinical significance (e.g. Pathogenic, Benign)"},
		{Name: "clinvar.clnrevstat", Description: "Review status (e.g. reviewed_by_expert_panel)"},
		{Name: "clinvar.clndn", Description: "Disease name(s)"},
		// SIGNAL
		{Name: "signal.mutation_status", Description: "Germline mutation status"},
		{Name: "signal.count_carriers", Description: "Number of carriers in SIGNAL cohort"},
		{Name: "signal.frequency", Description: "Overall allele frequency in SIGNAL cohort"},
	}
}

// Annotate performs a single point lookup and distributes results to annotations.
// AlphaMissense scores are only applied to missense annotations (existing behavior).
// ClinVar and SIGNAL data are applied to all annotations.
func (s *GenomicSource) Annotate(v *vcf.Variant, anns []*annotate.Annotation) {
	chrom := v.NormalizeChrom()
	pos, ref, alt := NormalizeAlleles(v.Pos, v.Ref, v.Alt)

	r, ok := s.store.Lookup(chrom, pos, ref, alt)
	if !ok {
		return
	}

	hasAM := r.AMScore > 0
	hasCV := r.CVClnSig != ""
	hasSig := r.SigMutStatus != ""

	for _, ann := range anns {
		// AlphaMissense: missense only
		if hasAM && isMissense(ann.Consequence) {
			ann.SetExtraKey(extraKeyAMScore, formatScore(r.AMScore))
			ann.SetExtraKey(extraKeyAMClass, r.AMClass)
		}

		// ClinVar: all annotations
		if hasCV {
			ann.SetExtraKey(extraKeyCVClnSig, r.CVClnSig)
			if r.CVClnRevStat != "" {
				ann.SetExtraKey(extraKeyCVRevStat, r.CVClnRevStat)
			}
			if r.CVClnDN != "" {
				ann.SetExtraKey(extraKeyCVClnDN, r.CVClnDN)
			}
		}

		// SIGNAL: all annotations
		if hasSig {
			ann.SetExtraKey(extraKeySigMut, r.SigMutStatus)
			if r.SigCount != "" {
				ann.SetExtraKey(extraKeySigCount, r.SigCount)
			}
			if r.SigFreq != "" {
				ann.SetExtraKey(extraKeySigFreq, r.SigFreq)
			}
		}
	}
}

// Store returns the underlying Store (for Close).
func (s *GenomicSource) Store() *Store {
	return s.store
}

// formatScore formats a float32 score as a 4-decimal string.
func formatScore(score float32) string {
	return strconv.FormatFloat(float64(score), 'f', 4, 32)
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
