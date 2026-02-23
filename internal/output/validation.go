// Package output provides output formatting for annotations.
package output

import (
	"fmt"
	"io"
	"sort"
	"strings"
	"text/tabwriter"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// ValidationWriter writes comparison output between MAF annotations and VEP predictions.
type ValidationWriter struct {
	w           *tabwriter.Writer
	matches     int
	mismatches  int
	total       int
	showAll     bool // if false, only show mismatches
}

// NewValidationWriter creates a new validation output writer.
func NewValidationWriter(w io.Writer, showAll bool) *ValidationWriter {
	return &ValidationWriter{
		w:       tabwriter.NewWriter(w, 0, 0, 2, ' ', 0),
		showAll: showAll,
	}
}

// WriteHeader writes the validation output header.
func (v *ValidationWriter) WriteHeader() error {
	_, err := fmt.Fprintln(v.w, "Variant\tGene\tMAF_Consequence\tVEP_Consequence\tMAF_HGVSp\tVEP_HGVSp\tMatch")
	return err
}

// WriteComparison writes a comparison between MAF annotation and VEP prediction.
func (v *ValidationWriter) WriteComparison(variant *vcf.Variant, mafAnn *maf.MAFAnnotation, vepAnns []*annotate.Annotation) error {
	v.total++

	// Find the best matching VEP annotation (prefer same transcript, then canonical)
	var bestAnn *annotate.Annotation
	for _, ann := range vepAnns {
		if mafAnn.TranscriptID != "" && strings.HasPrefix(ann.TranscriptID, mafAnn.TranscriptID) {
			bestAnn = ann
			break
		}
		if bestAnn == nil || ann.IsCanonical {
			bestAnn = ann
		}
	}

	// Build comparison
	variantStr := fmt.Sprintf("%s:%d %s>%s", variant.Chrom, variant.Pos, variant.Ref, variant.Alt)

	mafConseq := mafAnn.Consequence
	mafHGVSp := mafAnn.HGVSpShort

	var vepConseq, vepHGVSp string
	if bestAnn != nil {
		vepConseq = bestAnn.Consequence
		vepHGVSp = bestAnn.AminoAcidChange
	}

	// Check for match (consequence match is primary)
	conseqMatch := normalizeConsequence(mafConseq) == normalizeConsequence(vepConseq)

	var matchStr string
	if conseqMatch {
		v.matches++
		matchStr = "Y"
	} else {
		v.mismatches++
		matchStr = "N"
	}

	// Only write if showAll or mismatch
	if v.showAll || !conseqMatch {
		_, err := fmt.Fprintf(v.w, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
			variantStr,
			mafAnn.HugoSymbol,
			mafConseq,
			vepConseq,
			mafHGVSp,
			vepHGVSp,
			matchStr,
		)
		return err
	}

	return nil
}

// Flush flushes the writer.
func (v *ValidationWriter) Flush() error {
	return v.w.Flush()
}

// Summary returns match statistics.
func (v *ValidationWriter) Summary() (total, matches, mismatches int) {
	return v.total, v.matches, v.mismatches
}

// WriteSummary writes a summary of the validation results.
func (v *ValidationWriter) WriteSummary(w io.Writer) {
	matchRate := float64(0)
	if v.total > 0 {
		matchRate = float64(v.matches) / float64(v.total) * 100
	}
	fmt.Fprintf(w, "\nValidation Summary:\n")
	fmt.Fprintf(w, "  Total variants:  %d\n", v.total)
	fmt.Fprintf(w, "  Matches:         %d (%.1f%%)\n", v.matches, matchRate)
	fmt.Fprintf(w, "  Mismatches:      %d (%.1f%%)\n", v.mismatches, 100-matchRate)
}

// normalizeConsequence normalizes consequence terms for comparison.
// MAF and VEP may use slightly different terms for the same consequence.
// Comma-separated terms are individually normalized and sorted for consistent comparison.
func normalizeConsequence(conseq string) string {
	conseq = strings.ToLower(strings.TrimSpace(conseq))

	// Common mappings between MAF Variant_Classification and SO terms
	mappings := map[string]string{
		"missense_mutation":        "missense_variant",
		"nonsense_mutation":        "stop_gained",
		"silent":                   "synonymous_variant",
		"splice_site":              "splice_donor_variant",
		"frame_shift_del":          "frameshift_variant",
		"frame_shift_ins":          "frameshift_variant",
		"in_frame_del":             "inframe_deletion",
		"in_frame_ins":             "inframe_insertion",
		"nonstop_mutation":         "stop_lost",
		"translation_start_site":   "start_lost",
		"3'utr":                    "3_prime_utr_variant",
		"5'utr":                    "5_prime_utr_variant",
		"intron":                   "intron_variant",
		"igr":                      "intergenic_variant",
		"3'flank":                  "downstream_gene_variant",
		"5'flank":                  "upstream_gene_variant",
	}

	// Split on comma, normalize each term, then sort for consistent comparison
	terms := strings.Split(conseq, ",")
	for i, term := range terms {
		term = strings.TrimSpace(term)
		if mapped, ok := mappings[term]; ok {
			terms[i] = mapped
		} else {
			terms[i] = term
		}
	}
	sort.Strings(terms)
	return strings.Join(terms, ",")
}
