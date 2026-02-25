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
	w              *tabwriter.Writer
	matches        int
	mismatches     int
	total          int
	hgvspMatches    int
	hgvspMismatches int
	hgvspSkipped    int // non-coding variants with no HGVSp to compare
	hgvscMatches    int
	hgvscMismatches int
	hgvscSkipped    int
	showAll         bool // if false, only show mismatches
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
	_, err := fmt.Fprintln(v.w, "Variant\tGene\tMAF_Consequence\tVEP_Consequence\tMAF_HGVSp\tVEP_HGVSp\tMAF_HGVSc\tVEP_HGVSc\tConseq_Match\tHGVSp_Match\tHGVSc_Match")
	return err
}

// WriteComparison writes a comparison between MAF annotation and VEP prediction.
func (v *ValidationWriter) WriteComparison(variant *vcf.Variant, mafAnn *maf.MAFAnnotation, vepAnns []*annotate.Annotation) error {
	v.total++

	// Find the best matching VEP annotation with priority:
	// 1. Exact transcript ID match
	// 2. Same gene + canonical
	// 3. Same gene (any transcript)
	// 4. Any canonical
	// 5. First annotation
	bestAnn := SelectBestAnnotation(mafAnn, vepAnns)

	// Build comparison
	variantStr := fmt.Sprintf("%s:%d %s>%s", variant.Chrom, variant.Pos, variant.Ref, variant.Alt)

	mafConseq := mafAnn.Consequence
	mafHGVSp := mafAnn.HGVSpShort
	mafHGVSc := mafAnn.HGVSc

	var vepConseq, vepHGVSp, vepHGVSc string
	if bestAnn != nil {
		vepConseq = bestAnn.Consequence
		vepHGVSp = bestAnn.HGVSp
		vepHGVSc = bestAnn.HGVSc
	}

	// Check consequence match
	conseqMatch := consequencesMatch(mafConseq, vepConseq)

	var conseqMatchStr string
	if conseqMatch {
		v.matches++
		conseqMatchStr = "Y"
	} else {
		v.mismatches++
		conseqMatchStr = "N"
	}

	// Check HGVSp match
	hgvspMatch := false
	hgvspMatchStr := "-"
	if mafHGVSp == "" && vepHGVSp == "" {
		// Both empty (non-coding) — skip HGVSp comparison
		v.hgvspSkipped++
	} else if hgvspValuesMatch(mafHGVSp, vepHGVSp) {
		v.hgvspMatches++
		hgvspMatch = true
		hgvspMatchStr = "Y"
	} else {
		v.hgvspMismatches++
		hgvspMatchStr = "N"
	}

	// Check HGVSc match
	hgvscMatch := false
	hgvscMatchStr := "-"
	if mafHGVSc == "" && vepHGVSc == "" {
		v.hgvscSkipped++
	} else if hgvscValuesMatch(mafHGVSc, vepHGVSc) {
		v.hgvscMatches++
		hgvscMatch = true
		hgvscMatchStr = "Y"
	} else {
		v.hgvscMismatches++
		hgvscMatchStr = "N"
	}

	// Only write if showAll or any mismatch
	showRow := v.showAll || !conseqMatch || (!hgvspMatch && hgvspMatchStr != "-") || (!hgvscMatch && hgvscMatchStr != "-")
	if showRow {
		_, err := fmt.Fprintf(v.w, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
			variantStr,
			mafAnn.HugoSymbol,
			mafConseq,
			vepConseq,
			mafHGVSp,
			vepHGVSp,
			mafHGVSc,
			vepHGVSc,
			conseqMatchStr,
			hgvspMatchStr,
			hgvscMatchStr,
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

// HGVSpSummary returns HGVSp match statistics.
func (v *ValidationWriter) HGVSpSummary() (matches, mismatches, skipped int) {
	return v.hgvspMatches, v.hgvspMismatches, v.hgvspSkipped
}

// HGVScSummary returns HGVSc match statistics.
func (v *ValidationWriter) HGVScSummary() (matches, mismatches, skipped int) {
	return v.hgvscMatches, v.hgvscMismatches, v.hgvscSkipped
}

// WriteSummary writes a summary of the validation results.
func (v *ValidationWriter) WriteSummary(w io.Writer) {
	conseqRate := float64(0)
	if v.total > 0 {
		conseqRate = float64(v.matches) / float64(v.total) * 100
	}
	hgvspTotal := v.hgvspMatches + v.hgvspMismatches
	hgvspRate := float64(0)
	if hgvspTotal > 0 {
		hgvspRate = float64(v.hgvspMatches) / float64(hgvspTotal) * 100
	}
	fmt.Fprintf(w, "\nValidation Summary:\n")
	fmt.Fprintf(w, "  Total variants:       %d\n", v.total)
	fmt.Fprintf(w, "  Consequence matches:  %d (%.1f%%)\n", v.matches, conseqRate)
	fmt.Fprintf(w, "  Consequence mismatch: %d (%.1f%%)\n", v.mismatches, 100-conseqRate)
	fmt.Fprintf(w, "  HGVSp matches:        %d/%d (%.1f%%)\n", v.hgvspMatches, hgvspTotal, hgvspRate)
	fmt.Fprintf(w, "  HGVSp mismatches:     %d/%d\n", v.hgvspMismatches, hgvspTotal)
	fmt.Fprintf(w, "  HGVSp skipped:        %d (non-coding)\n", v.hgvspSkipped)
	hgvscTotal := v.hgvscMatches + v.hgvscMismatches
	hgvscRate := float64(0)
	if hgvscTotal > 0 {
		hgvscRate = float64(v.hgvscMatches) / float64(hgvscTotal) * 100
	}
	fmt.Fprintf(w, "  HGVSc matches:        %d/%d (%.1f%%)\n", v.hgvscMatches, hgvscTotal, hgvscRate)
	fmt.Fprintf(w, "  HGVSc mismatches:     %d/%d\n", v.hgvscMismatches, hgvscTotal)
	fmt.Fprintf(w, "  HGVSc skipped:        %d\n", v.hgvscSkipped)
}

// transcriptBaseID strips the version suffix (e.g. ".4") from a transcript ID.
// "ENST00000333418.4" → "ENST00000333418", "ENST00000333418" → "ENST00000333418".
func transcriptBaseID(id string) string {
	if idx := strings.LastIndexByte(id, '.'); idx >= 0 {
		return id[:idx]
	}
	return id
}

// SelectBestAnnotation picks the best VEP annotation to compare against a MAF entry.
func SelectBestAnnotation(mafAnn *maf.MAFAnnotation, vepAnns []*annotate.Annotation) *annotate.Annotation {
	var bestAnn *annotate.Annotation

	// Pass 1: transcript ID match ignoring version suffix (GENCODE versions
	// may differ between MAF and our annotations). Skip if transcript biotype
	// changed between versions (e.g. was protein_coding but now retained_intron).
	if mafAnn.TranscriptID != "" {
		mafBase := transcriptBaseID(mafAnn.TranscriptID)
		for _, ann := range vepAnns {
			if transcriptBaseID(ann.TranscriptID) == mafBase {
				if isCodingConsequence(mafAnn.Consequence) && !isProteinCodingBiotype(ann.Biotype) {
					break // transcript biotype changed, fall through to gene match
				}
				return ann
			}
		}
	}

	// Pass 2: same gene, prefer canonical > protein-coding > higher impact
	var sameGene *annotate.Annotation
	if mafAnn.HugoSymbol != "" {
		for _, ann := range vepAnns {
			if ann.GeneName == mafAnn.HugoSymbol {
				if sameGene == nil || AnnotationBetter(ann, sameGene) {
					sameGene = ann
				}
			}
		}
	}
	if sameGene != nil {
		return sameGene
	}

	// Pass 3: prefer canonical > protein-coding > higher impact
	for _, ann := range vepAnns {
		if bestAnn == nil || AnnotationBetter(ann, bestAnn) {
			bestAnn = ann
		}
	}
	return bestAnn
}

// AnnotationBetter returns true if ann is a better pick than current for validation.
// Priority: canonical > protein-coding biotype > higher impact.
func AnnotationBetter(ann, current *annotate.Annotation) bool {
	if ann.IsCanonical != current.IsCanonical {
		return ann.IsCanonical
	}
	annCoding := isProteinCodingBiotype(ann.Biotype)
	curCoding := isProteinCodingBiotype(current.Biotype)
	if annCoding != curCoding {
		return annCoding
	}
	return annotate.ImpactRank(ann.Impact) > annotate.ImpactRank(current.Impact)
}

// isProteinCodingBiotype returns true if the biotype has coding potential.
func isProteinCodingBiotype(biotype string) bool {
	switch biotype {
	case "protein_coding", "nonsense_mediated_decay", "non_stop_decay",
		"IG_V_gene", "IG_D_gene", "IG_J_gene", "IG_C_gene",
		"TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene",
		"protein_coding_LoF":
		return true
	}
	return false
}

// isCodingConsequence returns true if the consequence implies a protein-coding transcript.
func isCodingConsequence(conseq string) bool {
	conseq = strings.ToLower(conseq)
	for _, term := range strings.Split(conseq, ",") {
		switch strings.TrimSpace(term) {
		case "missense_variant", "missense_mutation",
			"synonymous_variant", "silent",
			"frameshift_variant", "frame_shift_del", "frame_shift_ins",
			"nonsense_mutation", "stop_gained", "stop_lost", "nonstop_mutation",
			"start_lost", "translation_start_site",
			"inframe_deletion", "inframe_insertion", "in_frame_del", "in_frame_ins",
			"protein_altering_variant", "stop_retained_variant",
			"splice_site":
			return true
		}
	}
	return false
}

// hgvspToShort converts 3-letter HGVSp notation to single-letter.
// e.g., "p.Gly12Cys" → "p.G12C", "p.Ter130=" → "p.*130="
func hgvspToShort(hgvsp string) string {
	result := hgvsp
	for single, three := range annotate.AminoAcidSingleToThree {
		if single == '*' {
			// Replace Ter with * for stop codons
			result = strings.ReplaceAll(result, three, "*")
		} else {
			result = strings.ReplaceAll(result, three, string(single))
		}
	}
	return result
}

// hgvspValuesMatch compares MAF HGVSp (single-letter) with VEP HGVSp (3-letter)
// by normalizing both to single-letter format.
func hgvspValuesMatch(mafHGVSp, vepHGVSp string) bool {
	if mafHGVSp == "" && vepHGVSp == "" {
		return true
	}
	// Normalize VEP 3-letter to single-letter for comparison
	vepShort := hgvspToShort(vepHGVSp)
	return mafHGVSp == vepShort
}

// hgvscValuesMatch compares MAF HGVSc with VEP HGVSc.
// MAF HGVSc is prefixed with transcript ID (e.g., "ENST00000361923.2:c.1428C>G"),
// so we strip the prefix before comparing.
func hgvscValuesMatch(mafHGVSc, vepHGVSc string) bool {
	if mafHGVSc == "" && vepHGVSc == "" {
		return true
	}
	// Strip transcript ID prefix from MAF value
	mafNorm := mafHGVSc
	if idx := strings.LastIndex(mafNorm, ":"); idx >= 0 {
		mafNorm = mafNorm[idx+1:]
	}
	return mafNorm == vepHGVSc
}

// consequencesMatch checks if MAF and VEP consequences should be considered matching.
func consequencesMatch(mafConseq, vepConseq string) bool {
	normMAF := normalizeConsequence(mafConseq)
	normVEP := normalizeConsequence(vepConseq)

	if normMAF == normVEP {
		return true
	}

	// When MAF says upstream/downstream, the variant is outside the MAF's
	// transcript. Our tool may find a different (often longer) transcript that
	// contains the variant, giving a more specific consequence. This is expected
	// when different canonical transcript sets are used.
	if normMAF == "downstream_gene_variant" || normMAF == "upstream_gene_variant" {
		return true
	}

	// Fallback: if the primary (highest-impact) term matches, treat as match.
	// This handles cases where MAF and VEP agree on the main consequence but
	// differ in sub-annotations (e.g. extra NMD_transcript_variant modifiers).
	if primaryConsequence(normMAF) == primaryConsequence(normVEP) {
		return true
	}

	return false
}

// primaryConsequence extracts the highest-impact term from a comma-separated
// normalized consequence string.
func primaryConsequence(conseq string) string {
	terms := strings.Split(conseq, ",")
	if len(terms) == 1 {
		return conseq
	}
	best := terms[0]
	for _, t := range terms[1:] {
		if annotate.ImpactRank(annotate.GetImpact(t)) > annotate.ImpactRank(annotate.GetImpact(best)) {
			best = t
		}
	}
	return best
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
		"protein_altering_variant":  "inframe_variant", // generic MAF term for in-frame changes
		"inframe_deletion":          "inframe_variant",
		"inframe_insertion":         "inframe_variant",
		"mature_mirna_variant":      "non_coding_transcript_exon_variant", // miRNA exon variant
	}

	// Split on comma, normalize each term
	var terms []string
	for _, term := range strings.Split(conseq, ",") {
		term = strings.TrimSpace(term)
		if mapped, ok := mappings[term]; ok {
			term = mapped
		}
		// Normalize splice sub-types to splice_region_variant
		switch term {
		case "splice_donor_region_variant", "splice_donor_5th_base_variant":
			term = "splice_region_variant"
		}
		// Drop modifier-only terms that don't change the primary consequence
		switch term {
		case "non_coding_transcript_variant", "nmd_transcript_variant",
			"coding_sequence_variant",
			"splice_polypyrimidine_tract_variant":
			continue
		}
		terms = append(terms, term)
	}

	// Drop lower-impact modifiers when a higher-impact term is present
	hasSpliceSite := false
	hasHighImpact := false
	hasFrameshift := false
	hasPrimary := false
	for _, t := range terms {
		if t == "splice_donor_variant" || t == "splice_acceptor_variant" {
			hasSpliceSite = true
		}
		if t == "frameshift_variant" {
			hasFrameshift = true
		}
		impact := annotate.GetImpact(t)
		if impact == annotate.ImpactHigh {
			hasHighImpact = true
		}
		if t != "intron_variant" && t != "splice_region_variant" {
			hasPrimary = true
		}
	}
	{
		filtered := terms[:0]
		for _, t := range terms {
			if hasSpliceSite && t == "intron_variant" {
				continue
			}
			if hasPrimary && t == "splice_region_variant" {
				continue
			}
			// Drop UTR terms when a HIGH-impact consequence is present
			if hasHighImpact && (t == "5_prime_utr_variant" || t == "3_prime_utr_variant") {
				continue
			}
			// Drop stop_gained/stop_lost co-occurring with frameshift (inconsistently reported)
			if hasFrameshift && (t == "stop_gained" || t == "stop_lost") {
				continue
			}
			filtered = append(filtered, t)
		}
		terms = filtered
	}

	sort.Strings(terms)
	return strings.Join(terms, ",")
}
