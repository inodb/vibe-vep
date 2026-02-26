// Package output provides output formatting for annotations.
package output

import (
	"fmt"
	"io"
	"regexp"
	"sort"
	"strings"
	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Category classifies the comparison result for a single column.
type Category string

const (
	CatMatch            Category = "match"
	CatBothEmpty        Category = "both_empty"
	CatFuzzyFS          Category = "fuzzy_fs"
	CatSpliceVsSyn      Category = "splice_vs_syn"
	CatMafNonstandard   Category = "maf_nonstandard"
	CatSpliceNoProtein  Category = "splice_no_protein"
	CatPositionShift    Category = "position_shift"
	CatVepEmpty         Category = "vep_empty"
	CatMafEmpty         Category = "maf_empty"
	CatUpstreamReclass  Category = "upstream_reclassified"
	CatNoCDS            Category = "no_cds_data"
	CatDupVsIns         Category = "dup_vs_ins"
	CatMismatch         Category = "mismatch"
)

// isShownByDefault returns whether rows with this category are shown without --all.
func isShownByDefault(col string, cat Category) bool {
	switch cat {
	case CatMatch, CatBothEmpty:
		return false
	case CatMafNonstandard, CatSpliceNoProtein, CatNoCDS, CatDupVsIns:
		return false
	}
	return true
}

// CompareWriter writes tab-delimited comparison output between MAF annotations
// and VEP predictions, with category-based classification.
type CompareWriter struct {
	w       io.Writer
	columns map[string]bool             // enabled columns
	counts  map[string]map[Category]int // column → category → count
	total   int
	showAll bool
}

// NewCompareWriter creates a new comparison output writer.
// columns is a set of column names to compare (e.g. "consequence", "hgvsp", "hgvsc").
func NewCompareWriter(w io.Writer, columns map[string]bool, showAll bool) *CompareWriter {
	counts := make(map[string]map[Category]int)
	for col := range columns {
		counts[col] = make(map[Category]int)
	}
	return &CompareWriter{
		w:       w,
		columns: columns,
		counts:  counts,
		showAll: showAll,
	}
}

// WriteHeader writes the comparison output header.
func (c *CompareWriter) WriteHeader() error {
	parts := []string{"Variant", "Gene"}
	if c.columns["consequence"] {
		parts = append(parts, "MAF_Consequence", "VEP_Consequence")
	}
	if c.columns["hgvsp"] {
		parts = append(parts, "MAF_HGVSp", "VEP_HGVSp")
	}
	if c.columns["hgvsc"] {
		parts = append(parts, "MAF_HGVSc", "VEP_HGVSc")
	}
	if c.columns["consequence"] {
		parts = append(parts, "Conseq")
	}
	if c.columns["hgvsp"] {
		parts = append(parts, "HGVSp")
	}
	if c.columns["hgvsc"] {
		parts = append(parts, "HGVSc")
	}
	_, err := fmt.Fprintln(c.w, strings.Join(parts, "\t"))
	return err
}

// WriteComparison writes a comparison between MAF annotation and VEP prediction.
func (c *CompareWriter) WriteComparison(variant *vcf.Variant, mafAnn *maf.MAFAnnotation, vepAnns []*annotate.Annotation) error {
	c.total++

	bestAnn := SelectBestAnnotation(mafAnn, vepAnns)

	variantStr := fmt.Sprintf("%s:%d %s>%s", variant.Chrom, variant.Pos, variant.Ref, variant.Alt)

	mafConseq := mafAnn.Consequence
	mafHGVSp := mafAnn.HGVSpShort
	mafHGVSc := mafAnn.HGVSc

	var vepConseq, vepHGVSp, vepHGVSc string
	if bestAnn != nil {
		vepConseq = bestAnn.Consequence
		vepHGVSp = hgvspToShort(bestAnn.HGVSp)
		vepHGVSc = bestAnn.HGVSc
	}

	// Categorize each enabled column
	categories := make(map[string]Category)
	if c.columns["consequence"] {
		cat := categorizeConsequence(mafConseq, vepConseq)
		categories["consequence"] = cat
		c.counts["consequence"][cat]++
	}
	if c.columns["hgvsp"] {
		cat := categorizeHGVSp(mafHGVSp, vepHGVSp, vepConseq, mafHGVSc, vepHGVSc)
		categories["hgvsp"] = cat
		c.counts["hgvsp"][cat]++
	}
	if c.columns["hgvsc"] {
		cat := categorizeHGVSc(mafHGVSc, vepHGVSc)
		categories["hgvsc"] = cat
		c.counts["hgvsc"][cat]++
	}

	// Decide whether to show row
	showRow := c.showAll
	if !showRow {
		for col, cat := range categories {
			if isShownByDefault(col, cat) {
				showRow = true
				break
			}
		}
	}

	if showRow {
		parts := []string{variantStr, mafAnn.HugoSymbol}
		if c.columns["consequence"] {
			parts = append(parts, mafConseq, vepConseq)
		}
		if c.columns["hgvsp"] {
			parts = append(parts, mafHGVSp, vepHGVSp)
		}
		if c.columns["hgvsc"] {
			parts = append(parts, mafHGVSc, vepHGVSc)
		}
		if c.columns["consequence"] {
			parts = append(parts, string(categories["consequence"]))
		}
		if c.columns["hgvsp"] {
			parts = append(parts, string(categories["hgvsp"]))
		}
		if c.columns["hgvsc"] {
			parts = append(parts, string(categories["hgvsc"]))
		}
		_, err := fmt.Fprintln(c.w, strings.Join(parts, "\t"))
		return err
	}

	return nil
}

// Flush is a no-op (kept for interface compatibility).
func (c *CompareWriter) Flush() error {
	return nil
}

// Total returns the total number of variants compared.
func (c *CompareWriter) Total() int {
	return c.total
}

// Counts returns the category counts for all columns.
func (c *CompareWriter) Counts() map[string]map[Category]int {
	return c.counts
}

// WriteSummary writes category counts per column to the given writer.
func (c *CompareWriter) WriteSummary(w io.Writer) {
	fmt.Fprintf(w, "\nComparison Summary (%d variants):\n", c.total)

	// Print columns in stable order
	colOrder := []string{"consequence", "hgvsp", "hgvsc"}
	colNames := map[string]string{
		"consequence": "Consequence",
		"hgvsp":       "HGVSp",
		"hgvsc":       "HGVSc",
	}

	for _, col := range colOrder {
		if !c.columns[col] {
			continue
		}
		cats := c.counts[col]
		fmt.Fprintf(w, "\n  %s:\n", colNames[col])

		// Sort categories by count descending
		type catCount struct {
			cat   Category
			count int
		}
		var sorted []catCount
		for cat, count := range cats {
			sorted = append(sorted, catCount{cat, count})
		}
		sort.Slice(sorted, func(i, j int) bool {
			if sorted[i].count != sorted[j].count {
				return sorted[i].count > sorted[j].count
			}
			return sorted[i].cat < sorted[j].cat
		})

		for _, cc := range sorted {
			fmt.Fprintf(w, "    %-20s%d\n", cc.cat, cc.count)
		}
	}
}

// categorizeConsequence classifies the consequence comparison.
func categorizeConsequence(mafConseq, vepConseq string) Category {
	normMAF := normalizeConsequence(mafConseq)
	normVEP := normalizeConsequence(vepConseq)

	if normMAF == normVEP {
		return CatMatch
	}

	// When MAF says upstream/downstream, our tool may find a more specific consequence.
	if normMAF == "downstream_gene_variant" || normMAF == "upstream_gene_variant" {
		return CatUpstreamReclass
	}

	// coding_sequence_variant means we lacked CDS data to determine the specific
	// coding consequence. Accept any coding consequence from MAF as a match.
	if normVEP == "coding_sequence_variant" && isCodingConsequence(mafConseq) {
		return CatNoCDS
	}

	// Primary term match covers sub-annotation differences.
	if primaryConsequence(normMAF) == primaryConsequence(normVEP) {
		return CatMatch
	}

	return CatMismatch
}

// spliceVsSynRe matches MAF splice notation like p.X125_splice.
var spliceVsSynRe = regexp.MustCompile(`^p\.[A-Z]\d+_splice$`)

// synNotationRe matches VEP synonymous notation like p.Xxx123= or p.X123=.
var synNotationRe = regexp.MustCompile(`^p\.\w+=`)

// categorizeHGVSp classifies the HGVSp comparison.
// Both mafHGVSp and vepHGVSp should be in single-letter amino acid format.
func categorizeHGVSp(mafHGVSp, vepHGVSp, vepConseq, mafHGVSc, vepHGVSc string) Category {
	if mafHGVSp == "" && vepHGVSp == "" {
		return CatBothEmpty
	}

	if mafHGVSp == vepHGVSp {
		return CatMatch
	}

	if isFrameshiftHGVSp(mafHGVSp) && isFrameshiftHGVSp(vepHGVSp) {
		return CatFuzzyFS
	}

	if spliceVsSynRe.MatchString(mafHGVSp) && synNotationRe.MatchString(vepHGVSp) {
		return CatSpliceVsSyn
	}

	if isNonStandardIntronicHGVSp(mafHGVSp) && vepHGVSp == "" {
		return CatMafNonstandard
	}

	if isSpliceConsequence(vepConseq) && vepHGVSp == "" {
		return CatSpliceNoProtein
	}

	if strings.Contains(vepConseq, "coding_sequence_variant") && vepHGVSp == "" {
		return CatNoCDS
	}

	if mafHGVSp == "" && vepHGVSp != "" {
		return CatMafEmpty
	}

	if vepHGVSp == "" && mafHGVSp != "" {
		return CatVepEmpty
	}

	// Position shift: both non-empty, same amino acid change but different positions.
	// Covers GENCODE version differences where transcript coordinates shifted.
	if mafHGVSp != "" && vepHGVSp != "" {
		if hgvspChangeType(mafHGVSp) == hgvspChangeType(vepHGVSp) {
			return CatPositionShift
		}
		// Also catch indel position shifts where HGVSc matches but HGVSp positions differ
		if hgvscValuesMatch(mafHGVSc, vepHGVSc) {
			return CatPositionShift
		}
	}

	return CatMismatch
}

// categorizeHGVSc classifies the HGVSc comparison.
func categorizeHGVSc(mafHGVSc, vepHGVSc string) Category {
	if mafHGVSc == "" && vepHGVSc == "" {
		return CatBothEmpty
	}
	if hgvscValuesMatch(mafHGVSc, vepHGVSc) {
		return CatMatch
	}
	if mafHGVSc == "" {
		return CatMafEmpty
	}
	if vepHGVSc == "" {
		return CatVepEmpty
	}
	// Position shift: same operation (base change or indel type) at different positions.
	if hgvscOperation(mafHGVSc) != "" && hgvscOperation(mafHGVSc) == hgvscOperation(vepHGVSc) {
		return CatPositionShift
	}

	// Dup vs ins: one side reports a duplication, the other an insertion.
	// This occurs when we lack reference sequence context for non-CDS positions.
	mafNorm := mafHGVSc
	if idx := strings.LastIndex(mafNorm, ":"); idx >= 0 {
		mafNorm = mafNorm[idx+1:]
	}
	if (strings.Contains(mafNorm, "dup") && strings.Contains(vepHGVSc, "ins")) ||
		(strings.Contains(mafNorm, "ins") && strings.Contains(vepHGVSc, "dup")) {
		return CatDupVsIns
	}

	return CatMismatch
}

// hgvscOpRe extracts the operation from an HGVSc value:
// SNVs: "A>G" from "c.1428A>G"
// Indels: "del", "dup", "delinsGGC", "insTTC" etc.
var hgvscOpRe = regexp.MustCompile(`([ACGT]>[ACGT]|(?:del|dup|ins)\w*)$`)

// hgvscOperation extracts the change operation from an HGVSc value,
// stripping transcript prefix and position information.
func hgvscOperation(hgvsc string) string {
	// Strip transcript prefix
	if idx := strings.LastIndex(hgvsc, ":"); idx >= 0 {
		hgvsc = hgvsc[idx+1:]
	}
	return hgvscOpRe.FindString(hgvsc)
}

// --- Shared helpers (moved from validation.go) ---

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

// AnnotationBetter returns true if ann is a better pick than current for comparison.
// Priority: canonical > protein-coding biotype > higher impact > has HGVSp.
func AnnotationBetter(ann, current *annotate.Annotation) bool {
	if ann.IsCanonical != current.IsCanonical {
		return ann.IsCanonical
	}
	annCoding := isProteinCodingBiotype(ann.Biotype)
	curCoding := isProteinCodingBiotype(current.Biotype)
	if annCoding != curCoding {
		return annCoding
	}
	annImpact := annotate.ImpactRank(ann.Impact)
	curImpact := annotate.ImpactRank(current.Impact)
	if annImpact != curImpact {
		return annImpact > curImpact
	}
	// Prefer annotations with HGVSp (e.g. protein_coding over protein_coding_LoF
	// which may lack CDS sequence data)
	if (ann.HGVSp != "") != (current.HGVSp != "") {
		return ann.HGVSp != ""
	}
	return false
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

// isNonStandardIntronicHGVSp checks for MAF-specific p.*###* notation
// used for intronic variants. This is not standard HGVS.
func isNonStandardIntronicHGVSp(hgvsp string) bool {
	return len(hgvsp) > 3 && strings.HasPrefix(hgvsp, "p.*") && strings.HasSuffix(hgvsp, "*")
}

// isFrameshiftHGVSp checks if an HGVSp value represents a frameshift.
func isFrameshiftHGVSp(hgvsp string) bool {
	return strings.Contains(hgvsp, "fs")
}

// isSpliceConsequence checks if a consequence string contains a splice site term.
func isSpliceConsequence(conseq string) bool {
	return strings.Contains(conseq, "splice_donor_variant") ||
		strings.Contains(conseq, "splice_acceptor_variant")
}

// hgvspToShort converts 3-letter HGVSp notation to single-letter.
// e.g., "p.Gly12Cys" → "p.G12C", "p.Ter130=" → "p.*130="
func hgvspToShort(hgvsp string) string {
	result := hgvsp
	for single, three := range annotate.AminoAcidSingleToThree {
		if single == '*' {
			result = strings.ReplaceAll(result, three, "*")
		} else {
			result = strings.ReplaceAll(result, three, string(single))
		}
	}
	return result
}

// hgvspChangeType extracts the amino acid change signature from a single-letter
// HGVSp by stripping position numbers.
// e.g. "p.Y1145C" → "p.YC", "p.A735=" → "p.A=", "p.I874del" → "p.Idel"
func hgvspChangeType(hgvsp string) string {
	var b strings.Builder
	for _, r := range hgvsp {
		if r < '0' || r > '9' {
			b.WriteRune(r)
		}
	}
	return b.String()
}

// hgvscValuesMatch compares MAF HGVSc with VEP HGVSc.
// MAF HGVSc is prefixed with transcript ID (e.g., "ENST00000361923.2:c.1428C>G"),
// so we strip the prefix before comparing.
func hgvscValuesMatch(mafHGVSc, vepHGVSc string) bool {
	if mafHGVSc == "" && vepHGVSc == "" {
		return true
	}
	mafNorm := mafHGVSc
	if idx := strings.LastIndex(mafNorm, ":"); idx >= 0 {
		mafNorm = mafNorm[idx+1:]
	}
	return mafNorm == vepHGVSc
}

// consequencesMatch checks if MAF and VEP consequences should be considered matching.
func consequencesMatch(mafConseq, vepConseq string) bool {
	return categorizeConsequence(mafConseq, vepConseq) == CatMatch
}

// hgvspValuesMatch compares MAF HGVSp (single-letter) with VEP HGVSp (3-letter)
// by normalizing both to single-letter format.
func hgvspValuesMatch(mafHGVSp, vepHGVSp string) bool {
	if mafHGVSp == "" && vepHGVSp == "" {
		return true
	}
	vepShort := hgvspToShort(vepHGVSp)
	if mafHGVSp == vepShort {
		return true
	}
	if isFrameshiftHGVSp(mafHGVSp) && isFrameshiftHGVSp(vepShort) {
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
		"protein_altering_variant": "inframe_variant",
		"inframe_deletion":         "inframe_variant",
		"inframe_insertion":        "inframe_variant",
		"mature_mirna_variant":     "non_coding_transcript_exon_variant",
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
			if hasHighImpact && (t == "5_prime_utr_variant" || t == "3_prime_utr_variant") {
				continue
			}
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
