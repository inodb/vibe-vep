package output

import (
	"regexp"
	"sort"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
)

// Category classifies the comparison result for a single column.
type Category string

const (
	CatMatch                 Category = "match"
	CatBothEmpty             Category = "both_empty"
	CatFuzzyFS               Category = "fuzzy_fs"
	CatSpliceVsSyn           Category = "splice_vs_syn"
	CatMafNonstandard        Category = "maf_nonstandard"
	CatSpliceNoProtein       Category = "splice_no_protein"
	CatPositionShift         Category = "position_shift"
	CatVepEmpty              Category = "vep_empty"
	CatMafEmpty              Category = "maf_empty"
	CatUpstreamReclass       Category = "upstream_reclassified"
	CatNoCDS                 Category = "no_cds_data"
	CatDupVsIns              Category = "dup_vs_ins"
	CatDelinsNorm            Category = "delins_normalized"
	CatSpliceVsPredicted     Category = "splice_vs_predicted"
	CatTranscriptModelChange Category = "transcript_model_change"
	CatGeneModelChange       Category = "gene_model_change"
	CatMismatch              Category = "mismatch"
)

// Categorizer applies semantic categorization to column value comparisons.
type Categorizer struct{}

// CategorizeRow categorizes all compared columns for a single variant,
// including cross-column corrections (e.g. HGVSc position_shift → reclassify consequence).
// columns are the display names, left/right map column names to values,
// leftNames/rightNames are the actual column names in each file (parallel to columns).
// Returns column display name → Category.
func (c *Categorizer) CategorizeRow(columns []string, left, right map[string]string, leftNames, rightNames []string) map[string]Category {
	categories := make(map[string]Category, len(columns))

	// First pass: categorize each column
	for i, col := range columns {
		lv := left[leftNames[i]]
		rv := right[rightNames[i]]
		categories[col] = c.categorizeColumn(col, lv, rv, left, right, leftNames, rightNames, columns)
	}

	// Second pass: cross-column corrections
	c.applyCrossColumnCorrections(categories, columns)

	return categories
}

// categorizeColumn applies semantic categorization to a single column comparison.
func (c *Categorizer) categorizeColumn(col, lv, rv string, left, right map[string]string, leftNames, rightNames, columns []string) Category {
	colLower := strings.ToLower(col)

	switch {
	case colLower == "consequence" || colLower == "variant_classification":
		return categorizeConsequence(lv, rv)

	case colLower == "hgvsp_short" || colLower == "hgvsp":
		// Need HGVSc and Consequence from both rows for context
		var leftHGVSc, rightHGVSc, rightConseq string
		for j, c2 := range columns {
			c2Lower := strings.ToLower(c2)
			if c2Lower == "hgvsc" {
				leftHGVSc = left[leftNames[j]]
				rightHGVSc = right[rightNames[j]]
			}
			if c2Lower == "consequence" || c2Lower == "variant_classification" {
				rightConseq = right[rightNames[j]]
			}
		}
		return categorizeHGVSp(lv, rv, rightConseq, leftHGVSc, rightHGVSc)

	case colLower == "hgvsc":
		return categorizeHGVSc(lv, rv)

	default:
		// Simple string equality for unknown columns
		if lv == rv {
			if lv == "" {
				return CatBothEmpty
			}
			return CatMatch
		}
		return CatMismatch
	}
}

// applyCrossColumnCorrections applies cross-column corrections to categories.
func (c *Categorizer) applyCrossColumnCorrections(categories map[string]Category, columns []string) {
	// Find HGVSc category
	var hgvscCat Category
	for _, col := range columns {
		if strings.ToLower(col) == "hgvsc" {
			hgvscCat = categories[col]
			break
		}
	}

	if hgvscCat == "" {
		return
	}

	// Determine which categories to reclassify based on HGVSc
	var reclassTo Category
	var affectedCols []string

	switch hgvscCat {
	case CatPositionShift:
		reclassTo = CatPositionShift
		affectedCols = consequenceAndHGVSpCols(columns)
	case CatDelinsNorm:
		reclassTo = CatDelinsNorm
		affectedCols = consequenceAndHGVSpCols(columns)
	case CatDupVsIns:
		reclassTo = CatDupVsIns
		affectedCols = hgvspCols(columns)
	case CatTranscriptModelChange:
		reclassTo = CatTranscriptModelChange
		affectedCols = consequenceAndHGVSpCols(columns)
	}

	if reclassTo == "" {
		return
	}

	for _, col := range affectedCols {
		if categories[col] == CatMismatch {
			categories[col] = reclassTo
		}
	}
}

// consequenceAndHGVSpCols returns column names matching consequence or hgvsp.
func consequenceAndHGVSpCols(columns []string) []string {
	var result []string
	for _, col := range columns {
		lower := strings.ToLower(col)
		if lower == "consequence" || lower == "variant_classification" ||
			lower == "hgvsp" || lower == "hgvsp_short" {
			result = append(result, col)
		}
	}
	return result
}

// hgvspCols returns column names matching hgvsp.
func hgvspCols(columns []string) []string {
	var result []string
	for _, col := range columns {
		lower := strings.ToLower(col)
		if lower == "hgvsp" || lower == "hgvsp_short" {
			result = append(result, col)
		}
	}
	return result
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

	// Transcript biotype changed between GENCODE versions: a transcript that was
	// protein_coding in the MAF's version may be non-coding in ours (or vice
	// versa). When one side has a coding consequence and the other has a
	// non-coding consequence, this is a transcript model difference.
	if isCodingConsequence(mafConseq) && isNonCodingConsequence(normVEP) {
		return CatTranscriptModelChange
	}
	if isNonCodingConsequence(normMAF) && isCodingConsequence(vepConseq) {
		return CatTranscriptModelChange
	}

	// Non-coding ↔ non-coding: exon boundary shifts in non-coding transcripts
	// between GENCODE versions (intron_variant ↔ non_coding_transcript_exon_variant).
	if isNonCodingConsequence(normMAF) && isNonCodingConsequence(normVEP) {
		return CatTranscriptModelChange
	}

	// Non-coding exon ↔ UTR: transcript was non-coding (non_coding_transcript_exon_variant)
	// in one version and coding (UTR) in another. This is a biotype change.
	// Note: intron ↔ UTR is handled below as an exon boundary shift.
	if normMAF == "non_coding_transcript_exon_variant" && isUTRConsequence(normVEP) {
		return CatTranscriptModelChange
	}
	if isUTRConsequence(normMAF) && normVEP == "non_coding_transcript_exon_variant" {
		return CatTranscriptModelChange
	}

	// Non-coding ↔ intergenic: gene model boundary change for non-coding
	// transcripts.
	if isNonCodingConsequence(normMAF) && normVEP == "intergenic_variant" {
		return CatGeneModelChange
	}
	if normMAF == "intergenic_variant" && isNonCodingConsequence(normVEP) {
		return CatGeneModelChange
	}

	// Gene model boundary change: one side has a coding consequence and the
	// other is intergenic. The gene may exist in one GENCODE version but not
	// the other.
	if isCodingConsequence(mafConseq) && normVEP == "intergenic_variant" {
		return CatGeneModelChange
	}
	if normMAF == "intergenic_variant" && isCodingConsequence(vepConseq) {
		return CatGeneModelChange
	}

	// UTR ↔ intron: exon boundary shifts between GENCODE versions can move a
	// position from UTR to intron or vice versa. Reclassify similarly to
	// upstream/downstream reclassification.
	if (isUTRConsequence(normMAF) || normMAF == "intron_variant") &&
		(isUTRConsequence(normVEP) || normVEP == "intron_variant") {
		return CatUpstreamReclass
	}

	// Primary term match covers sub-annotation differences.
	if primaryConsequence(normMAF) == primaryConsequence(normVEP) {
		return CatMatch
	}

	// Frameshift ↔ stop_gained/start_lost: when a frameshift creates an
	// immediate stop codon or occurs at the start codon boundary, VEP
	// convention may differ from MAF.
	mafPrimary := primaryConsequence(normMAF)
	vepPrimary := primaryConsequence(normVEP)
	if mafPrimary == "frameshift_variant" || vepPrimary == "frameshift_variant" {
		other := vepPrimary
		if vepPrimary == "frameshift_variant" {
			other = mafPrimary
		}
		if other == "stop_gained" || other == "start_lost" {
			return CatMatch
		}
	}

	// One annotation is more specific: one consequence contains the other's
	// primary term (e.g., inframe_deletion → stop_gained,inframe_deletion,
	// or frameshift_variant,start_lost → start_lost).
	if strings.Contains(normVEP, mafPrimary) || strings.Contains(normMAF, vepPrimary) {
		return CatMatch
	}

	// Splice site reclassification: when an indel at a splice boundary is
	// classified as frameshift+splice_region or splice_region+intron by one
	// annotator but as splice_acceptor/donor by another, these represent
	// the same variant effect. 3' normalization can shift an insertion from
	// the exonic side (frameshift) to the intronic side (splice site).
	if spliceReclassMatch(normMAF, normVEP) || spliceReclassMatch(normVEP, normMAF) {
		return CatMatch
	}

	// Inframe indel ↔ stop_gained: GENCODE version differences can change
	// whether a junction codon translates to a stop. Accept when the primary
	// consequence is inframe_variant on one side and stop_gained on the other.
	if (mafPrimary == "inframe_variant" && vepPrimary == "stop_gained") ||
		(mafPrimary == "stop_gained" && vepPrimary == "inframe_variant") {
		return CatMatch
	}

	// Inframe indel ↔ stop_lost: GENCODE version differences at the stop
	// codon boundary. Accept inframe_variant ↔ stop_lost.
	if (mafPrimary == "inframe_variant" && vepPrimary == "stop_lost") ||
		(mafPrimary == "stop_lost" && vepPrimary == "inframe_variant") {
		return CatMatch
	}

	// stop_lost ↔ stop_retained: CDS boundary differences can change whether
	// the stop codon is lost or retained. Both represent stop codon effects.
	if (mafPrimary == "stop_lost" && vepPrimary == "stop_retained_variant") ||
		(mafPrimary == "stop_retained_variant" && vepPrimary == "stop_lost") {
		return CatMatch
	}

	// Synonymous ↔ stop_retained: CDS boundary differences between GENCODE
	// versions can place a variant at the stop codon in one version but not
	// the other. Both are LOW impact silent changes.
	if (mafPrimary == "synonymous_variant" && vepPrimary == "stop_retained_variant") ||
		(mafPrimary == "stop_retained_variant" && vepPrimary == "synonymous_variant") {
		return CatMatch
	}

	// start_lost ↔ coding: CDS start boundary differences between GENCODE
	// versions can move a variant into or out of the start codon region.
	if mafPrimary == "start_lost" || vepPrimary == "start_lost" {
		other := vepPrimary
		if mafPrimary != "start_lost" {
			other = mafPrimary
		}
		switch other {
		case "synonymous_variant", "missense_variant",
			"inframe_variant", "inframe_insertion", "inframe_deletion":
			return CatMatch
		}
	}

	return CatMismatch
}

// spliceReclassMatch returns true if conseqA has a splice_region or frameshift
// annotation and conseqB is a splice_donor/acceptor reclassification of the
// same variant.
func spliceReclassMatch(conseqA, conseqB string) bool {
	if conseqB != "splice_acceptor_variant" && conseqB != "splice_donor_variant" {
		return false
	}
	return strings.Contains(conseqA, "splice_region_variant") ||
		strings.Contains(conseqA, "frameshift_variant") ||
		strings.Contains(conseqA, "inframe_variant") ||
		conseqA == "start_lost"
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

	// Frameshift ↔ stop_gained: a frameshift that immediately creates a stop
	// codon can be reported as either frameshift or stop_gained.
	if (isFrameshiftHGVSp(mafHGVSp) && isStopGainedHGVSp(vepHGVSp)) ||
		(isStopGainedHGVSp(mafHGVSp) && isFrameshiftHGVSp(vepHGVSp)) {
		return CatFuzzyFS
	}

	if spliceVsSynRe.MatchString(mafHGVSp) && synNotationRe.MatchString(vepHGVSp) {
		return CatSpliceVsSyn
	}

	if isNonStandardIntronicHGVSp(mafHGVSp) {
		if vepHGVSp == "" {
			return CatMafNonstandard
		}
		return CatTranscriptModelChange
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
	if mafHGVSp != "" && vepHGVSp != "" {
		if hgvspChangeType(mafHGVSp) == hgvspChangeType(vepHGVSp) {
			return CatPositionShift
		}
		if mafHGVSc != "" && vepHGVSc != "" && hgvscValuesMatch(mafHGVSc, vepHGVSc) {
			return CatPositionShift
		}
	}

	// Splice vs predicted: MAF has non-standard p.X###_splice notation
	// while we provide an actual protein-level prediction.
	if spliceVsSynRe.MatchString(mafHGVSp) && vepHGVSp != "" {
		return CatSpliceVsPredicted
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

	// Strip transcript prefix for further comparisons.
	mafNorm := mafHGVSc
	if idx := strings.LastIndex(mafNorm, ":"); idx >= 0 {
		mafNorm = mafNorm[idx+1:]
	}

	// Transcript model change: one side uses n. notation (non-coding) and the
	// other uses c. notation (coding).
	mafIsNonCodingHGVSc := strings.HasPrefix(mafNorm, "n.")
	vepIsNonCodingHGVSc := strings.HasPrefix(vepHGVSc, "n.")
	if mafIsNonCodingHGVSc != vepIsNonCodingHGVSc {
		return CatTranscriptModelChange
	}

	// Both non-coding (n.) but different coordinates.
	if mafIsNonCodingHGVSc && vepIsNonCodingHGVSc {
		return CatTranscriptModelChange
	}

	// Position shift: same operation at different positions.
	if hgvscOperation(mafNorm) != "" && hgvscOperation(mafNorm) == hgvscOperation(vepHGVSc) {
		return CatPositionShift
	}

	// Dup vs ins.
	if (strings.Contains(mafNorm, "dup") && strings.Contains(vepHGVSc, "ins")) ||
		(strings.Contains(mafNorm, "ins") && strings.Contains(vepHGVSc, "dup")) {
		return CatDupVsIns
	}

	// Delins normalization.
	if strings.Contains(mafNorm, "delins") || strings.Contains(vepHGVSc, "delins") {
		return CatDelinsNorm
	}

	// Insertion cyclic rotation.
	if isCyclicRotationIns(mafNorm, vepHGVSc) {
		return CatPositionShift
	}

	return CatMismatch
}

// hgvscInsSeqRe extracts the inserted sequence from an HGVSc insertion.
var hgvscInsSeqRe = regexp.MustCompile(`ins([ACGT]+)$`)

// isCyclicRotationIns checks if two HGVSc insertions have inserted sequences
// that are cyclic rotations of each other.
func isCyclicRotationIns(a, b string) bool {
	mA := hgvscInsSeqRe.FindStringSubmatch(a)
	mB := hgvscInsSeqRe.FindStringSubmatch(b)
	if len(mA) < 2 || len(mB) < 2 {
		return false
	}
	seqA, seqB := mA[1], mB[1]
	if len(seqA) != len(seqB) || len(seqA) == 0 {
		return false
	}
	return strings.Contains(seqA+seqA, seqB)
}

// hgvscOpRe extracts the operation from an HGVSc value.
var hgvscOpRe = regexp.MustCompile(`([ACGT]>[ACGT]|(?:del|dup|ins)\w*)$`)

// hgvscOperation extracts the change operation from an HGVSc value.
func hgvscOperation(hgvsc string) string {
	if idx := strings.LastIndex(hgvsc, ":"); idx >= 0 {
		hgvsc = hgvsc[idx+1:]
	}
	return hgvscOpRe.FindString(hgvsc)
}

// --- Shared helpers ---

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

// isNonStandardIntronicHGVSp checks for MAF-specific p.*###* notation.
func isNonStandardIntronicHGVSp(hgvsp string) bool {
	return len(hgvsp) > 3 && strings.HasPrefix(hgvsp, "p.*") && strings.HasSuffix(hgvsp, "*")
}

// isFrameshiftHGVSp checks if an HGVSp value represents a frameshift.
func isFrameshiftHGVSp(hgvsp string) bool {
	return strings.Contains(hgvsp, "fs")
}

// isStopGainedHGVSp detects stop_gained notation like p.K453*.
func isStopGainedHGVSp(hgvsp string) bool {
	return strings.HasSuffix(hgvsp, "*") && !strings.Contains(hgvsp, "fs")
}

// isSpliceConsequence checks if a consequence string contains a splice site term.
func isSpliceConsequence(conseq string) bool {
	return strings.Contains(conseq, "splice_donor_variant") ||
		strings.Contains(conseq, "splice_acceptor_variant")
}

// isUTRConsequence checks if a consequence is a UTR variant.
func isUTRConsequence(conseq string) bool {
	return strings.Contains(conseq, "5_prime_utr_variant") ||
		strings.Contains(conseq, "3_prime_utr_variant")
}

// isNonCodingConsequence returns true if the consequence indicates a non-coding transcript region.
func isNonCodingConsequence(conseq string) bool {
	for _, term := range strings.Split(conseq, ",") {
		switch strings.TrimSpace(term) {
		case "non_coding_transcript_exon_variant", "intron_variant":
			return true
		}
	}
	return false
}

// hgvspChangeType extracts the amino acid change signature by stripping position numbers.
func hgvspChangeType(hgvsp string) string {
	var b strings.Builder
	for _, r := range hgvsp {
		if r < '0' || r > '9' {
			b.WriteRune(r)
		}
	}
	return b.String()
}

// hgvscValuesMatch compares MAF HGVSc with VEP HGVSc, stripping transcript prefix.
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
func normalizeConsequence(conseq string) string {
	conseq = strings.ToLower(strings.TrimSpace(conseq))

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

	var terms []string
	for _, term := range strings.Split(conseq, ",") {
		term = strings.TrimSpace(term)
		if mapped, ok := mappings[term]; ok {
			term = mapped
		}
		switch term {
		case "splice_donor_region_variant", "splice_donor_5th_base_variant":
			term = "splice_region_variant"
		}
		switch term {
		case "non_coding_transcript_variant", "nmd_transcript_variant",
			"splice_polypyrimidine_tract_variant":
			continue
		}
		terms = append(terms, term)
	}

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
			if hasHighImpact && (t == "start_retained_variant" || t == "stop_retained_variant" || t == "coding_sequence_variant") {
				continue
			}
			filtered = append(filtered, t)
		}
		terms = filtered
	}

	sort.Strings(terms)
	return strings.Join(terms, ",")
}

// consequencesMatch checks if MAF and VEP consequences should be considered matching.
func consequencesMatch(mafConseq, vepConseq string) bool {
	return categorizeConsequence(mafConseq, vepConseq) == CatMatch
}

// hgvspValuesMatch compares MAF HGVSp (single-letter) with VEP HGVSp (3-letter).
func hgvspValuesMatch(mafHGVSp, vepHGVSp string) bool {
	if mafHGVSp == "" && vepHGVSp == "" {
		return true
	}
	vepShort := HGVSpToShort(vepHGVSp)
	if mafHGVSp == vepShort {
		return true
	}
	if isFrameshiftHGVSp(mafHGVSp) && isFrameshiftHGVSp(vepShort) {
		return true
	}
	return false
}
