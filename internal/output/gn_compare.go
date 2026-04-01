package output

import (
	"encoding/json"
	"fmt"
	"math"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
)

// GNAnnotation represents a genome-nexus annotation response.
type GNAnnotation struct {
	Variant               string                    `json:"variant"`
	OriginalVariantQuery  string                    `json:"originalVariantQuery"`
	HGVSg                 string                    `json:"hgvsg"`
	ID                    string                    `json:"id"`
	AssemblyName          string                    `json:"assembly_name"`
	SeqRegionName         string                    `json:"seq_region_name"`
	Start                 int64                     `json:"start"`
	End                   int64                     `json:"end"`
	AlleleString          string                    `json:"allele_string"`
	Strand                int                       `json:"strand"`
	MostSevereConsequence string                    `json:"most_severe_consequence"`
	TranscriptConsequences []GNTranscriptConsequence `json:"transcript_consequences"`
	SuccessfullyAnnotated bool                      `json:"successfully_annotated"`
	AnnotationSummary     *GNAnnotationSummary      `json:"annotation_summary,omitempty"`
	ClinVar               *GNClinVar                `json:"clinvar,omitempty"`
	ColocatedVariants     []GNColocatedVariant       `json:"colocatedVariants,omitempty"`
	Hotspots              *GNHotspots               `json:"hotspots,omitempty"`
	SignalAnnotation      *GNSignalAnnotation        `json:"signalAnnotation,omitempty"`
	MyVariantInfo         *GNMyVariantInfoAnnotation `json:"my_variant_info,omitempty"`
}

// GNMyVariantInfoAnnotation wraps the myvariant.info annotation in genome-nexus format.
type GNMyVariantInfoAnnotation struct {
	Annotation *GNMyVariantInfo `json:"annotation,omitempty"`
}

// GNMyVariantInfo holds the transformed myvariant.info data.
type GNMyVariantInfo struct {
	Dbsnp        *GNDbsnp               `json:"dbsnp,omitempty"`
	GnomadExome  *GNGnomad              `json:"gnomadExome,omitempty"`
	GnomadGenome *GNGnomad              `json:"gnomadGenome,omitempty"`
	Vcf          *GNVcf                 `json:"vcf,omitempty"`
	Variant      string                 `json:"variant,omitempty"`
	Query        string                 `json:"query,omitempty"`
	Hgvs         string                 `json:"hgvs,omitempty"`
}

// GNDbsnp holds dbSNP rsid.
type GNDbsnp struct {
	Rsid string `json:"rsid,omitempty"`
}

// GNGnomad holds gnomAD frequency data.
type GNGnomad struct {
	AlleleCount     map[string]interface{} `json:"alleleCount,omitempty"`
	AlleleFrequency map[string]interface{} `json:"alleleFrequency,omitempty"`
	AlleleNumber    map[string]interface{} `json:"alleleNumber,omitempty"`
	Homozygotes     map[string]interface{} `json:"homozygotes,omitempty"`
}

// GNVcf holds VCF-level variant info from myvariant.info.
type GNVcf struct {
	Ref      string `json:"ref,omitempty"`
	Alt      string `json:"alt,omitempty"`
	Position string `json:"position,omitempty"`
}

// GNSignalAnnotation represents SIGNAL annotation in the GN response.
type GNSignalAnnotation struct {
	License    string            `json:"license"`
	Annotation []GNSignalMutation `json:"annotation"`
}

// GNSignalMutation represents a SIGNAL mutation entry.
type GNSignalMutation struct {
	HugoGeneSymbol  string              `json:"hugoGeneSymbol,omitempty"`
	Chromosome      string              `json:"chromosome,omitempty"`
	StartPosition   int64               `json:"startPosition,omitempty"`
	EndPosition     int64               `json:"endPosition,omitempty"`
	ReferenceAllele string              `json:"referenceAllele,omitempty"`
	VariantAllele   string              `json:"variantAllele,omitempty"`
	MutationStatus  string              `json:"mutationStatus,omitempty"`
	CountsByTumorType []GNCountByTumorType `json:"countsByTumorType,omitempty"`
}

// GNCountByTumorType represents a count by tumor type entry.
type GNCountByTumorType struct {
	TumorType      string `json:"tumorType"`
	TumorTypeCount int    `json:"tumorTypeCount"`
	VariantCount   int    `json:"variantCount"`
}

// GNClinVar represents ClinVar annotation in the GN response.
type GNClinVar struct {
	Annotation *GNClinVarAnnotation `json:"annotation,omitempty"`
}

// GNClinVarAnnotation contains ClinVar details.
type GNClinVarAnnotation struct {
	Chromosome                      string `json:"chromosome,omitempty"`
	StartPosition                   int64  `json:"startPosition,omitempty"`
	EndPosition                     int64  `json:"endPosition,omitempty"`
	ReferenceAllele                 string `json:"referenceAllele,omitempty"`
	AlternateAllele                 string `json:"alternateAllele,omitempty"`
	ClinicalSignificance            string `json:"clinicalSignificance,omitempty"`
	ConflictingClinicalSignificance string `json:"conflictingClinicalSignificance,omitempty"`
	ReviewStatus                    string `json:"reviewStatus,omitempty"`
	DiseaseName                     string `json:"diseaseName,omitempty"`
}

// GNColocatedVariant represents a dbSNP/COSMIC ID at the same position.
type GNColocatedVariant struct {
	DbSnpID string `json:"dbSnpId,omitempty"`
}

// GNHotspots represents cancer hotspot annotation.
type GNHotspots struct {
	Annotation []GNHotspotEntry `json:"annotation,omitempty"`
}

// GNHotspotEntry represents a single hotspot hit.
type GNHotspotEntry struct {
	Type   string  `json:"type,omitempty"`
	QValue float64 `json:"qValue,omitempty"`
}

// GNAnnotationSummary is the enriched annotation summary returned when
// ?fields=annotation_summary is requested.
type GNAnnotationSummary struct {
	Variant                        string                              `json:"variant"`
	GenomicLocation                GNGenomicLocation                   `json:"genomicLocation"`
	StrandSign                     string                              `json:"strandSign"`
	VariantType                    string                              `json:"variantType"`
	AssemblyName                   string                              `json:"assemblyName"`
	CanonicalTranscriptID          string                              `json:"canonicalTranscriptId"`
	TranscriptConsequenceSummary   *GNTranscriptConsequenceSummary     `json:"transcriptConsequenceSummary"`
	TranscriptConsequenceSummaries []GNTranscriptConsequenceSummary    `json:"transcriptConsequenceSummaries"`
	TranscriptConsequences         []GNTranscriptConsequenceSummary    `json:"transcriptConsequences"`
}

// GNGenomicLocation represents a genomic location in the annotation summary.
type GNGenomicLocation struct {
	Chromosome      string `json:"chromosome"`
	Start           int64  `json:"start"`
	End             int64  `json:"end"`
	ReferenceAllele string `json:"referenceAllele"`
	VariantAllele   string `json:"variantAllele"`
}

// GNTranscriptConsequenceSummary is the enriched per-transcript summary.
type GNTranscriptConsequenceSummary struct {
	TranscriptID          string            `json:"transcriptId"`
	CodonChange           string            `json:"codonChange,omitempty"`
	AminoAcids            string            `json:"aminoAcids,omitempty"`
	AminoAcidRef          string            `json:"aminoAcidRef,omitempty"`
	AminoAcidAlt          string            `json:"aminoAcidAlt,omitempty"`
	EntrezGeneID          string            `json:"entrezGeneId,omitempty"`
	HugoGeneSymbol        string            `json:"hugoGeneSymbol,omitempty"`
	HGVSpShort            string            `json:"hgvspShort,omitempty"`
	HGVSp                 string            `json:"hgvsp,omitempty"`
	HGVSc                 string            `json:"hgvsc,omitempty"`
	ConsequenceTerms      string            `json:"consequenceTerms,omitempty"`
	VariantClassification string            `json:"variantClassification,omitempty"`
	Exon                  string            `json:"exon,omitempty"`
	ProteinPosition       *GNIntegerRange   `json:"proteinPosition,omitempty"`
	PolyphenScore         *float64          `json:"polyphenScore,omitempty"`
	PolyphenPrediction    string            `json:"polyphenPrediction,omitempty"`
	SIFTScore             *float64          `json:"siftScore,omitempty"`
	SIFTPrediction        string            `json:"siftPrediction,omitempty"`
}

// GNIntegerRange represents a start/end integer range.
type GNIntegerRange struct {
	Start int64 `json:"start"`
	End   int64 `json:"end"`
}

// GNTranscriptConsequence represents a transcript consequence from genome-nexus.
type GNTranscriptConsequence struct {
	TranscriptID       string   `json:"transcript_id"`
	GeneSymbol         string   `json:"gene_symbol"`
	GeneID             string   `json:"gene_id"`
	HGNCId             string   `json:"hgnc_id,omitempty"`
	ProteinID          string   `json:"protein_id,omitempty"`
	ConsequenceTerms   []string `json:"consequence_terms"`
	Impact             string   `json:"impact"`
	VariantAllele      string   `json:"variant_allele"`
	AminoAcids         string   `json:"amino_acids"`
	Codons             string   `json:"codons"`
	ProteinStart       int64    `json:"protein_start"`
	ProteinEnd         int64    `json:"protein_end"`
	CDSStart           int64    `json:"cds_start"`
	CDSEnd             int64    `json:"cds_end"`
	CDNAStart          int64    `json:"cdna_start"`
	CDNAEnd            int64    `json:"cdna_end"`
	HGVSp              string   `json:"hgvsp"`
	HGVSc              string   `json:"hgvsc"`
	Exon               string   `json:"exon"`
	Intron             string   `json:"intron"`
	Biotype            string   `json:"biotype"`
	Canonical          string   `json:"canonical"`
	SIFTScore          *float64 `json:"sift_score"`
	SIFTPrediction     string   `json:"sift_prediction"`
	PolyPhenScore      *float64 `json:"polyphen_score"`
	PolyPhenPrediction string   `json:"polyphen_prediction"`
}

// ParseGNAnnotation parses a JSON line into a GNAnnotation.
func ParseGNAnnotation(line []byte) (*GNAnnotation, error) {
	var gn GNAnnotation
	if err := json.Unmarshal(line, &gn); err != nil {
		return nil, fmt.Errorf("parse genome-nexus annotation: %w", err)
	}
	return &gn, nil
}

// GNFieldComparison holds the result of comparing a single field.
type GNFieldComparison struct {
	Field    string   `json:"field"`
	GNValue  string   `json:"gn_value"`
	VEPValue string   `json:"vep_value"`
	Match    bool     `json:"match"`
	Category Category `json:"category"`
}

// GNVariantComparison holds the comparison result for one variant.
type GNVariantComparison struct {
	VariantID    string              `json:"variant_id"`
	GNInput      string              `json:"gn_input"`
	TranscriptID string              `json:"transcript_id"`
	Fields       []GNFieldComparison `json:"fields"`
	AllMatch     bool                `json:"all_match"`
	Error        string              `json:"error,omitempty"`
}

// GNComparisonReport holds the summary of all variant comparisons.
type GNComparisonReport struct {
	TotalVariants    int                       `json:"total_variants"`
	FullMatches      int                       `json:"full_matches"`
	PartialMatches   int                       `json:"partial_matches"`
	Errors           int                       `json:"errors"`
	FieldStats       map[string]*GNFieldStats  `json:"field_stats"`
	CategoryCounts   map[Category]int          `json:"category_counts"`
	Comparisons      []GNVariantComparison     `json:"comparisons,omitempty"`
}

// GNFieldStats tracks match/mismatch counts per field.
type GNFieldStats struct {
	Total      int            `json:"total"`
	Matches    int            `json:"matches"`
	Mismatches int            `json:"mismatches"`
	BothEmpty  int            `json:"both_empty"`
	Categories map[Category]int `json:"categories"`
}

// NewGNComparisonReport creates a new comparison report.
func NewGNComparisonReport() *GNComparisonReport {
	return &GNComparisonReport{
		FieldStats:     make(map[string]*GNFieldStats),
		CategoryCounts: make(map[Category]int),
	}
}

// AddComparison adds a variant comparison to the report.
func (r *GNComparisonReport) AddComparison(cmp GNVariantComparison) {
	r.TotalVariants++
	if cmp.Error != "" {
		r.Errors++
		r.Comparisons = append(r.Comparisons, cmp)
		return
	}

	if cmp.AllMatch {
		r.FullMatches++
	} else {
		r.PartialMatches++
	}

	for _, f := range cmp.Fields {
		stats, ok := r.FieldStats[f.Field]
		if !ok {
			stats = &GNFieldStats{Categories: make(map[Category]int)}
			r.FieldStats[f.Field] = stats
		}
		stats.Total++
		if f.Match {
			stats.Matches++
		} else if f.Category == CatBothEmpty {
			stats.BothEmpty++
		} else {
			stats.Mismatches++
		}
		stats.Categories[f.Category]++
		r.CategoryCounts[f.Category]++
	}

	r.Comparisons = append(r.Comparisons, cmp)
}

// MatchRate returns the match rate as a percentage for a field.
func (s *GNFieldStats) MatchRate() float64 {
	applicable := s.Total - s.BothEmpty
	if applicable == 0 {
		return 100.0
	}
	return float64(s.Matches) / float64(applicable) * 100.0
}

// CompareGNFields is the list of fields to compare.
var CompareGNFields = []string{
	"gene_symbol",
	"consequence_terms",
	"amino_acids",
	"protein_start",
	"hgvsc",
	"hgvsp",
	"sift_score",
	"sift_prediction",
	"polyphen_score",
	"polyphen_prediction",
}

// CompareGNToVEP compares a genome-nexus annotation against a vibe-vep annotation
// on the canonical/best transcript.
func CompareGNToVEP(gn *GNAnnotation, vepAnns []*annotate.Annotation) GNVariantComparison {
	cmp := GNVariantComparison{
		VariantID: gn.Variant,
		GNInput:   gn.OriginalVariantQuery,
	}

	if len(gn.TranscriptConsequences) == 0 {
		cmp.Error = "no transcript consequences in genome-nexus response"
		return cmp
	}
	if len(vepAnns) == 0 {
		cmp.Error = "no annotations from vibe-vep"
		return cmp
	}

	// Pick canonical transcript from GN response.
	gnTC := pickGNCanonical(gn.TranscriptConsequences)

	// Find matching transcript in vibe-vep annotations.
	vepAnn := matchVEPTranscript(gnTC.TranscriptID, vepAnns)
	if vepAnn == nil {
		// Fall back to best annotation.
		vepAnn = PickBestAnnotation(vepAnns)
	}

	cmp.TranscriptID = gnTC.TranscriptID

	cmp.Fields = compareTranscriptFields(gnTC, vepAnn)

	cmp.AllMatch = true
	for _, f := range cmp.Fields {
		if !f.Match && f.Category != CatBothEmpty {
			cmp.AllMatch = false
			break
		}
	}

	return cmp
}

// pickGNCanonical selects the canonical transcript from GN consequences.
// Prefers canonical=="1", then protein_coding biotype, then first.
func pickGNCanonical(tcs []GNTranscriptConsequence) *GNTranscriptConsequence {
	var best *GNTranscriptConsequence
	for i := range tcs {
		tc := &tcs[i]
		if tc.Canonical == "1" {
			if best == nil || best.Canonical != "1" {
				best = tc
				continue
			}
			// Both canonical — prefer protein_coding.
			if isProteinCodingBiotype(tc.Biotype) && !isProteinCodingBiotype(best.Biotype) {
				best = tc
			}
			continue
		}
		if best == nil {
			best = tc
			continue
		}
		if best.Canonical != "1" && isProteinCodingBiotype(tc.Biotype) && !isProteinCodingBiotype(best.Biotype) {
			best = tc
		}
	}
	return best
}

// matchVEPTranscript finds a vibe-vep annotation matching the given transcript ID
// (ignoring version suffix).
func matchVEPTranscript(gnTranscriptID string, anns []*annotate.Annotation) *annotate.Annotation {
	gnBase := transcriptBaseID(gnTranscriptID)
	for _, ann := range anns {
		if transcriptBaseID(ann.TranscriptID) == gnBase {
			return ann
		}
	}
	return nil
}

// compareTranscriptFields compares individual fields between GN and vibe-vep.
func compareTranscriptFields(gn *GNTranscriptConsequence, vep *annotate.Annotation) []GNFieldComparison {
	fields := make([]GNFieldComparison, 0, len(CompareGNFields))

	for _, field := range CompareGNFields {
		gnVal, vepVal := extractFieldValues(field, gn, vep)
		fc := GNFieldComparison{
			Field:    field,
			GNValue:  gnVal,
			VEPValue: vepVal,
		}
		fc.Match, fc.Category = classifyFieldMatch(field, gnVal, vepVal)
		fields = append(fields, fc)
	}
	return fields
}

// extractFieldValues extracts the string values for a field from both sources.
func extractFieldValues(field string, gn *GNTranscriptConsequence, vep *annotate.Annotation) (string, string) {
	switch field {
	case "gene_symbol":
		return gn.GeneSymbol, vep.GeneName
	case "consequence_terms":
		return strings.Join(gn.ConsequenceTerms, ","), vep.Consequence
	case "amino_acids":
		return gn.AminoAcids, formatAminoAcidsVEP(vep.AminoAcidChange)
	case "protein_start":
		gnVal := ""
		if gn.ProteinStart > 0 {
			gnVal = fmt.Sprintf("%d", gn.ProteinStart)
		}
		vepVal := ""
		if vep.ProteinPosition > 0 {
			vepVal = fmt.Sprintf("%d", vep.ProteinPosition)
		}
		return gnVal, vepVal
	case "hgvsc":
		// GN includes transcript prefix (ENST...:c.xxx), vibe-vep stores just c.xxx
		gnVal := gn.HGVSc
		if idx := strings.LastIndex(gnVal, ":"); idx >= 0 {
			gnVal = gnVal[idx+1:]
		}
		return gnVal, vep.HGVSc
	case "hgvsp":
		// GN includes protein prefix (ENSP...:p.xxx), normalize both to short form.
		gnVal := gn.HGVSp
		if idx := strings.LastIndex(gnVal, ":"); idx >= 0 {
			gnVal = gnVal[idx+1:]
		}
		gnVal = HGVSpToShort(gnVal)
		vepVal := HGVSpToShort(vep.HGVSp)
		return gnVal, vepVal
	case "sift_score":
		gnVal := formatOptionalFloat(gn.SIFTScore)
		vepVal := vep.GetExtraKey("sift.score")
		return gnVal, vepVal
	case "sift_prediction":
		return gn.SIFTPrediction, vep.GetExtraKey("sift.prediction")
	case "polyphen_score":
		gnVal := formatOptionalFloat(gn.PolyPhenScore)
		vepVal := vep.GetExtraKey("polyphen.score")
		return gnVal, vepVal
	case "polyphen_prediction":
		return gn.PolyPhenPrediction, vep.GetExtraKey("polyphen.prediction")
	default:
		return "", ""
	}
}

// classifyFieldMatch determines if two field values match and categorizes mismatches.
func classifyFieldMatch(field, gnVal, vepVal string) (bool, Category) {
	if gnVal == "" && vepVal == "" {
		return true, CatBothEmpty
	}
	if gnVal == vepVal {
		return true, CatMatch
	}

	switch field {
	case "consequence_terms":
		cat := categorizeConsequence(gnVal, vepVal)
		return cat == CatMatch, cat
	case "hgvsc":
		cat := categorizeHGVSc(gnVal, vepVal)
		return cat == CatMatch, cat
	case "hgvsp":
		cat := categorizeHGVSp(gnVal, vepVal, vepVal, "", "")
		return cat == CatMatch || cat == CatFuzzyFS, cat
	case "sift_score", "polyphen_score":
		return floatScoresMatch(gnVal, vepVal), CatMismatch
	case "sift_prediction", "polyphen_prediction":
		// Exact match only.
		return false, CatMismatch
	case "amino_acids":
		// Normalize: VEP uses "/" separator, check both orderings.
		return false, CatMismatch
	default:
		return false, CatMismatch
	}
}

// floatScoresMatch compares two float score strings with tolerance for rounding.
func floatScoresMatch(a, b string) bool {
	if a == "" || b == "" {
		return a == b
	}
	var fa, fb float64
	if _, err := fmt.Sscanf(a, "%f", &fa); err != nil {
		return false
	}
	if _, err := fmt.Sscanf(b, "%f", &fb); err != nil {
		return false
	}
	return math.Abs(fa-fb) < 0.005
}

// formatOptionalFloat formats a *float64 as a string, empty if nil.
func formatOptionalFloat(f *float64) string {
	if f == nil {
		return ""
	}
	return fmt.Sprintf("%g", *f)
}

// FormatGNReport formats the comparison report as markdown.
func FormatGNReport(report *GNComparisonReport) string {
	var b strings.Builder

	b.WriteString("# Genome Nexus vs vibe-vep Comparison Report\n\n")
	b.WriteString(fmt.Sprintf("**Total variants:** %d\n", report.TotalVariants))
	b.WriteString(fmt.Sprintf("**Full matches:** %d (%.1f%%)\n",
		report.FullMatches, pct(report.FullMatches, report.TotalVariants)))
	b.WriteString(fmt.Sprintf("**Partial matches:** %d\n", report.PartialMatches))
	b.WriteString(fmt.Sprintf("**Errors:** %d\n\n", report.Errors))

	// Field-level stats table.
	b.WriteString("## Field Match Rates\n\n")
	b.WriteString("| Field | Total | Matches | Mismatches | Both Empty | Match Rate |\n")
	b.WriteString("|-------|------:|--------:|-----------:|-----------:|-----------:|\n")
	for _, field := range CompareGNFields {
		stats, ok := report.FieldStats[field]
		if !ok {
			continue
		}
		b.WriteString(fmt.Sprintf("| %s | %d | %d | %d | %d | %.1f%% |\n",
			field, stats.Total, stats.Matches, stats.Mismatches, stats.BothEmpty, stats.MatchRate()))
	}

	// Category breakdown.
	b.WriteString("\n## Mismatch Categories\n\n")
	b.WriteString("| Category | Count |\n")
	b.WriteString("|----------|------:|\n")
	for cat, count := range report.CategoryCounts {
		if cat == CatMatch || cat == CatBothEmpty {
			continue
		}
		b.WriteString(fmt.Sprintf("| %s | %d |\n", cat, count))
	}

	// Sample mismatches (first 20).
	b.WriteString("\n## Sample Mismatches\n\n")
	shown := 0
	for _, cmp := range report.Comparisons {
		if cmp.AllMatch || cmp.Error != "" {
			continue
		}
		if shown >= 20 {
			break
		}
		b.WriteString(fmt.Sprintf("### %s (transcript: %s)\n\n", cmp.VariantID, cmp.TranscriptID))
		for _, f := range cmp.Fields {
			if f.Match || f.Category == CatBothEmpty {
				continue
			}
			b.WriteString(fmt.Sprintf("- **%s**: GN=`%s` vibe-vep=`%s` [%s]\n",
				f.Field, f.GNValue, f.VEPValue, f.Category))
		}
		b.WriteString("\n")
		shown++
	}

	return b.String()
}

func pct(num, denom int) float64 {
	if denom == 0 {
		return 0
	}
	return float64(num) / float64(denom) * 100.0
}
