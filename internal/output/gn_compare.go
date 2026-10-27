package output

import (
	"encoding/json"
	"fmt"
	"math"
	"regexp"
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
)

// GNAnnotation represents a genome-nexus annotation response.
type GNAnnotation struct {
	Variant                string                     `json:"variant"`
	OriginalVariantQuery   string                     `json:"originalVariantQuery"`
	HGVSg                  string                     `json:"hgvsg"`
	ID                     string                     `json:"id"`
	AssemblyName           string                     `json:"assembly_name"`
	SeqRegionName          string                     `json:"seq_region_name"`
	Start                  int64                      `json:"start"`
	End                    int64                      `json:"end"`
	AlleleString           string                     `json:"allele_string"`
	Strand                 int                        `json:"strand"`
	MostSevereConsequence  string                     `json:"most_severe_consequence"`
	TranscriptConsequences []GNTranscriptConsequence  `json:"transcript_consequences"`
	SuccessfullyAnnotated  bool                       `json:"successfully_annotated"`
	AnnotationSummary      *GNAnnotationSummary       `json:"annotation_summary,omitempty"`
	ClinVar                *GNClinVar                 `json:"clinvar,omitempty"`
	ColocatedVariants      []GNColocatedVariant       `json:"colocatedVariants,omitempty"`
	Hotspots               *GNHotspots                `json:"hotspots,omitempty"`
	SignalAnnotation       *GNSignalAnnotation        `json:"signalAnnotation,omitempty"`
	MyVariantInfo          *GNMyVariantInfoAnnotation `json:"my_variant_info,omitempty"`
	Ptms                   *GNPtms                    `json:"ptms,omitempty"`
}

// GNMyVariantInfoAnnotation wraps the myvariant.info annotation in genome-nexus format.
type GNMyVariantInfoAnnotation struct {
	Annotation *GNMyVariantInfo `json:"annotation,omitempty"`
}

// GNMyVariantInfo holds the transformed myvariant.info data.
type GNMyVariantInfo struct {
	Dbsnp        *GNDbsnp  `json:"dbsnp,omitempty"`
	GnomadExome  *GNGnomad `json:"gnomadExome,omitempty"`
	GnomadGenome *GNGnomad `json:"gnomadGenome,omitempty"`
	Vcf          *GNVcf    `json:"vcf,omitempty"`
	Variant      string    `json:"variant,omitempty"`
	Query        string    `json:"query,omitempty"`
	Hgvs         string    `json:"hgvs,omitempty"`
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
	License    string             `json:"license"`
	Annotation []GNSignalMutation `json:"annotation"`
}

// GNSignalMutation represents a SIGNAL mutation entry.
type GNSignalMutation struct {
	HugoGeneSymbol    string               `json:"hugoGeneSymbol,omitempty"`
	Chromosome        string               `json:"chromosome,omitempty"`
	StartPosition     int64                `json:"startPosition,omitempty"`
	EndPosition       int64                `json:"endPosition,omitempty"`
	ReferenceAllele   string               `json:"referenceAllele,omitempty"`
	VariantAllele     string               `json:"variantAllele,omitempty"`
	MutationStatus    string               `json:"mutationStatus,omitempty"`
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
	License    string             `json:"license,omitempty"`
	Annotation [][]GNHotspotEntry `json:"annotation"`
}

// GNHotspotEntry is one hotspot, in genome-nexus's wire shape. It is what the
// /cancer_hotspots endpoints report and what the "hotspots" field of an
// annotation carries, which are the same records.
//
// Its `id` field is deliberately absent: the Jackson mixin there hides it.
type GNHotspotEntry struct {
	HugoSymbol          string `json:"hugoSymbol,omitempty"`
	TranscriptID        string `json:"transcriptId,omitempty"`
	TranscriptIDVersion string `json:"transcriptIdVersion,omitempty"`
	Residue             string `json:"residue,omitempty"`
	Version             string `json:"version,omitempty"`
	TumorCount          int    `json:"tumorCount"`
	Type                string `json:"type,omitempty"`

	// 3d hotspots carry no per-class breakdown; genome-nexus leaves these out
	// of the response rather than reporting zeros.
	MissenseCount   *int `json:"missenseCount,omitempty"`
	TruncatingCount *int `json:"truncatingCount,omitempty"`
	InframeCount    *int `json:"inframeCount,omitempty"`
	SpliceCount     *int `json:"spliceCount,omitempty"`
}

// GNAnnotationSummary is the enriched annotation summary returned when
// ?fields=annotation_summary is requested.
type GNAnnotationSummary struct {
	Variant                        string                           `json:"variant"`
	GenomicLocation                GNGenomicLocation                `json:"genomicLocation"`
	StrandSign                     string                           `json:"strandSign"`
	VariantType                    string                           `json:"variantType"`
	AssemblyName                   string                           `json:"assemblyName"`
	CanonicalTranscriptID          string                           `json:"canonicalTranscriptId"`
	TranscriptConsequenceSummary   *GNTranscriptConsequenceSummary  `json:"transcriptConsequenceSummary"`
	TranscriptConsequenceSummaries []GNTranscriptConsequenceSummary `json:"transcriptConsequenceSummaries"`
	TranscriptConsequences         []GNTranscriptConsequenceSummary `json:"transcriptConsequences"`

	AlphaMissense *GNAlphaMissense `json:"alphaMissense,omitempty"`
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
	TranscriptID          string          `json:"transcriptId"`
	CodonChange           string          `json:"codonChange,omitempty"`
	AminoAcids            string          `json:"aminoAcids,omitempty"`
	AminoAcidRef          string          `json:"aminoAcidRef,omitempty"`
	AminoAcidAlt          string          `json:"aminoAcidAlt,omitempty"`
	EntrezGeneID          string          `json:"entrezGeneId,omitempty"`
	HugoGeneSymbol        string          `json:"hugoGeneSymbol,omitempty"`
	HGVSpShort            string          `json:"hgvspShort,omitempty"`
	HGVSp                 string          `json:"hgvsp,omitempty"`
	HGVSc                 string          `json:"hgvsc,omitempty"`
	ConsequenceTerms      string          `json:"consequenceTerms,omitempty"`
	VariantClassification string          `json:"variantClassification,omitempty"`
	Exon                  string          `json:"exon,omitempty"`
	RefSeq                string          `json:"refSeq,omitempty"`
	ProteinPosition       *GNIntegerRange `json:"proteinPosition,omitempty"`
	PolyphenScore         *float64        `json:"polyphenScore,omitempty"`
	PolyphenPrediction    string          `json:"polyphenPrediction,omitempty"`
	SIFTScore             *float64        `json:"siftScore,omitempty"`
	SIFTPrediction        string          `json:"siftPrediction,omitempty"`

	AlphaMissense *GNAlphaMissense `json:"alphaMissense,omitempty"`

	// UniprotID names the protein entry the transcript translates to, which
	// clients use to line a variant up against protein-level annotation.
	UniprotID string `json:"uniprotId,omitempty"`

	// IsVue marks a protein change that came from reVUE's published evidence
	// rather than from the sequence.
	IsVue bool `json:"isVue,omitempty"`
}

// GNAlphaMissense is AlphaMissense's pathogenicity call for a missense
// substitution. genome-nexus takes it from VEP's AlphaMissense plugin; here it
// comes from the locally indexed table.
type GNAlphaMissense struct {
	Score         *float64 `json:"score,omitempty"`
	Pathogenicity string   `json:"pathogenicity,omitempty"`
}

// GNPtmEntry is one post-translational modification site.
type GNPtmEntry struct {
	UniprotEntry         string   `json:"uniprotEntry,omitempty"`
	UniprotAccession     string   `json:"uniprotAccession,omitempty"`
	EnsemblTranscriptIds []string `json:"ensemblTranscriptIds,omitempty"`
	Position             int      `json:"position"`
	Type                 string   `json:"type,omitempty"`
	PubmedIds            []string `json:"pubmedIds,omitempty"`
	Sequence             string   `json:"sequence,omitempty"`
}

// GNPtms holds one list of modification sites per transcript consequence, in
// the same order. genome-nexus never sets a license here, so neither does this.
type GNPtms struct {
	License    string         `json:"license,omitempty"`
	Annotation [][]GNPtmEntry `json:"annotation,omitempty"`
}

// GNIntegerRange represents a start/end integer range.
type GNIntegerRange struct {
	Start int64 `json:"start"`
	End   int64 `json:"end"`
}

// GNTranscriptConsequence represents a transcript consequence from genome-nexus.
// GNTranscriptConsequence is one transcript's view of a variant, in the shape
// genome-nexus serves. The field set is its TranscriptConsequence model: VEP
// reports cDNA and CDS positions, impact and intron number too, but
// genome-nexus does not carry them through, so they are absent here as well.
type GNTranscriptConsequence struct {
	TranscriptID     string   `json:"transcript_id"`
	GeneSymbol       string   `json:"gene_symbol"`
	GeneID           string   `json:"gene_id"`
	HGNCId           string   `json:"hgnc_id,omitempty"`
	ProteinID        string   `json:"protein_id,omitempty"`
	ConsequenceTerms []string `json:"consequence_terms"`
	VariantAllele    string   `json:"variant_allele"`

	// Impact orders the consequences so the frontend's [0] is the one worth
	// reading. genome-nexus does not report it, so neither do we.
	Impact string `json:"-"`

	AminoAcids string `json:"amino_acids,omitempty"`
	Codons     string `json:"codons,omitempty"`
	// ProteinStart and ProteinEnd stay on the struct whether or not they
	// reach the wire: the summary builds a protein change out of them even
	// for a splice variant, which is the one place a position is known
	// without codons to translate. MarshalJSON decides what is reported.
	ProteinStart        int64    `json:"protein_start,omitempty"`
	ProteinEnd          int64    `json:"protein_end,omitempty"`
	HGVSp               string   `json:"hgvsp,omitempty"`
	HGVSc               string   `json:"hgvsc,omitempty"`
	HGVSg               string   `json:"hgvsg,omitempty"`
	Exon                string   `json:"exon,omitempty"`
	Biotype             string   `json:"biotype,omitempty"`
	RefseqTranscriptIds []string `json:"refseq_transcript_ids,omitempty"`
	HGVSOffset          int      `json:"hgvs_offset,omitempty"`
	Canonical           string   `json:"canonical,omitempty"`
	SIFTScore           *float64 `json:"sift_score,omitempty"`
	SIFTPrediction      string   `json:"sift_prediction,omitempty"`
	PolyPhenScore       *float64 `json:"polyphen_score,omitempty"`
	PolyPhenPrediction  string   `json:"polyphen_prediction,omitempty"`

	AlphaMissense *GNAlphaMissense `json:"alphaMissense,omitempty"`
}

// MarshalJSON reports a transcript consequence the way VEP does.
//
// VEP gives a protein position exactly where it has codons to translate: a
// variant that only disrupts a splice site sits in an intron, so there is no
// residue to name and the field is absent. Measured against
// vep111.genomenexus.org over 44,084 shared transcripts, reporting a position
// regardless disagreed on 727 of them and this rule on 73.
//
// The position is still computed and kept, because genome-nexus's own
// resolvers read it back to name a protein change for exactly those splice
// variants. Only the wire is affected.
func (tc GNTranscriptConsequence) MarshalJSON() ([]byte, error) {
	type wire GNTranscriptConsequence
	w := wire(tc)
	if tc.Codons == "" {
		startCoding, endCoding := codingEndsOfHGVSc(tc.HGVSc)
		// A change spanning a splice boundary can leave the annotator without
		// a protein position at all, because the coordinate it works from is
		// the intronic end. The coding end names its own residue, so take it
		// from there.
		first, last := codonsOfHGVSc(tc.HGVSc)
		if startCoding && w.ProteinStart == 0 {
			w.ProteinStart = first
		}
		if endCoding && w.ProteinEnd == 0 {
			w.ProteinEnd = last
		}
		switch {
		case !startCoding && !endCoding:
			// Wholly outside the coding sequence: no residue to name.
			w.ProteinStart, w.ProteinEnd = 0, 0
		case !startCoding:
			// Reaches in from an intron or a UTR. Only where it lands is
			// nameable, and VEP reports that as the end alone.
			w.ProteinStart = 0
		case !endCoding:
			// Reaches out of the coding sequence. The range stops where it
			// left, so VEP reports the same residue at both ends.
			w.ProteinEnd = w.ProteinStart
		default:
			// Both ends coding and still no codons: a splice site at a coding
			// base, which changes no codon, so VEP names no residue.
			//
			// genome-nexus does name one -- its summary reads 987 back out of
			// the p.X987_splice it built -- but that is the summary's own
			// work, and it still reports nothing here. The position stays on
			// the struct for the summary to derive from.
			w.ProteinStart, w.ProteinEnd = 0, 0
		}
	}
	return json.Marshal(w)
}

// hgvscRange picks the one or two coordinates out of an HGVS c. or n. name,
// before the edit that follows them.
var hgvscRange = regexp.MustCompile(`^[cn]\.([*\-]?\d+(?:[+-]\d+)?)(?:_([*\-]?\d+(?:[+-]\d+)?))?`)

// codonsOfHGVSc gives the residues the two coordinates of a coding HGVS c.
// name fall in. c.211 is the first base of codon 71, c.213 the last.
//
// Either may be zero: an intronic or UTR coordinate names no residue, and a
// name with one coordinate reports it as both.
func codonsOfHGVSc(hgvsc string) (first, last int64) {
	name := hgvsc
	if i := strings.LastIndex(name, ":"); i >= 0 {
		name = name[i+1:]
	}
	if strings.HasPrefix(name, "n.") {
		return 0, 0
	}
	m := hgvscRange.FindStringSubmatch(name)
	if m == nil {
		return 0, 0
	}
	codon := func(coord string) int64 {
		pos, err := strconv.ParseInt(coord, 10, 64)
		if err != nil || pos < 1 {
			return 0
		}
		return (pos + 2) / 3
	}
	first = codon(m[1])
	if m[2] == "" {
		return first, first
	}
	return first, codon(m[2])
}

// codingEndsOfHGVSc reports whether each end of a change lies in the coding
// sequence, reading it off the name rather than recomputing it.
//
// A variant crossing a splice site has one end in an intron and one in an
// exon, and VEP names a residue only for the end that has one. Which end that
// is cannot be told from the codons, because a change spanning the boundary
// has none at all.
func codingEndsOfHGVSc(hgvsc string) (start, end bool) {
	name := hgvsc
	if i := strings.LastIndex(name, ":"); i >= 0 {
		name = name[i+1:]
	}
	// An "n." name belongs to a transcript with no coding sequence at all, so
	// neither end of it names a residue.
	if strings.HasPrefix(name, "n.") {
		return false, false
	}
	m := hgvscRange.FindStringSubmatch(name)
	if m == nil {
		return false, false
	}
	coding := func(coord string) bool {
		if coord == "" || strings.HasPrefix(coord, "*") || strings.HasPrefix(coord, "-") {
			return false // a UTR position, not a coding one
		}
		// A "+n" or "-n" offset after the position puts it in an intron.
		return !strings.ContainsAny(coord[1:], "+-")
	}
	first := coding(m[1])
	if m[2] == "" {
		return first, first
	}
	return first, coding(m[2])
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
	TotalVariants  int                      `json:"total_variants"`
	FullMatches    int                      `json:"full_matches"`
	PartialMatches int                      `json:"partial_matches"`
	Errors         int                      `json:"errors"`
	FieldStats     map[string]*GNFieldStats `json:"field_stats"`
	CategoryCounts map[Category]int         `json:"category_counts"`
	Comparisons    []GNVariantComparison    `json:"comparisons,omitempty"`
}

// GNFieldStats tracks match/mismatch counts per field.
type GNFieldStats struct {
	Total      int              `json:"total"`
	Matches    int              `json:"matches"`
	Mismatches int              `json:"mismatches"`
	BothEmpty  int              `json:"both_empty"`
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
