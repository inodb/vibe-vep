package output

import (
	"bufio"
	"encoding/json"
	"io"
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// VEPTranscriptConsequence represents a transcript consequence in VEP JSON format.
type VEPTranscriptConsequence struct {
	TranscriptID     string   `json:"transcript_id"`
	GeneID           string   `json:"gene_id,omitempty"`
	GeneSymbol       string   `json:"gene_symbol,omitempty"`
	GeneSymbolSource string   `json:"gene_symbol_source,omitempty"`
	Biotype          string   `json:"biotype,omitempty"`
	ConsequenceTerms []string `json:"consequence_terms"`
	Impact           string   `json:"impact"`
	VariantAllele    string   `json:"variant_allele"`

	// Protein fields
	AminoAcids   string `json:"amino_acids,omitempty"`
	Codons       string `json:"codons,omitempty"`
	ProteinStart int64  `json:"protein_start,omitempty"`
	ProteinEnd   int64  `json:"protein_end,omitempty"`

	// Position fields
	CDSStart  int64 `json:"cds_start,omitempty"`
	CDSEnd    int64 `json:"cds_end,omitempty"`
	CDNAStart int64 `json:"cdna_start,omitempty"`
	CDNAEnd   int64 `json:"cdna_end,omitempty"`

	// HGVS
	HGVSc string `json:"hgvsc,omitempty"`
	HGVSp string `json:"hgvsp,omitempty"`

	// Exon/Intron
	Exon   string `json:"exon,omitempty"`
	Intron string `json:"intron,omitempty"`

	// Predictions
	SIFTScore          *float64 `json:"sift_score,omitempty"`
	SIFTPrediction     string   `json:"sift_prediction,omitempty"`
	PolyPhenScore      *float64 `json:"polyphen_score,omitempty"`
	PolyPhenPrediction string   `json:"polyphen_prediction,omitempty"`
}

// VEPVariantAnnotation represents the top-level VEP JSON output for one variant.
type VEPVariantAnnotation struct {
	Input                  string                     `json:"input"`
	ID                     string                     `json:"id,omitempty"`
	SeqRegionName          string                     `json:"seq_region_name"`
	Start                  int64                      `json:"start"`
	End                    int64                      `json:"end"`
	AlleleString           string                     `json:"allele_string"`
	Strand                 int                        `json:"strand"`
	AssemblyName           string                     `json:"assembly_name"`
	MostSevereConsequence  string                     `json:"most_severe_consequence"`
	TranscriptConsequences []VEPTranscriptConsequence `json:"transcript_consequences"`
}

// VibeVepTranscriptConsequence represents a transcript consequence in vibe-vep native format.
type VibeVepTranscriptConsequence struct {
	TranscriptID          string            `json:"transcript_id"`
	HugoSymbol            string            `json:"hugo_symbol,omitempty"`
	GeneID                string            `json:"gene_id,omitempty"`
	Consequence           string            `json:"consequence"`
	Impact                string            `json:"impact"`
	VariantClassification string            `json:"variant_classification,omitempty"`
	HGVSc                 string            `json:"hgvsc,omitempty"`
	HGVSp                 string            `json:"hgvsp,omitempty"`
	HGVSpShort            string            `json:"hgvsp_short,omitempty"`
	ProteinPosition       int64             `json:"protein_position,omitempty"`
	CDSPosition           int64             `json:"cds_position,omitempty"`
	CDNAPosition          int64             `json:"cdna_position,omitempty"`
	AminoAcidChange       string            `json:"amino_acid_change,omitempty"`
	Codons                string            `json:"codons,omitempty"`
	Exon                  string            `json:"exon,omitempty"`
	Intron                string            `json:"intron,omitempty"`
	Biotype               string            `json:"biotype,omitempty"`
	CanonicalMSKCC        bool              `json:"canonical_mskcc,omitempty"`
	CanonicalEnsembl      bool              `json:"canonical_ensembl,omitempty"`
	CanonicalMANE         bool              `json:"canonical_mane,omitempty"`
	Extra                 map[string]string `json:"extra,omitempty"`
}

// VibeVepVariantAnnotation is the top-level vibe-vep native JSON output.
type VibeVepVariantAnnotation struct {
	Input                  string                          `json:"input"`
	Chromosome             string                          `json:"chromosome"`
	Start                  int64                           `json:"start"`
	End                    int64                           `json:"end"`
	ReferenceAllele        string                          `json:"reference_allele"`
	VariantAllele          string                          `json:"variant_allele"`
	Assembly               string                          `json:"assembly"`
	MostSevereConsequence  string                          `json:"most_severe_consequence"`
	TranscriptConsequences []VibeVepTranscriptConsequence  `json:"transcript_consequences"`
}

// JSONLWriter writes annotations in JSONL format (one JSON line per variant).
type JSONLWriter struct {
	w        *bufio.Writer
	format   string // "ensembl-vep-jsonl" or "vibe-vep-jsonl"
	assembly string
	input    string // current input line for context

	// Buffer annotations for current variant.
	curVariant *vcf.Variant
	curAnns    []*annotate.Annotation
}

// NewJSONLWriter creates a JSONL writer.
// format is "ensembl-vep-jsonl" or "vibe-vep-jsonl".
func NewJSONLWriter(w io.Writer, format, assembly string) *JSONLWriter {
	return &JSONLWriter{
		w:        bufio.NewWriter(w),
		format:   format,
		assembly: assembly,
	}
}

// SetInput sets the original input string for the current variant.
func (j *JSONLWriter) SetInput(input string) {
	j.input = input
}

// WriteHeader is a no-op for JSONL.
func (j *JSONLWriter) WriteHeader() error { return nil }

// Write buffers an annotation for the current variant.
func (j *JSONLWriter) Write(v *vcf.Variant, ann *annotate.Annotation) error {
	// If this is a new variant, flush the previous one.
	if j.curVariant != nil && (v.Chrom != j.curVariant.Chrom || v.Pos != j.curVariant.Pos ||
		v.Ref != j.curVariant.Ref || v.Alt != j.curVariant.Alt) {
		if err := j.flushVariant(); err != nil {
			return err
		}
	}
	j.curVariant = v
	j.curAnns = append(j.curAnns, ann)
	return nil
}

// Flush writes any remaining buffered variant and flushes the underlying writer.
func (j *JSONLWriter) Flush() error {
	if j.curVariant != nil {
		if err := j.flushVariant(); err != nil {
			return err
		}
	}
	return j.w.Flush()
}

func (j *JSONLWriter) flushVariant() error {
	var line []byte
	var err error

	switch j.format {
	case "vibe-vep-jsonl":
		line, err = j.marshalVibeVep()
	default: // ensembl-vep-jsonl
		line, err = j.marshalVEP()
	}
	if err != nil {
		return err
	}

	line = append(line, '\n')
	if _, err := j.w.Write(line); err != nil {
		return err
	}
	// Flush after each variant for real-time pipe response.
	if err := j.w.Flush(); err != nil {
		return err
	}

	j.curVariant = nil
	j.curAnns = nil
	return nil
}

func (j *JSONLWriter) marshalVEP() ([]byte, error) {
	v := j.curVariant
	ref, alt := alleleStrings(v)

	end := v.Pos + int64(len(v.Ref)) - 1
	if end < v.Pos {
		end = v.Pos // insertions
	}

	result := VEPVariantAnnotation{
		Input:         j.input,
		ID:            annotate.FormatVariantID(v.Chrom, v.Pos, v.Ref, v.Alt),
		SeqRegionName: v.Chrom,
		Start:         v.Pos,
		End:           end,
		AlleleString:  ref + "/" + alt,
		Strand:        1,
		AssemblyName:  j.assembly,
	}

	// Find most severe consequence.
	bestImpact := -1
	for _, ann := range j.curAnns {
		impact := annotate.ImpactRank(ann.Impact)
		if impact > bestImpact {
			bestImpact = impact
			result.MostSevereConsequence = firstConsequence(ann.Consequence)
		}
	}

	for _, ann := range j.curAnns {
		tc := VEPTranscriptConsequence{
			TranscriptID:     stripVersion(ann.TranscriptID),
			GeneID:           ann.GeneID,
			GeneSymbol:       ann.GeneName,
			GeneSymbolSource: "HGNC",
			Biotype:          ann.Biotype,
			ConsequenceTerms: splitConsequence(ann.Consequence),
			Impact:           ann.Impact,
			VariantAllele:    ann.Allele,
			AminoAcids:       formatAminoAcidsVEP(ann.AminoAcidChange),
			Codons:           ann.CodonChange,
			ProteinStart:     ann.ProteinPosition,
			ProteinEnd:       ann.ProteinPosition,
			CDSStart:         ann.CDSPosition,
			CDSEnd:           ann.CDSPosition,
			CDNAStart:        ann.CDNAPosition,
			CDNAEnd:          ann.CDNAPosition,
			HGVSc:            ann.HGVSc,
			HGVSp:            ann.HGVSp,
			Exon:             ann.ExonNumber,
			Intron:           ann.IntronNumber,
		}

		// SIFT/PolyPhen from annotation source extras.
		if s := ann.GetExtraKey("sift.score"); s != "" {
			if f, err := strconv.ParseFloat(s, 64); err == nil {
				tc.SIFTScore = &f
			}
		}
		tc.SIFTPrediction = ann.GetExtraKey("sift.prediction")
		if s := ann.GetExtraKey("polyphen.score"); s != "" {
			if f, err := strconv.ParseFloat(s, 64); err == nil {
				tc.PolyPhenScore = &f
			}
		}
		tc.PolyPhenPrediction = ann.GetExtraKey("polyphen.prediction")

		result.TranscriptConsequences = append(result.TranscriptConsequences, tc)
	}

	return json.Marshal(result)
}

func (j *JSONLWriter) marshalVibeVep() ([]byte, error) {
	v := j.curVariant
	ref, alt := alleleStrings(v)

	end := v.Pos + int64(len(v.Ref)) - 1
	if end < v.Pos {
		end = v.Pos
	}

	result := VibeVepVariantAnnotation{
		Input:           j.input,
		Chromosome:      v.Chrom,
		Start:           v.Pos,
		End:             end,
		ReferenceAllele: ref,
		VariantAllele:   alt,
		Assembly:        j.assembly,
	}

	bestImpact := -1
	for _, ann := range j.curAnns {
		impact := annotate.ImpactRank(ann.Impact)
		if impact > bestImpact {
			bestImpact = impact
			result.MostSevereConsequence = firstConsequence(ann.Consequence)
		}
	}

	for _, ann := range j.curAnns {
		tc := VibeVepTranscriptConsequence{
			TranscriptID:          ann.TranscriptID,
			HugoSymbol:            ann.GeneName,
			GeneID:                ann.GeneID,
			Consequence:           ann.Consequence,
			Impact:                ann.Impact,
			VariantClassification: SOToMAFClassification(ann.Consequence, v),
			HGVSc:                ann.HGVSc,
			HGVSp:                ann.HGVSp,
			HGVSpShort:           HGVSpToShort(ann.HGVSp),
			ProteinPosition:      ann.ProteinPosition,
			CDSPosition:          ann.CDSPosition,
			CDNAPosition:         ann.CDNAPosition,
			AminoAcidChange:      ann.AminoAcidChange,
			Codons:               ann.CodonChange,
			Exon:                 ann.ExonNumber,
			Intron:               ann.IntronNumber,
			Biotype:              ann.Biotype,
			CanonicalMSKCC:       ann.IsCanonicalMSK,
			CanonicalEnsembl:     ann.IsCanonicalEnsembl,
			CanonicalMANE:        ann.IsMANESelect,
			Extra:                ann.Extra,
		}
		result.TranscriptConsequences = append(result.TranscriptConsequences, tc)
	}

	return json.Marshal(result)
}

// WriteError writes an error as a JSON line to stdout.
func (j *JSONLWriter) WriteError(input string, errMsg string) error {
	obj := map[string]string{"error": errMsg, "input": input}
	line, err := json.Marshal(obj)
	if err != nil {
		return err
	}
	line = append(line, '\n')
	if _, err := j.w.Write(line); err != nil {
		return err
	}
	return j.w.Flush()
}

// alleleStrings returns ref and alt with "-" for empty alleles (MAF convention).
func alleleStrings(v *vcf.Variant) (string, string) {
	ref := v.Ref
	if ref == "" {
		ref = "-"
	}
	alt := v.Alt
	if alt == "" {
		alt = "-"
	}
	return ref, alt
}

func splitConsequence(c string) []string {
	if c == "" {
		return nil
	}
	return strings.Split(c, ",")
}

func firstConsequence(c string) string {
	if i := strings.IndexByte(c, ','); i >= 0 {
		return c[:i]
	}
	return c
}

func stripVersion(id string) string {
	if i := strings.IndexByte(id, '.'); i >= 0 {
		return id[:i]
	}
	return id
}

// formatAminoAcidsVEP converts "G12C" to "G/C" (VEP format).
func formatAminoAcidsVEP(change string) string {
	if len(change) < 2 {
		return ""
	}
	return change[:1] + "/" + change[len(change)-1:]
}

