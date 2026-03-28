package output

import (
	"encoding/json"
	"fmt"
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// GNMarshalOptions controls which fields to include in the GN response.
type GNMarshalOptions struct {
	IncludeAnnotationSummary bool
}

// MarshalGNAnnotation builds a GNAnnotation from a variant and its annotations,
// then marshals it to JSON. This produces a genome-nexus compatible response.
func MarshalGNAnnotation(input string, v *vcf.Variant, anns []*annotate.Annotation, assembly string, opts ...GNMarshalOptions) ([]byte, error) {
	ref, alt := alleleStrings(v)

	end := v.Pos + int64(len(v.Ref)) - 1
	if end < v.Pos {
		end = v.Pos // insertions
	}

	// Build HGVSg-style variant notation for the "variant" field.
	variant := fmt.Sprintf("%s:g.%d%s>%s", v.Chrom, v.Pos, ref, alt)

	result := GNAnnotation{
		Variant:              variant,
		OriginalVariantQuery: input,
		HGVSg:               variant,
		ID:                   annotate.FormatVariantID(v.Chrom, v.Pos, v.Ref, v.Alt),
		AssemblyName:         assembly,
		SeqRegionName:        v.Chrom,
		Start:                v.Pos,
		End:                  end,
		AlleleString:         ref + "/" + alt,
		Strand:               1,
		SuccessfullyAnnotated: true,
	}

	// Find most severe consequence and canonical transcript.
	bestImpact := -1
	var canonicalAnn *annotate.Annotation
	for _, ann := range anns {
		impact := annotate.ImpactRank(ann.Impact)
		if impact > bestImpact {
			bestImpact = impact
			result.MostSevereConsequence = firstConsequence(ann.Consequence)
		}
		if ann.IsCanonicalEnsembl && canonicalAnn == nil {
			canonicalAnn = ann
		}
	}

	for _, ann := range anns {
		canonical := ""
		if ann.IsCanonicalEnsembl {
			canonical = "1"
		}

		tc := GNTranscriptConsequence{
			TranscriptID:     stripVersion(ann.TranscriptID),
			GeneSymbol:       ann.GeneName,
			GeneID:           ann.GeneID,
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
			HGVSp:            ann.HGVSp,
			HGVSc:            ann.HGVSc,
			Exon:             ann.ExonNumber,
			Intron:           ann.IntronNumber,
			Biotype:          ann.Biotype,
			Canonical:        canonical,
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

	// Build annotation_summary if requested.
	var opt GNMarshalOptions
	if len(opts) > 0 {
		opt = opts[0]
	}
	if opt.IncludeAnnotationSummary {
		result.AnnotationSummary = buildAnnotationSummary(variant, v, anns, canonicalAnn, assembly)
	}

	return json.Marshal(result)
}

// buildAnnotationSummary constructs the annotation_summary enrichment.
func buildAnnotationSummary(variant string, v *vcf.Variant, anns []*annotate.Annotation, canonical *annotate.Annotation, assembly string) *GNAnnotationSummary {
	ref, alt := alleleStrings(v)
	end := v.Pos + int64(len(v.Ref)) - 1
	if end < v.Pos {
		end = v.Pos
	}

	summary := &GNAnnotationSummary{
		Variant: variant,
		GenomicLocation: GNGenomicLocation{
			Chromosome:      v.Chrom,
			Start:           v.Pos,
			End:             end,
			ReferenceAllele: ref,
			VariantAllele:   alt,
		},
		StrandSign:   "+",
		VariantType:  resolveVariantType(ref, alt),
		AssemblyName: assembly,
	}

	if canonical != nil {
		summary.CanonicalTranscriptID = stripVersion(canonical.TranscriptID)
	}

	// Build transcript consequence summaries for all transcripts.
	summaries := make([]GNTranscriptConsequenceSummary, 0, len(anns))
	var canonicalSummary *GNTranscriptConsequenceSummary
	for _, ann := range anns {
		tcs := buildTranscriptConsequenceSummary(ann, v)
		summaries = append(summaries, tcs)
		if ann.IsCanonicalEnsembl && canonicalSummary == nil {
			cp := tcs
			canonicalSummary = &cp
		}
	}

	summary.TranscriptConsequenceSummaries = summaries
	summary.TranscriptConsequenceSummary = canonicalSummary
	// Deprecated field: only canonical transcript.
	if canonicalSummary != nil {
		summary.TranscriptConsequences = []GNTranscriptConsequenceSummary{*canonicalSummary}
	}

	return summary
}

// buildTranscriptConsequenceSummary builds a single transcript consequence summary.
func buildTranscriptConsequenceSummary(ann *annotate.Annotation, v *vcf.Variant) GNTranscriptConsequenceSummary {
	aminoAcids := formatAminoAcidsVEP(ann.AminoAcidChange)
	aaRef, aaAlt := splitAminoAcids(aminoAcids)

	tcs := GNTranscriptConsequenceSummary{
		TranscriptID:          stripVersion(ann.TranscriptID),
		CodonChange:           ann.CodonChange,
		AminoAcids:            aminoAcids,
		AminoAcidRef:          aaRef,
		AminoAcidAlt:          aaAlt,
		HugoGeneSymbol:        ann.GeneName,
		HGVSpShort:            hgvspToShort(ann.HGVSp),
		HGVSp:                 hgvspStripTranscript(ann.HGVSp),
		HGVSc:                 ann.HGVSc,
		ConsequenceTerms:      firstConsequence(ann.Consequence),
		VariantClassification: SOToMAFClassification(ann.Consequence, v),
		Exon:                  ann.ExonNumber,
	}

	if ann.ProteinPosition > 0 {
		tcs.ProteinPosition = &GNIntegerRange{
			Start: ann.ProteinPosition,
			End:   ann.ProteinPosition,
		}
	}

	// SIFT/PolyPhen
	if s := ann.GetExtraKey("sift.score"); s != "" {
		if f, err := strconv.ParseFloat(s, 64); err == nil {
			tcs.SIFTScore = &f
		}
	}
	tcs.SIFTPrediction = ann.GetExtraKey("sift.prediction")
	if s := ann.GetExtraKey("polyphen.score"); s != "" {
		if f, err := strconv.ParseFloat(s, 64); err == nil {
			tcs.PolyphenScore = &f
		}
	}
	tcs.PolyphenPrediction = ann.GetExtraKey("polyphen.prediction")

	return tcs
}

// resolveVariantType determines the variant type from ref/alt alleles.
func resolveVariantType(ref, alt string) string {
	refLen := len(ref)
	altLen := len(alt)
	if ref == "-" {
		refLen = 0
	}
	if alt == "-" {
		altLen = 0
	}

	if refLen == 0 || altLen > refLen {
		return "INS"
	}
	if altLen == 0 || refLen > altLen {
		return "DEL"
	}
	switch refLen {
	case 1:
		return "SNP"
	case 2:
		return "DNP"
	case 3:
		return "TNP"
	default:
		return "ONP"
	}
}

// hgvspToShort converts a long-form HGVSp (e.g. "ENSP00000288602.7:p.Val640Glu")
// to short form (e.g. "p.V640E") using single-letter amino acid codes.
func hgvspToShort(hgvsp string) string {
	if hgvsp == "" {
		return ""
	}
	// Strip transcript prefix if present.
	if idx := strings.Index(hgvsp, ":p."); idx >= 0 {
		hgvsp = hgvsp[idx+1:]
	} else if !strings.HasPrefix(hgvsp, "p.") {
		return hgvsp
	}

	// Replace three-letter amino acid codes with single-letter.
	result := hgvsp
	for three, single := range annotate.AminoAcidThreeToSingle {
		result = strings.ReplaceAll(result, three, string(single))
	}
	// Handle Ter→*
	result = strings.ReplaceAll(result, "Ter", "*")
	return result
}

// hgvspStripTranscript strips the transcript prefix from HGVSp,
// returning just the protein-level notation (e.g. "p.Val640Glu").
func hgvspStripTranscript(hgvsp string) string {
	if idx := strings.Index(hgvsp, ":p."); idx >= 0 {
		return hgvsp[idx+1:]
	}
	return hgvsp
}

// splitAminoAcids splits "V/E" into ref "V" and alt "E".
func splitAminoAcids(aa string) (ref, alt string) {
	parts := strings.SplitN(aa, "/", 2)
	if len(parts) == 2 {
		return parts[0], parts[1]
	}
	return aa, ""
}
