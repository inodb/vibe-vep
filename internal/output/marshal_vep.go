package output

import (
	"encoding/json"
	"strconv"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// MarshalVEPAnnotation builds a VEPVariantAnnotation from a variant and its annotations,
// then marshals it to JSON. This is the standalone equivalent of JSONLWriter.marshalVEP().
func MarshalVEPAnnotation(input string, v *vcf.Variant, anns []*annotate.Annotation, assembly string) ([]byte, error) {
	ref, alt := alleleStrings(v)

	end := v.Pos + int64(len(v.Ref)) - 1
	if end < v.Pos {
		end = v.Pos // insertions
	}

	result := VEPVariantAnnotation{
		Input:         input,
		ID:            annotate.FormatVariantID(v.Chrom, v.Pos, v.Ref, v.Alt),
		SeqRegionName: v.Chrom,
		Start:         v.Pos,
		End:           end,
		AlleleString:  ref + "/" + alt,
		Strand:        1,
		AssemblyName:  assembly,
	}

	// Find most severe consequence.
	bestImpact := -1
	for _, ann := range anns {
		impact := annotate.ImpactRank(ann.Impact)
		if impact > bestImpact {
			bestImpact = impact
			result.MostSevereConsequence = firstConsequence(ann.Consequence)
		}
	}

	for _, ann := range anns {
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
