package output

import (
	"encoding/json"
	"fmt"
	"sort"
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// GNMarshalOptions controls which fields to include in the GN response.
type GNMarshalOptions struct {
	IncludeAnnotationSummary bool
	IncludeClinVar           bool
	IncludeHotspots          bool
	IncludeSignal            bool
	IncludeMyVariantInfo     bool

	// IsoformOverrideSource names the isoform override table that decides
	// which transcript is canonical ("mskcc"; anything else falls back to
	// Ensembl's own flag). It is the API's isoformOverrideSource parameter.
	IsoformOverrideSource string

	// MyVariantInfoData holds pre-fetched myvariant.info data to include in the response.
	// The handler populates this before calling MarshalGNAnnotation.
	MyVariantInfoData *GNMyVariantInfoAnnotation
}

// MarshalGNAnnotation builds a GNAnnotation from a variant and its annotations,
// then marshals it to JSON. This produces a genome-nexus compatible response.
func MarshalGNAnnotation(input string, v *vcf.Variant, anns []*annotate.Annotation, assembly string, opts ...GNMarshalOptions) ([]byte, error) {
	ref, alt := alleleStrings(v)

	start, end := mafCoordinates(v)

	// Build HGVSg-style variant notation for the "variant" field.
	variant := fmt.Sprintf("%s:g.%d%s>%s", v.Chrom, v.Pos, ref, alt)

	result := GNAnnotation{
		Variant:               variant,
		OriginalVariantQuery:  input,
		HGVSg:                 variant,
		ID:                    annotate.FormatVariantID(v.Chrom, v.Pos, v.Ref, v.Alt),
		AssemblyName:          assembly,
		SeqRegionName:         v.Chrom,
		Start:                 start,
		End:                   end,
		AlleleString:          ref + "/" + alt,
		Strand:                1,
		SuccessfullyAnnotated: true,
	}

	var opt GNMarshalOptions
	if len(opts) > 0 {
		opt = opts[0]
	}

	// Every annotation column below is read off one transcript, chosen the way
	// genome-nexus chooses it: the isoform override table names the candidates,
	// and the variant's most severe consequence breaks a tie between them.
	result.MostSevereConsequence = mostSevereConsequence(anns)
	canonicalAnn, canonicalIDs := resolveCanonical(anns, opt.IsoformOverrideSource, result.MostSevereConsequence)

	for _, ann := range anns {
		canonical := ""
		if canonicalIDs[ann.TranscriptID] {
			canonical = "1"
		}

		pStart, pEnd := proteinRange(ann)
		tc := GNTranscriptConsequence{
			TranscriptID:        stripVersion(ann.TranscriptID),
			GeneSymbol:          ann.GeneName,
			GeneID:              ann.GeneID,
			HGNCId:              ann.HGNCId,
			ProteinID:           ann.ProteinID,
			ConsequenceTerms:    splitConsequence(ann.Consequence),
			Impact:              ann.Impact,
			VariantAllele:       ann.Allele,
			AminoAcids:          aminoAcidsField(ann),
			Codons:              ann.CodonChange,
			ProteinStart:        pStart,
			ProteinEnd:          pEnd,
			CDSStart:            ann.CDSPosition,
			CDSEnd:              ann.CDSPosition,
			CDNAStart:           ann.CDNAPosition,
			CDNAEnd:             ann.CDNAPosition,
			HGVSp:               ann.HGVSp,
			HGVSc:               prefixTranscript(ann.TranscriptID, ann.HGVSc),
			Exon:                ann.ExonNumber,
			Intron:              ann.IntronNumber,
			Biotype:             ann.Biotype,
			Canonical:           canonical,
			RefseqTranscriptIds: ann.RefSeqIDs,
			HGVSOffset:          ann.HGVSOffset,
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

	if opt.IncludeAnnotationSummary {
		result.AnnotationSummary = buildAnnotationSummary(&result, anns, canonicalAnn, assembly)
	}

	// Sort transcript consequences: canonical first (genome-nexus convention).
	sortCanonicalFirst(result.TranscriptConsequences)

	// Source fields are populated from annotation extras (same for all transcripts).
	src := firstAnnotationWithExtras(anns)

	// dbSNP → colocatedVariants (from local index or myvariant.info).
	if src != nil {
		if rsID := src.GetExtraKey("dbsnp.id"); rsID != "" {
			result.ColocatedVariants = []GNColocatedVariant{{DbSnpID: rsID}}
		}
	}
	// Fall back to myvariant.info dbSNP rsid for colocatedVariants.
	if len(result.ColocatedVariants) == 0 && opt.MyVariantInfoData != nil &&
		opt.MyVariantInfoData.Annotation != nil && opt.MyVariantInfoData.Annotation.Dbsnp != nil {
		rsid := opt.MyVariantInfoData.Annotation.Dbsnp.Rsid
		if rsid != "" {
			result.ColocatedVariants = []GNColocatedVariant{{DbSnpID: rsid}}
		}
	}

	// ClinVar
	if opt.IncludeClinVar && src != nil {
		if clnsig := src.GetExtraKey("clinvar.clnsig"); clnsig != "" {
			result.ClinVar = &GNClinVar{
				Annotation: &GNClinVarAnnotation{
					Chromosome:           v.Chrom,
					StartPosition:        v.Pos,
					EndPosition:          end,
					ReferenceAllele:      ref,
					AlternateAllele:      alt,
					ClinicalSignificance: clnsig,
					ReviewStatus:         src.GetExtraKey("clinvar.clnrevstat"),
					DiseaseName:          src.GetExtraKey("clinvar.clndn"),
				},
			}
		}
	}

	// Hotspots — build per-transcript arrays matching genome-nexus format.
	if opt.IncludeHotspots {
		hotspotAnnotation := make([][]GNHotspotEntry, 0, len(anns))
		for _, ann := range anns {
			if ann.GetExtraKey("hotspots.hotspot") == "Y" {
				hotspotAnnotation = append(hotspotAnnotation, []GNHotspotEntry{{
					HugoSymbol:   ann.GeneName,
					TranscriptID: stripVersion(ann.TranscriptID),
					Residue:      formatResidue(ann.AminoAcidChange, ann.ProteinPosition),
					Type:         ann.GetExtraKey("hotspots.type"),
				}})
			} else {
				hotspotAnnotation = append(hotspotAnnotation, []GNHotspotEntry{})
			}
		}
		result.Hotspots = &GNHotspots{
			License:    "https://opendatacommons.org/licenses/odbl/1.0/",
			Annotation: hotspotAnnotation,
		}
	}

	// SIGNAL
	if opt.IncludeSignal && src != nil {
		if mutStatus := src.GetExtraKey("signal.mutation_status"); mutStatus != "" {
			// Find gene name from canonical annotation.
			geneName := ""
			if canonicalAnn != nil {
				geneName = canonicalAnn.GeneName
			}
			result.SignalAnnotation = &GNSignalAnnotation{
				License: "https://www.signaldb.org/about",
				Annotation: []GNSignalMutation{{
					HugoGeneSymbol:  geneName,
					Chromosome:      v.Chrom,
					StartPosition:   v.Pos,
					EndPosition:     end,
					ReferenceAllele: ref,
					VariantAllele:   alt,
					MutationStatus:  mutStatus,
				}},
			}
		} else {
			// Return empty annotation array (not null) so frontend doesn't show N/A.
			result.SignalAnnotation = &GNSignalAnnotation{
				License:    "https://www.signaldb.org/about",
				Annotation: []GNSignalMutation{},
			}
		}
	}

	// MyVariantInfo (pre-fetched by handler).
	if opt.IncludeMyVariantInfo && opt.MyVariantInfoData != nil {
		result.MyVariantInfo = opt.MyVariantInfoData
	}

	return json.Marshal(result)
}

// buildAnnotationSummary constructs the annotation_summary enrichment.
func buildAnnotationSummary(result *GNAnnotation, anns []*annotate.Annotation, canonical *annotate.Annotation, assembly string) *GNAnnotationSummary {
	ref, alt := splitAlleleString(result.AlleleString)
	variantType := gnVariantType(ref, alt)
	inframe := gnIsInframe(ref, alt)

	summary := &GNAnnotationSummary{
		Variant: result.Variant,
		GenomicLocation: GNGenomicLocation{
			Chromosome:      result.SeqRegionName,
			Start:           result.Start,
			End:             result.End,
			ReferenceAllele: ref,
			VariantAllele:   alt,
		},
		StrandSign:   "+",
		VariantType:  variantType,
		AssemblyName: assembly,
	}

	if canonical != nil {
		summary.CanonicalTranscriptID = stripVersion(canonical.TranscriptID)
	}

	// result.TranscriptConsequences and anns are still index-aligned here: the
	// list is only reordered for display after the summary is built. The
	// pairing is what lets the Entrez gene id, which VEP does not report and
	// genome-nexus looks up for itself, reach the summary.
	summaries := make([]GNTranscriptConsequenceSummary, 0, len(result.TranscriptConsequences))
	var canonicalSummary *GNTranscriptConsequenceSummary
	for i := range result.TranscriptConsequences {
		tc := &result.TranscriptConsequences[i]
		var ann *annotate.Annotation
		if i < len(anns) {
			ann = anns[i]
		}
		tcs := buildTranscriptConsequenceSummary(tc, ann, variantType, inframe)
		summaries = append(summaries, tcs)
		if canonical != nil && ann != nil && ann.TranscriptID == canonical.TranscriptID && canonicalSummary == nil {
			cp := tcs
			canonicalSummary = &cp
		}
	}

	// Sort by transcript ID to match genome-nexus convention.
	sort.Slice(summaries, func(i, j int) bool {
		return summaries[i].TranscriptID < summaries[j].TranscriptID
	})
	summary.TranscriptConsequenceSummaries = summaries
	summary.TranscriptConsequenceSummary = canonicalSummary
	// Deprecated field: only canonical transcript.
	if canonicalSummary != nil {
		summary.TranscriptConsequences = []GNTranscriptConsequenceSummary{*canonicalSummary}
	}

	return summary
}

// splitAlleleString splits the "ref/alt" form genome-nexus reads the variant's
// alleles back out of.
func splitAlleleString(alleleString string) (ref, alt string) {
	parts := strings.SplitN(alleleString, "/", 2)
	if len(parts) != 2 {
		return "", ""
	}
	return parts[0], parts[1]
}

func buildTranscriptConsequenceSummary(tc *GNTranscriptConsequence, ann *annotate.Annotation, variantType string, inframe bool) GNTranscriptConsequenceSummary {
	aaRef, aaAlt := splitAminoAcids(tc.AminoAcids)
	hgvspShort := gnResolveHgvspShort(tc, variantType, inframe)

	tcs := GNTranscriptConsequenceSummary{
		TranscriptID:   tc.TranscriptID,
		CodonChange:    tc.Codons,
		AminoAcids:     tc.AminoAcids,
		AminoAcidRef:   aaRef,
		AminoAcidAlt:   aaAlt,
		HugoGeneSymbol: tc.GeneSymbol,
		HGVSpShort:     hgvspShort,
		HGVSp:          gnResolveHgvsp(tc),
		HGVSc:          tc.HGVSc,
		// Every term, comma-joined: this is the MAF Consequence column, which
		// names all of them and not just the most severe.
		ConsequenceTerms:      strings.Join(tc.ConsequenceTerms, ","),
		VariantClassification: gnVariantClassification(pickHighestPriorityConsequence(tc.ConsequenceTerms), variantType, inframe),
		Exon:                  tc.Exon,
		ProteinPosition:       gnResolveProteinPosition(tc, hgvspShort),
		SIFTScore:             tc.SIFTScore,
		SIFTPrediction:        tc.SIFTPrediction,
		PolyphenScore:         tc.PolyPhenScore,
		PolyphenPrediction:    tc.PolyPhenPrediction,
	}

	// genome-nexus's RefSeqResolver reports the first RefSeq accession only.
	if len(tc.RefseqTranscriptIds) > 0 {
		tcs.RefSeq = tc.RefseqTranscriptIds[0]
	}

	// The Entrez gene id has no VEP-level field to come from; genome-nexus
	// looks it up by gene, and vibe-vep carries it on the annotation.
	if ann != nil {
		tcs.EntrezGeneID = ann.EntrezGeneID
	}

	return tcs
}

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

// formatResidue builds a residue string like "T790" from amino acid change and position.
func formatResidue(aaChange string, proteinPos int64) string {
	if aaChange == "" || proteinPos == 0 {
		return ""
	}
	// aaChange is like "T790M" — take the first character (ref AA) + position.
	if len(aaChange) > 0 {
		return string(aaChange[0]) + strconv.FormatInt(proteinPos, 10)
	}
	return ""
}

// prefixTranscript prepends the transcript ID to an HGVSc/HGVSp notation
// (e.g. "c.2369C>T" → "ENST00000275493.7:c.2369C>T").
func prefixTranscript(txID, notation string) string {
	if notation == "" || txID == "" {
		return notation
	}
	// Don't double-prefix if already has a transcript prefix.
	if strings.Contains(notation, ":") {
		return notation
	}
	return txID + ":" + notation
}

// hgvspStripTranscript strips the transcript prefix from HGVSp,
// returning just the protein-level notation (e.g. "p.Val640Glu").
func hgvspStripTranscript(hgvsp string) string {
	if idx := strings.Index(hgvsp, ":p."); idx >= 0 {
		return hgvsp[idx+1:]
	}
	return hgvsp
}

// firstAnnotationWithExtras returns the first annotation that has extras populated.
func firstAnnotationWithExtras(anns []*annotate.Annotation) *annotate.Annotation {
	for _, ann := range anns {
		if len(ann.Extra) > 0 {
			return ann
		}
	}
	if len(anns) > 0 {
		return anns[0]
	}
	return nil
}

// sortCanonicalFirst sorts transcript consequences so the highest-impact
// canonical transcript comes first, matching genome-nexus convention.
// The frontend reads [0] for SIFT/PolyPhen scores and gene symbol.
func sortCanonicalFirst(tcs []GNTranscriptConsequence) {
	sort.SliceStable(tcs, func(i, j int) bool {
		ci := tcs[i].Canonical == "1"
		cj := tcs[j].Canonical == "1"
		if ci != cj {
			return ci
		}
		// Among canonicals, prefer higher impact (protein-coding > non-coding).
		if ci && cj {
			return annotate.ImpactRank(tcs[i].Impact) > annotate.ImpactRank(tcs[j].Impact)
		}
		return false
	})
}

// splitAminoAcids splits "V/E" into ref "V" and alt "E".
func splitAminoAcids(aa string) (ref, alt string) {
	parts := strings.SplitN(aa, "/", 2)
	if len(parts) == 2 {
		return parts[0], parts[1]
	}
	return aa, ""
}
