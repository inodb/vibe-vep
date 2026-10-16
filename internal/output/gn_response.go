package output

import (
	"encoding/json"
	"sort"
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/datasource/revue"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// GNMarshalOptions controls which fields to include in the GN response.
type GNMarshalOptions struct {
	IncludeAnnotationSummary bool
	IncludeClinVar           bool
	IncludeHotspots          bool
	IncludeSignal            bool
	IncludeMyVariantInfo     bool

	// OverrideTranscripts is the set of transcript ids the isoform override
	// source names canonical, across all genes. Empty leaves the annotations'
	// own canonical flags alone.
	OverrideTranscripts map[string]bool

	// CuratedGenes and CancerGenes are the OncoKB curated genes and the whole
	// cancer gene list. They decide which gene a variant falling in two of
	// them is reported against.
	CuratedGenes map[string]bool
	CancerGenes  map[string]bool

	// OfficialGeneSymbols maps a retired gene symbol to the approved one.
	// The MAF names the gene by its approved symbol even where the GENCODE
	// release the annotation came from still uses the old one.
	OfficialGeneSymbols map[string]string

	// UniprotIDs maps a transcript to the protein entry it translates to.
	UniprotIDs func(transcriptID string) string

	// IncludePtms adds the post-translational modification sites on each
	// transcript, and PtmLookup supplies them.
	IncludePtms bool
	PtmLookup   func(transcriptID string) []GNPtmEntry

	// Revue is the catalogue of variants whose published protein effect is
	// not the one the sequence predicts. Nil leaves every prediction alone.
	Revue *revue.Store

	// RevueConfirmedOnly restricts those replacements to the entries
	// published evidence supports, leaving the ones still under review out.
	RevueConfirmedOnly bool

	// MyVariantInfoData holds pre-fetched myvariant.info data to include in the response.
	// The handler populates this before calling MarshalGNAnnotation.
	MyVariantInfoData *GNMyVariantInfoAnnotation
}

// MarshalGNAnnotation builds a GNAnnotation from a variant and its annotations,
// then marshals it to JSON. This produces a genome-nexus compatible response.
func MarshalGNAnnotation(input string, v *vcf.Variant, anns []*annotate.Annotation, assembly string, opts ...GNMarshalOptions) ([]byte, error) {
	// Order the annotations the way VEP orders transcript consequences,
	// because the choice of canonical transcript reads back to it: where two
	// candidates are equally severe, the first one listed wins. CDKN2A
	// 9:21971036 reaches both ENST00000304494 and ENST00000579755 as
	// missense, both are MSKCC canonical, and the order alone decides which
	// protein change the MAF reports.
	anns = sortedByTranscriptID(anns)

	ref, alt := alleleStrings(v)

	start, end := mafCoordinates(v)

	variant := formatHGVSg(v.Chrom, start, end, ref, alt)
	hgvsg := formatHGVSgWithReference(v.Chrom, start, end, ref, alt)

	result := GNAnnotation{
		Variant:               variant,
		OriginalVariantQuery:  input,
		HGVSg:                 hgvsg,
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
	canonicalAnn, canonicalIDs := resolveCanonical(anns, opt, result.MostSevereConsequence)

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
		tc.AlphaMissense = alphaMissenseFrom(ann)

		result.TranscriptConsequences = append(result.TranscriptConsequences, tc)
	}

	if opt.IncludeAnnotationSummary {
		result.AnnotationSummary = buildAnnotationSummary(&result, anns, canonicalAnn, assembly, opt)
	}

	// PTM sites, one list per transcript consequence in the same order.
	if opt.IncludePtms && opt.PtmLookup != nil {
		sites := make([][]GNPtmEntry, 0, len(result.TranscriptConsequences))
		for _, tc := range result.TranscriptConsequences {
			found := opt.PtmLookup(tc.TranscriptID)
			if found == nil {
				found = []GNPtmEntry{}
			}
			sites = append(sites, found)
		}
		result.PtmAnnotation = &GNPtms{Annotation: sites}
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
func buildAnnotationSummary(result *GNAnnotation, anns []*annotate.Annotation, canonical *annotate.Annotation, assembly string, opt GNMarshalOptions) *GNAnnotationSummary {
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
	// The reVUE catalogue is keyed by this, and it is the same for every
	// transcript the variant reaches.
	revueKey := ""
	if opt.Revue != nil {
		revueKey = genomicLocationKey(summary.GenomicLocation)
	}

	summaries := make([]GNTranscriptConsequenceSummary, 0, len(result.TranscriptConsequences))
	var canonicalSummary *GNTranscriptConsequenceSummary
	for i := range result.TranscriptConsequences {
		tc := &result.TranscriptConsequences[i]
		var ann *annotate.Annotation
		if i < len(anns) {
			ann = anns[i]
		}
		tcs := buildTranscriptConsequenceSummary(tc, ann, variantType, inframe, opt)
		if revueKey != "" {
			applyRevisedEffect(&tcs, revueKey, opt)
		}
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
	if canonicalSummary != nil {
		summary.AlphaMissense = canonicalSummary.AlphaMissense
	}
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

func buildTranscriptConsequenceSummary(tc *GNTranscriptConsequence, ann *annotate.Annotation, variantType string, inframe bool, opt GNMarshalOptions) GNTranscriptConsequenceSummary {
	aaRef, aaAlt := splitAminoAcids(tc.AminoAcids)

	// A splice site past the stop codon is named for the last residue of the
	// protein it sits beyond, so the length has to come from the transcript.
	proteinLength := 0
	if ann != nil {
		proteinLength = ann.ProteinLength
	}
	hgvspShort := gnResolveHgvspShort(tc, proteinLength, variantType, inframe)

	tcs := GNTranscriptConsequenceSummary{
		TranscriptID:   tc.TranscriptID,
		CodonChange:    tc.Codons,
		AminoAcids:     tc.AminoAcids,
		AminoAcidRef:   aaRef,
		AminoAcidAlt:   aaAlt,
		HugoGeneSymbol: officialGeneSymbol(tc.GeneSymbol, opt.OfficialGeneSymbols),
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
		AlphaMissense:         tc.AlphaMissense,
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

	if opt.UniprotIDs != nil {
		tcs.UniprotID = opt.UniprotIDs(tc.TranscriptID)
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

// officialGeneSymbol reports a gene by its approved symbol. A symbol with no
// known replacement, which is nearly all of them, is returned unchanged.
func officialGeneSymbol(symbol string, officialGeneSymbols map[string]string) string {
	if official, ok := officialGeneSymbols[symbol]; ok {
		return official
	}
	return symbol
}

// sortedByTranscriptID returns the annotations in ascending transcript ID
// order, which is the order Ensembl VEP emits them in. The input is left
// alone.
func sortedByTranscriptID(anns []*annotate.Annotation) []*annotate.Annotation {
	sorted := make([]*annotate.Annotation, len(anns))
	copy(sorted, anns)
	sort.SliceStable(sorted, func(i, j int) bool {
		return sorted[i].TranscriptID < sorted[j].TranscriptID
	})
	return sorted
}

// genomicLocationKey spells a variant the way the reVUE catalogue spells it.
func genomicLocationKey(loc GNGenomicLocation) string {
	return strings.Join([]string{
		loc.Chromosome,
		strconv.FormatInt(loc.Start, 10),
		strconv.FormatInt(loc.End, 10),
		loc.ReferenceAllele,
		loc.VariantAllele,
	}, ",")
}

// applyRevisedEffect replaces a predicted protein change with the published
// one where reVUE has an entry for this variant on this transcript. The
// position is read back out of the revised name, since that is now the only
// place it is recorded.
func applyRevisedEffect(tcs *GNTranscriptConsequenceSummary, genomicLocation string, opt GNMarshalOptions) {
	entry, ok := opt.Revue.Lookup(tcs.TranscriptID, genomicLocation)
	if !ok || (opt.RevueConfirmedOnly && !entry.Confirmed) {
		return
	}
	tcs.VariantClassification = entry.VariantClassification
	tcs.HGVSpShort = entry.ProteinEffect
	tcs.ProteinPosition = extractProteinPos(entry.ProteinEffect)
	tcs.IsVue = true
}

// formatHGVSg names a variant in genomic HGVS. This is the "variant" field,
// which names the change without restating the bases it removes: a deletion
// is g.100_102del rather than g.100_102delACG.
func formatHGVSg(chrom string, start, end int64, ref, alt string) string {
	prefix := chrom + ":g."
	switch {
	case start < 1:
		// Ensembl counts from one, so there is nothing to name here.
		return ""
	case isEmptyAllele(ref):
		return prefix + strconv.FormatInt(start, 10) + "_" + strconv.FormatInt(start+1, 10) + "ins" + alt
	case isEmptyAllele(alt):
		if len(ref) == 1 {
			return prefix + strconv.FormatInt(start, 10) + "del"
		}
		return prefix + strconv.FormatInt(start, 10) + "_" + strconv.FormatInt(end, 10) + "del"
	case len(ref) > 1:
		return prefix + strconv.FormatInt(start, 10) + "_" + strconv.FormatInt(end, 10) + "delins" + alt
	case len(alt) > 1:
		return prefix + strconv.FormatInt(start, 10) + "delins" + alt
	default:
		return prefix + strconv.FormatInt(start, 10) + ref + ">" + alt
	}
}

// formatHGVSgWithReference is the "hgvsg" field, which does restate the bases
// a change removes, so that the notation carries the reference allele the
// annotator worked from.
func formatHGVSgWithReference(chrom string, start, end int64, ref, alt string) string {
	hgvs := formatHGVSg(chrom, start, end, ref, alt)
	switch {
	case hgvs == "" || isEmptyAllele(ref):
		// An insertion removes nothing, so nothing is left out.
		return hgvs
	case isEmptyAllele(alt):
		return hgvs + ref
	default:
		if i := strings.Index(hgvs, "delins"); i >= 0 {
			return hgvs[:i] + "del" + ref + hgvs[i+3:]
		}
		// A substitution already carries both alleles.
		return hgvs
	}
}

// isEmptyAllele reports whether a MAF allele stands for no sequence at all.
func isEmptyAllele(allele string) bool {
	switch allele {
	case "", "-", "NA":
		return true
	}
	return strings.Contains(allele, "--")
}

// alphaMissenseFrom reads AlphaMissense's call off an annotation, or nothing
// where the variant is not a missense substitution it scores.
func alphaMissenseFrom(ann *annotate.Annotation) *GNAlphaMissense {
	class := ann.GetExtraKey("alphamissense.class")
	raw := ann.GetExtraKey("alphamissense.score")
	if class == "" && raw == "" {
		return nil
	}

	am := &GNAlphaMissense{Pathogenicity: class}
	if score, err := strconv.ParseFloat(raw, 64); err == nil {
		am.Score = &score
	}
	return am
}
