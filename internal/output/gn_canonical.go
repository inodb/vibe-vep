package output

import (
	"math"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
)

// A MAF row describes exactly one transcript, and every annotation column is
// read off that transcript alone, so picking it is the single most consequential
// decision in the response. genome-nexus makes that decision in two steps, and
// this file reproduces both: IsoformAnnotationEnricher, which re-marks the
// canonical flag from an isoform override table, and CanonicalTranscriptResolver
// with TranscriptConsequencePrioritizer, which breaks a tie between several
// canonical transcripts. Serving the genome-nexus API directly means making the
// same choice genome-nexus would have made from the same annotations.

// effectPriority ranks Sequence Ontology terms by severity, lowest number most
// severe. The table is genome-nexus's TranscriptConsequencePrioritizer map,
// which in turn follows vcf2maf's GetEffectPriority. It is deliberately not
// annotate.MoreSevereTerm: that ranks terms the way Ensembl does, for
// most_severe_consequence, while this one only ever breaks ties between
// transcripts and has its own idea of the order.
var effectPriority = map[string]int{
	// high
	"transcript_ablation":      1,
	"exon_loss_variant":        1,
	"splice_donor_variant":     2,
	"splice_acceptor_variant":  2,
	"stop_gained":              3,
	"frameshift_variant":       3,
	"stop_lost":                3,
	"start_lost":               4,
	"initiator_codon_variant":  4,
	"transcript_amplification": 5,
	"feature_elongation":       6,
	"feature_truncation":       6,
	// moderate
	"disruptive_inframe_insertion":  7,
	"disruptive_inframe_deletion":   7,
	"inframe_insertion":             7,
	"inframe_deletion":              7,
	"protein_altering_variant":      8,
	"missense_variant":              9,
	"conservative_missense_variant": 9,
	"rare_amino_acid_variant":       9,
	// low
	"splice_donor_5th_base_variant":       10,
	"splice_region_variant":               10,
	"splice_donor_region_variant":         10,
	"splice_polypyrimidine_tract_variant": 10,
	"synonymous_variant":                  11,
	"start_retained_variant":              11,
	"stop_retained_variant":               11,
	"incomplete_terminal_codon_variant":   12,
	// modifier
	"coding_sequence_variant":                        13,
	"mature_mirna_variant":                           13,
	"exon_variant":                                   13,
	"5_prime_utr_variant":                            14,
	"5_prime_utr_premature_start_codon_gain_variant": 14,
	"3_prime_utr_variant":                            15,
	"non_coding_exon_variant":                        15,
	"non_coding_transcript_exon_variant":             15,
	"non_coding_transcript_variant":                  16,
	"nc_transcript_variant":                          16,
	"intron_variant":                                 16,
	"intragenic_variant":                             16,
	"intragenic":                                     16,
	"nmd_transcript_variant":                         17,
	"coding_transcript_variant":                      18,
	"upstream_gene_variant":                          19,
	"downstream_gene_variant":                        19,
	"tfbs_ablation":                                  20,
	"tfbs_amplification":                             20,
	"tf_binding_site_variant":                        20,
	"regulatory_region_ablation":                     20,
	"regulatory_region_amplification":                20,
	"regulatory_region_variant":                      20,
	"regulatory_region":                              20,
	"intergenic_variant":                             21,
	"intergenic_region":                              21,
	"sequence_variant":                               21,
	"":                                               22,
}

// rankEffect returns the severity rank of a consequence term. An unknown term
// ranks below every known one, which is what genome-nexus's getOrDefault does.
func rankEffect(term string) int {
	if p, ok := effectPriority[strings.ToLower(term)]; ok {
		return p
	}
	return math.MaxInt32
}

// isoformCanonical reports whether the isoform override table named by source
// marks this transcript canonical. genome-nexus defaults to uniprot when no
// source is given; vibe-vep only carries the MSKCC and Ensembl tables, so every
// other source falls back to Ensembl's own canonical flag.
func isoformCanonical(ann *annotate.Annotation, source string) bool {
	if strings.EqualFold(source, "mskcc") {
		return ann.IsCanonicalMSK
	}
	return ann.IsCanonicalEnsembl
}

// resolveCanonical picks the transcript a MAF row will describe, and reports
// which transcripts the isoform override marked canonical along the way.
//
// Where the override table knows at least one of the transcripts, its choice
// replaces the canonical flags wholesale; where it knows none of them, the
// annotations keep whatever canonical flag they arrived with.
func resolveCanonical(anns []*annotate.Annotation, source, mostSevere string) (*annotate.Annotation, map[string]bool) {
	var candidates []*annotate.Annotation
	for _, ann := range anns {
		if isoformCanonical(ann, source) {
			candidates = append(candidates, ann)
		}
	}

	marked := make(map[string]bool, len(candidates))
	if len(candidates) == 0 {
		// No override: fall back to the canonical flags already on the
		// annotations, exactly as genome-nexus falls back to VEP's.
		for _, ann := range anns {
			if ann.IsCanonicalEnsembl {
				candidates = append(candidates, ann)
			}
		}
	}
	for _, ann := range candidates {
		marked[ann.TranscriptID] = true
	}

	switch len(candidates) {
	case 1:
		return candidates[0], marked
	case 0:
		return transcriptWithMostSevereConsequence(anns, mostSevere), marked
	default:
		return transcriptWithMostSevereConsequence(candidates, mostSevere), marked
	}
}

// transcriptWithMostSevereConsequence breaks a tie between transcripts. A
// transcript that carries the variant's most severe consequence term wins
// outright; otherwise the most severe term anyone carries decides, and the
// first transcript to reach it wins. Order therefore matters, and it is the
// order the annotations were produced in.
func transcriptWithMostSevereConsequence(anns []*annotate.Annotation, mostSevere string) *annotate.Annotation {
	highest := math.MaxInt32
	var best *annotate.Annotation

	for _, ann := range anns {
		for _, term := range splitConsequence(ann.Consequence) {
			if p := rankEffect(term); p < highest {
				best = ann
				highest = p
			}
			if strings.TrimSpace(term) == mostSevere {
				return ann
			}
		}
	}

	if best != nil {
		return best
	}
	if len(anns) > 0 {
		return anns[0]
	}
	return nil
}
