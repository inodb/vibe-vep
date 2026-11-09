package output

import (
	"math"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
)

// A MAF row describes exactly one transcript, and every annotation column is
// read off that transcript alone, so picking it is the single most consequential
// decision in the response. This file reproduces how genome-nexus picks it:
// IsoformAnnotationEnricher narrows the field to the transcripts worth
// considering and marks them canonical, then CanonicalTranscriptResolver and
// TranscriptConsequencePrioritizer separate whatever is left. Serving the
// genome-nexus API directly means making the same choice it would have made
// from the same annotations.

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
	"splice_polypyrimidine_tract_variant": 11,
	"synonymous_variant":                  12,
	"start_retained_variant":              12,
	"stop_retained_variant":               12,
	"incomplete_terminal_codon_variant":   13,
	// modifier
	"coding_sequence_variant":                        14,
	"mature_mirna_variant":                           14,
	"exon_variant":                                   14,
	"5_prime_utr_variant":                            15,
	"5_prime_utr_premature_start_codon_gain_variant": 15,
	"3_prime_utr_variant":                            16,
	"non_coding_exon_variant":                        16,
	"non_coding_transcript_exon_variant":             16,
	"non_coding_transcript_variant":                  17,
	"nc_transcript_variant":                          17,
	"intron_variant":                                 17,
	"intragenic_variant":                             17,
	"intragenic":                                     17,
	"nmd_transcript_variant":                         18,
	"coding_transcript_variant":                      19,
	"upstream_gene_variant":                          20,
	"downstream_gene_variant":                        20,
	"tfbs_ablation":                                  21,
	"tfbs_amplification":                             21,
	"tf_binding_site_variant":                        21,
	"regulatory_region_ablation":                     21,
	"regulatory_region_amplification":                21,
	"regulatory_region_variant":                      21,
	"regulatory_region":                              21,
	"intergenic_variant":                             22,
	"intergenic_region":                              22,
	"sequence_variant":                               22,
	"":                                               23,
}

// rankEffect returns the severity rank of a consequence term. An unknown term
// ranks below every known one, which is what genome-nexus's getOrDefault does.
func rankEffect(term string) int {
	if p, ok := effectPriority[strings.ToLower(term)]; ok {
		return p
	}
	return math.MaxInt32
}

// biotypePriority ranks transcript biotypes, lowest first. It is vcf2maf's
// GetBiotypePriority table, which genome-nexus uses to prefer a protein-coding
// transcript over a retained intron or a nonsense-mediated decay copy of the
// same region. An unlisted biotype ranks below every listed one.
var biotypePriority = map[string]int{
	"protein_coding": 1,

	"LRG_gene": 2, "IG_C_gene": 2, "IG_D_gene": 2, "IG_J_gene": 2, "IG_LV_gene": 2,
	"IG_V_gene": 2, "TR_C_gene": 2, "TR_D_gene": 2, "TR_J_gene": 2, "TR_V_gene": 2,

	"miRNA": 3, "snRNA": 3, "snoRNA": 3, "ribozyme": 3, "tRNA": 3, "sRNA": 3,
	"scaRNA": 3, "rRNA": 3, "scRNA": 3, "lincRNA": 3, "lncRNA": 3,
	"bidirectional_promoter_lncrna": 3, "bidirectional_promoter_lncRNA": 3,

	"known_ncrna": 4, "vaultRNA": 4, "vault_RNA": 4, "macro_lncRNA": 4,
	"Mt_tRNA": 4, "Mt_rRNA": 4,

	"antisense": 5, "antisense_RNA": 5, "sense_intronic": 5, "sense_overlapping": 5,
	"3prime_overlapping_ncrna": 5, "3prime_overlapping_ncRNA": 5, "misc_RNA": 5,
	"non_coding": 5,

	"regulatory_region": 6, "disrupted_domain": 6, "processed_transcript": 6,
	"protein_coding_CDS_not_defined": 6, "TEC": 6,

	"TF_binding_site": 7, "CTCF_binding_site": 7, "promoter_flanking_region": 7,
	"enhancer": 7, "promoter": 7, "open_chromatin_region": 7, "retained_intron": 7,
	"nonsense_mediated_decay": 7, "non_stop_decay": 7, "ambiguous_orf": 7,

	"pseudogene": 8, "processed_pseudogene": 8, "polymorphic_pseudogene": 8,
	"protein_coding_LoF": 8, "retrotransposed": 8, "translated_processed_pseudogene": 8,
	"translated_unprocessed_pseudogene": 8, "transcribed_processed_pseudogene": 8,
	"transcribed_unprocessed_pseudogene": 8, "transcribed_unitary_pseudogene": 8,
	"unitary_pseudogene": 8, "unprocessed_pseudogene": 8, "Mt_tRNA_pseudogene": 8,
	"tRNA_pseudogene": 8, "snoRNA_pseudogene": 8, "snRNA_pseudogene": 8,
	"scRNA_pseudogene": 8, "rRNA_pseudogene": 8, "misc_RNA_pseudogene": 8,
	"miRNA_pseudogene": 8, "IG_pseudogene": 8, "IG_C_pseudogene": 8,
	"IG_D_pseudogene": 8, "IG_J_pseudogene": 8, "IG_V_pseudogene": 8,
	"TR_J_pseudogene": 8, "TR_V_pseudogene": 8,

	"artifact": 9,
}

// rankBiotype returns a transcript biotype's priority.
func rankBiotype(biotype string) int {
	if p, ok := biotypePriority[biotype]; ok {
		return p
	}
	return 10
}

// geneKey identifies the gene a transcript belongs to, by id where there is
// one and by symbol otherwise.
func geneKey(ann *annotate.Annotation) string {
	if ann.GeneID != "" {
		return ann.GeneID
	}
	return ann.GeneName
}

// resolveCanonical picks the transcript a MAF row will describe, and reports
// which transcripts ended up marked canonical along the way.
//
// Without an override table there is nothing to decide and the annotations
// keep the canonical flags they arrived with. With one, the field is narrowed
// four times -- to the genes with the best biotype, then to the cancer genes
// among those, then to the best biotype within them, then to the override's
// own choice -- and whatever survives is canonical. Each step only ever
// narrows, and a step that would empty the field is skipped, so there is
// always something left to choose from.
func resolveCanonical(anns []*annotate.Annotation, opt GNMarshalOptions, mostSevere string) (*annotate.Annotation, map[string]bool) {
	if len(opt.OverrideTranscripts) == 0 {
		var flagged []*annotate.Annotation
		for _, ann := range anns {
			if ann.IsCanonicalEnsembl {
				flagged = append(flagged, ann)
			}
		}
		return pickAmongCandidates(anns, flagged, mostSevere)
	}

	candidates := bestBiotypeGenes(anns)
	candidates = preferCancerGenes(candidates, opt.CuratedGenes, opt.CancerGenes)
	candidates = bestBiotypeTranscripts(candidates)

	var overridden []*annotate.Annotation
	for _, ann := range candidates {
		if opt.OverrideTranscripts[stripVersion(ann.TranscriptID)] {
			overridden = append(overridden, ann)
		}
	}
	if len(overridden) > 0 {
		candidates = overridden
	}

	return pickAmongCandidates(anns, candidates, mostSevere)
}

// pickAmongCandidates is genome-nexus's CanonicalTranscriptResolver: one
// candidate decides by itself, several are separated by the variant's most
// severe consequence, and none falls back to the whole field.
func pickAmongCandidates(anns, candidates []*annotate.Annotation, mostSevere string) (*annotate.Annotation, map[string]bool) {
	marked := make(map[string]bool, len(candidates))
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

// bestBiotypeGenes keeps the transcripts of the genes that reach the best
// biotype anyone reaches, so a protein-coding gene is preferred over a gene
// that only has non-coding transcripts here.
func bestBiotypeGenes(anns []*annotate.Annotation) []*annotate.Annotation {
	best := make(map[string]int)
	for _, ann := range anns {
		key := geneKey(ann)
		if key == "" {
			continue
		}
		p := rankBiotype(ann.Biotype)
		if cur, ok := best[key]; !ok || p < cur {
			best[key] = p
		}
	}
	if len(best) == 0 {
		return anns
	}

	globalBest := 10
	first := true
	for _, p := range best {
		if first || p < globalBest {
			globalBest = p
			first = false
		}
	}

	var kept []*annotate.Annotation
	for _, ann := range anns {
		if best[geneKey(ann)] == globalBest {
			kept = append(kept, ann)
		}
	}
	return kept
}

// bestBiotypeTranscripts keeps the transcripts whose own biotype is the best
// among those left.
func bestBiotypeTranscripts(anns []*annotate.Annotation) []*annotate.Annotation {
	if len(anns) == 0 {
		return anns
	}
	best := 10
	for i, ann := range anns {
		if p := rankBiotype(ann.Biotype); i == 0 || p < best {
			best = p
		}
	}
	var kept []*annotate.Annotation
	for _, ann := range anns {
		if rankBiotype(ann.Biotype) == best {
			kept = append(kept, ann)
		}
	}
	return kept
}

// preferCancerGenes narrows a variant that falls in more than one gene to the
// gene a cancer panel is actually interested in -- a curated OncoKB gene
// first, and failing that any gene on the cancer gene list. RTEL1 and the
// RTEL1-TNFRSF6B readthrough both carry 20:62319056, and this is what leaves
// the row describing RTEL1.
func preferCancerGenes(anns []*annotate.Annotation, curated, cancer map[string]bool) []*annotate.Annotation {
	for _, genes := range []map[string]bool{curated, cancer} {
		if len(genes) == 0 {
			continue
		}
		var kept []*annotate.Annotation
		for _, ann := range anns {
			if ann.GeneName != "" && genes[ann.GeneName] {
				kept = append(kept, ann)
			}
		}
		if len(kept) > 0 {
			return kept
		}
	}
	return anns
}

// transcriptWithMostSevereConsequence breaks a tie between transcripts. The
// transcripts carrying the variant's most severe consequence term are
// preferred, and among those the ones reaching the most severe term anyone
// reaches. Only if that still leaves several does the first one win, so the
// order the annotations arrive in decides as little as possible.
func transcriptWithMostSevereConsequence(anns []*annotate.Annotation, mostSevere string) *annotate.Annotation {
	if len(anns) == 0 {
		return nil
	}

	target := strings.ToLower(strings.TrimSpace(mostSevere))
	var candidates []*annotate.Annotation
	if target != "" {
		for _, ann := range anns {
			for _, term := range splitConsequence(ann.Consequence) {
				if strings.ToLower(strings.TrimSpace(term)) == target {
					candidates = append(candidates, ann)
					break
				}
			}
		}
	}
	if len(candidates) == 0 {
		candidates = anns
	}

	best := math.MaxInt32
	for _, ann := range candidates {
		for _, term := range splitConsequence(ann.Consequence) {
			if p := rankEffect(term); p < best {
				best = p
			}
		}
	}
	for _, ann := range candidates {
		for _, term := range splitConsequence(ann.Consequence) {
			if rankEffect(term) == best {
				return ann
			}
		}
	}
	return candidates[0]
}
