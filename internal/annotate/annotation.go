// Package annotate provides variant effect prediction functionality.
package annotate

import "strings"

// Impact levels for variant consequences.
const (
	ImpactHigh     = "HIGH"
	ImpactModerate = "MODERATE"
	ImpactLow      = "LOW"
	ImpactModifier = "MODIFIER"
)

// Consequence types (Sequence Ontology terms).
const (
	// HIGH impact
	ConsequenceStopGained          = "stop_gained"
	ConsequenceFrameshiftVariant   = "frameshift_variant"
	ConsequenceStopLost            = "stop_lost"
	ConsequenceStartLost           = "start_lost"
	ConsequenceSpliceAcceptor      = "splice_acceptor_variant"
	ConsequenceSpliceDonor         = "splice_donor_variant"

	// MODERATE impact
	ConsequenceMissenseVariant     = "missense_variant"
	ConsequenceInframeInsertion    = "inframe_insertion"
	ConsequenceInframeDeletion     = "inframe_deletion"

	// LOW impact
	ConsequenceSynonymousVariant   = "synonymous_variant"
	ConsequenceSpliceRegion        = "splice_region_variant"
	ConsequenceStopRetained        = "stop_retained_variant"
	ConsequenceStartRetained       = "start_retained_variant"

	// LOW impact (generic coding)
	ConsequenceCodingSequenceVariant = "coding_sequence_variant"

	// Compound consequences (pre-built to avoid runtime concatenation).
	ConsequenceSpliceRegionIntron  = "splice_region_variant,intron_variant"
	ConsequenceStopLost3PrimeUTR   = "stop_lost,3_prime_UTR_variant"
	ConsequenceStopGainedInframeDel = "stop_gained,inframe_deletion"
	ConsequenceFrameshiftStopLost  = "frameshift_variant,stop_lost"

	// MODIFIER impact
	ConsequenceIntronVariant       = "intron_variant"
	Consequence5PrimeUTR           = "5_prime_UTR_variant"
	Consequence3PrimeUTR           = "3_prime_UTR_variant"
	ConsequenceUpstreamGene        = "upstream_gene_variant"
	ConsequenceDownstreamGene      = "downstream_gene_variant"
	ConsequenceIntergenicVariant   = "intergenic_variant"
	ConsequenceNonCodingExon       = "non_coding_transcript_exon_variant"
	ConsequenceMatureMiRNA         = "mature_miRNA_variant"
)

// Annotation represents the predicted effect of a variant on a transcript.
type Annotation struct {
	VariantID       string // Source variant identifier (chrom_pos_ref/alt)
	TranscriptID    string // Affected transcript
	GeneName        string // Gene symbol
	GeneID          string // Gene identifier
	Consequence     string // SO consequence term
	Impact          string // HIGH, MODERATE, LOW, MODIFIER
	CDSPosition     int64  // Position in CDS, 0 if not in CDS
	ProteinPosition int64  // Amino acid position, 0 if not in CDS
	AminoAcidChange string // e.g., "G12C", empty if not missense
	CodonChange     string // e.g., "GGT/TGT", empty if not coding
	IsCanonical     bool   // Annotation on canonical transcript
	Allele          string // The alternate allele
	Biotype         string // Transcript biotype
	ExonNumber      string // Exon number (e.g., "2/5")
	IntronNumber    string // Intron number (e.g., "1/4")
	CDNAPosition    int64  // Position in cDNA
	HGVSp           string // HGVS protein notation (e.g., "p.Gly12Cys")
	HGVSc              string  // HGVS coding DNA notation (e.g., "c.34G>T")
	GeneType           string  // Gene classification (e.g., "ONCOGENE", "TSG"), empty if not available
	AlphaMissenseScore float64 // AlphaMissense pathogenicity score (0.0-1.0), 0 if not available
	AlphaMissenseClass string  // "likely_benign", "ambiguous", "likely_pathogenic", "" if N/A
}

// GetImpact returns the impact level for a given consequence type.
// For comma-separated consequences, returns the highest impact among all terms.
func GetImpact(consequence string) string {
	best := ImpactModifier
	for rest := consequence; rest != ""; {
		term := rest
		if i := strings.IndexByte(rest, ','); i >= 0 {
			term = rest[:i]
			rest = rest[i+1:]
		} else {
			rest = ""
		}
		var impact string
		switch term {
		case ConsequenceStopGained, ConsequenceFrameshiftVariant,
			ConsequenceStopLost, ConsequenceStartLost,
			ConsequenceSpliceAcceptor, ConsequenceSpliceDonor:
			impact = ImpactHigh
		case ConsequenceMissenseVariant, ConsequenceInframeInsertion,
			ConsequenceInframeDeletion, "inframe_variant":
			impact = ImpactModerate
		case ConsequenceSynonymousVariant, ConsequenceSpliceRegion,
			ConsequenceStopRetained, ConsequenceStartRetained,
			ConsequenceCodingSequenceVariant:
			impact = ImpactLow
		default:
			impact = ImpactModifier
		}
		if ImpactRank(impact) > ImpactRank(best) {
			best = impact
		}
	}
	return best
}

// ImpactRank returns numeric rank for impact comparison (higher = more severe).
func ImpactRank(impact string) int {
	switch impact {
	case ImpactHigh:
		return 3
	case ImpactModerate:
		return 2
	case ImpactLow:
		return 1
	default:
		return 0
	}
}

// FormatVariantID creates a variant identifier from components.
func FormatVariantID(chrom string, pos int64, ref, alt string) string {
	return chrom + "_" + formatInt64(pos) + "_" + ref + "/" + alt
}

func formatInt64(n int64) string {
	if n == 0 {
		return "0"
	}

	var digits []byte
	negative := n < 0
	if negative {
		n = -n
	}

	for n > 0 {
		digits = append([]byte{byte('0' + n%10)}, digits...)
		n /= 10
	}

	if negative {
		digits = append([]byte{'-'}, digits...)
	}

	return string(digits)
}
