// Package annotate provides variant effect prediction functionality.
package annotate

import (
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// ConsequenceResult holds the result of consequence prediction for a variant.
type ConsequenceResult struct {
	Consequence        string
	Impact             string
	CDSPosition        int64
	ProteinPosition    int64
	ProteinEndPosition int64 // End position for multi-AA deletions (0 = single AA)
	RefCodon           string
	AltCodon           string
	RefAA              byte
	AltAA              byte
	EndAA              byte // Amino acid at end position (for multi-AA deletions)
	AminoAcidChange    string
	// AminoAcidsVEP overrides the amino_acids field when VEP writes a bare
	// reference residue rather than a ref/alt pair, which is what it does
	// when it declines to describe the protein. genome-nexus builds its
	// placeholder from this, so p.R700fs needs "R" here; without it the
	// reference residue is unknown and it renders p.*700fs*.
	AminoAcidsVEP      string
	CodonChange        string
	ExonNumber         string
	IntronNumber       string
	CDNAPosition       int64
	HGVSp              string
	HGVSc              string
	FrameshiftStopDist int    // Distance to new stop codon in frameshift (1=at variant pos, 0=unknown)
	StopLostExtDist    int    // Distance to next stop codon in stop-lost (0=unknown)
	InsertedAAs        string // Inserted amino acids for inframe indels (single-letter codes)
	IsDup              bool   // True if inframe insertion is a protein-level duplication
	IsDelIns           bool   // True if insertion also modifies the anchor codon (delins format)
	// ProteinStart/ProteinEnd are the positions VEP reports as protein_start
	// and protein_end. They are the UNSHIFTED position of the change, whereas
	// ProteinPosition is the 3'-shifted position the HGVSp is written at. VEP
	// reports the two differently -- PBRM1 gives protein_start 1169 alongside
	// its own p.I1170Sfs*23 -- so they cannot share one field.
	// SuppressHGVSp marks a variant VEP declines to describe at the protein
	// level, so the final HGVSp is emitted empty and genome-nexus renders its
	// own placeholder from the protein position.
	SuppressHGVSp bool
	ProteinStart  int64
	ProteinEnd    int64
	// HGVSOffset is how far the variant was 3'-shifted to reach its HGVS
	// representation, in bases. VEP reports this as hgvs_offset and omits it
	// when zero; genome-nexus surfaces it as the MAF HGVS_Offset column.
	HGVSOffset int
}

// PredictConsequence determines the effect of a variant on a transcript.
func PredictConsequence(v *vcf.Variant, t *cache.Transcript) *ConsequenceResult {
	result := &ConsequenceResult{}

	// VEP names the exons a variant's genomic span overlaps, whether or not it
	// starts in one, so this is derived from the span rather than from the
	// exon lookup below. It is set before the boundary check because a
	// deletion can start outside the transcript and still reach an exon:
	// RAD51C 17:56769933 is upstream of its own start and VEP reports 1/9.
	result.ExonNumber = spannedExonNumber(v, t)

	// Check if variant is within transcript boundaries. The span decides it,
	// not the start: a deletion can begin before the transcript and still
	// delete part of its first exon, and calling that upstream_gene_variant
	// contradicted the exon this very function had just named for it. VEP
	// reports coding_sequence_variant plus the UTR it reaches.
	spanLast := v.Pos + int64(len(v.Ref)) - 1
	if spanLast < v.Pos {
		spanLast = v.Pos
	}
	if !t.Contains(v.Pos) && !t.Contains(spanLast) && !(v.Pos < t.Start && spanLast > t.End) {
		// Check upstream/downstream
		if v.Pos < t.Start {
			if t.IsForwardStrand() {
				result.Consequence = ConsequenceUpstreamGene
			} else {
				result.Consequence = ConsequenceDownstreamGene
			}
		} else {
			if t.IsForwardStrand() {
				result.Consequence = ConsequenceDownstreamGene
			} else {
				result.Consequence = ConsequenceUpstreamGene
			}
		}
		result.Impact = GetImpact(result.Consequence)
		return result
	}

	// Check if in exon
	exon := t.FindExon(v.Pos)
	if exon == nil && len(v.Ref) == 0 {
		// A pure insertion sits BETWEEN two bases and v.Pos names only one of
		// them; at an exon/intron junction the other may be the coding one. On
		// a reverse-strand transcript v.Pos is the intronic flank, so the
		// variant looked purely intronic and was called splice_donor_variant.
		// VEP treats it as an exonic change plus splice_region_variant --
		// c.990_990+1insT is frameshift_variant / p.R331*, not p.X330_splice --
		// because inserting between the two leaves the essential donor
		// dinucleotide intact. A deletion removes those bases and genuinely
		// does break the site, so this applies to insertions only.
		if other := t.FindExon(v.Pos + 1); other != nil {
			exon = other
		}
	}
	if exon == nil {
		// Intronic - check splice sites (±1-2bp), then splice region (±3-8bp)
		// For indels, check the entire span for splice site overlap
		// Collect the SO terms for every intronic position the variant covers.
		// VEP distinguishes the donor 5th base, the donor region and the
		// polypyrimidine tract, which a single splice_region_variant hides.
		terms := spanSpliceTerms(v, t)
		if insTerms, ok := insertionSpliceTerms(v, t); ok {
			// An insertion covers no base, so the span enumeration would
			// describe whichever of its two flanks v.Pos names.
			terms = insTerms
		}
		if joined := joinConsequenceTerms(terms); joined != "" {
			result.Consequence = joined
			result.Impact = GetImpact(firstTerm(joined))
		} else {
			result.Consequence = ConsequenceIntronVariant
			result.Impact = GetImpact(ConsequenceIntronVariant)
		}
		// For splice donor/acceptor variants on protein-coding transcripts,
		// compute p.X###_splice notation using the nearest exon boundary.
		if t.IsProteinCoding() &&
			(result.Consequence == ConsequenceSpliceDonor || result.Consequence == ConsequenceSpliceAcceptor) {
			pos := int64(0)
			// A deletion that 3'-shifts across the junction is named for the
			// codon it lands on, not the exon boundary: POLD1 c.2959del is
			// p.X987_splice where the boundary codon is 985. The same applies
			// to a deletion that starts in the exon and runs into the intron,
			// which shifts using the stored flank (MET c.3016_3028+7del is
			// p.X1006_splice, not p.X1005_splice).
			if sStart, _, _ := shiftedIntronicDelCDS(v, t); sStart > 0 {
				pos, _ = CDSToCodonPosition(sStart)
			}
			if pos == 0 {
				_, ns := straddlingDeletionShift(v, t)
				if ns == 0 {
					_, ns = acceptorSideDeletionShift(v, t)
				}
				if ns > 0 {
					pos, _ = CDSToCodonPosition(ns)
				}
			}
			if pos == 0 {
				pos = nearestSpliceBoundaryProteinPos(v.Pos, t)
			}
			if pos > 0 {
				result.ProteinPosition = pos
				result.HGVSp = FormatHGVSp(result)
			}
		}
		return result
	}

	// Check if transcript is protein coding
	if !t.IsProteinCoding() {
		if t.Biotype == "miRNA" {
			result.Consequence = ConsequenceMatureMiRNA
		} else {
			result.Consequence = ConsequenceNonCodingExon
		}
		result.Impact = GetImpact(result.Consequence)
		return result
	}

	// Check UTR regions
	inUTR := false
	utrConsequence := ""
	if t.IsForwardStrand() {
		if v.Pos < t.CDSStart {
			inUTR = true
			utrConsequence = Consequence5PrimeUTR
		} else if v.Pos > t.CDSEnd {
			inUTR = true
			utrConsequence = Consequence3PrimeUTR
		}
	} else {
		// Reverse strand: CDSEnd is where the start codon is (higher genomic coord),
		// CDSStart is where the stop codon is (lower genomic coord)
		if v.Pos > t.CDSEnd {
			inUTR = true
			utrConsequence = Consequence5PrimeUTR
		} else if v.Pos < t.CDSStart {
			inUTR = true
			utrConsequence = Consequence3PrimeUTR
		}
	}

	if inUTR {
		// For multi-base variants starting in UTR, check if they span into
		// higher-impact regions (splice sites, start codon, CDS).
		// This applies to both indels (different ref/alt length) and MNVs (same length).
		variantEnd := v.Pos + int64(len(v.Ref)) - 1
		if len(v.Ref) > 1 {
			if v.IsIndel() {
				// Check splice site overlap first (highest priority).
				// Name every region the span reaches, as the coding branch
				// does: a deletion running from the 3'UTR through the
				// acceptor into the intron is
				// splice_acceptor_variant,coding_sequence_variant,
				// 3_prime_UTR_variant,intron_variant, not the bare site.
				if spliceSite := indelSpliceSiteType(v, t); spliceSite != "" {
					result.Consequence = joinConsequenceTerms(
						append(spanSpliceTerms(v, t), spliceSite))
					result.Impact = GetImpact(spliceSite)
					return result
				}
			}
			// Check if variant spans into start codon
			startCodonStart, startCodonEnd := t.CDSStart, t.CDSStart+2
			if !t.IsForwardStrand() {
				startCodonStart, startCodonEnd = t.CDSEnd-2, t.CDSEnd
			}
			if utrConsequence == Consequence5PrimeUTR &&
				variantEnd >= startCodonStart && v.Pos <= startCodonEnd &&
				!startCodonPlaceable(v, t) {
				// VEP needs a CDS coordinate for the variant's own 5' end
				// before it will say the start is lost, and a deletion
				// reaching the ATG from an intron or from outside the
				// transcript has none. It names the regions the span covers
				// instead: NOTCH3 19:15311647-15311872 removes the ATG and is
				// still only coding_sequence_variant,5_prime_UTR_variant. Of
				// the 51 variants in a VEP111 MSK-IMPACT MAF whose 5' end
				// cannot be placed in an exon, not one carries start_lost.
				if joined := joinConsequenceTerms(spanSpliceTerms(v, t)); joined != "" {
					result.Consequence = joined
					result.Impact = GetImpact(firstTerm(joined))
					return result
				}
			}
			if utrConsequence == Consequence5PrimeUTR &&
				variantEnd >= startCodonStart && v.Pos <= startCodonEnd {
				result.Consequence = ConsequenceStartLost
				result.ProteinPosition = 1
				result.Impact = GetImpact(result.Consequence)
				result.HGVSp = FormatHGVSp(result)
				// A variant that begins in the 5'UTR and runs into the start
				// codon has no well-defined protein description; VEP emits none
				// and genome-nexus renders the p.*N* placeholder from the
				// protein position. f60b7a4 did this for variants reaching the
				// start codon from inside the CDS, but one starting in the UTR
				// is caught here first and returned before that point.
				if variantLeavesCDS(v, t, true) {
					result.Consequence = joinConsequenceTerms(
						append(spanSpliceTerms(v, t), ConsequenceStartLost))
					result.Impact = GetImpact(firstTerm(result.Consequence))
					result.HGVSp = ""
				}
				return result
			}
			// Check if variant from 3'UTR spans into the stop codon
			if utrConsequence == Consequence3PrimeUTR {
				var stopCodonStart, stopCodonEnd int64
				if t.IsForwardStrand() {
					stopCodonStart, stopCodonEnd = t.CDSEnd-2, t.CDSEnd
				} else {
					stopCodonStart, stopCodonEnd = t.CDSStart, t.CDSStart+2
				}
				if variantEnd >= stopCodonStart && v.Pos <= stopCodonEnd {
					result.Consequence = ConsequenceStopLost3PrimeUTR
					result.Impact = GetImpact(ConsequenceStopLost)
					return result
				}
			}
		}
		result.Consequence = utrConsequence
		// A UTR base can sit in a splice region like any other exonic base.
		// This branch returned before the check that the coding path makes,
		// so GLI1 c.-26G>C came out as 5_prime_UTR_variant where VEP has
		// splice_region_variant,5_prime_UTR_variant.
		if spanTouchesSpliceRegion(v, t) {
			result.Consequence = joinConsequenceTerms(
				[]string{ConsequenceSpliceRegion, utrConsequence})
		}
		result.Impact = GetImpact(firstTerm(result.Consequence))
		return result
	}

	// Variant is in CDS - calculate coding effect
	result = predictCodingConsequence(v, t, exon, result)

	// For indels, check for higher-impact consequences. The test is whether the
	// variant REMOVES start-codon bases, so it keys on deleted length rather
	// than ref length: a single-base deletion counts (c.1del is start_lost /
	// p.M1? to VEP, but the old len(ref)>1 guard skipped it and left a plain
	// frameshift, p.M1*), while an insertion written VCF-style with an anchor
	// base (ref=G, alt=GG) deletes nothing and must not.
	if v.IsIndel() && len(v.Ref) > len(v.Alt) {
		indelEnd := v.Pos + int64(len(v.Ref)) - 1
		// Check if indel spans the start codon → start_lost
		startCodonStart, startCodonEnd := t.CDSStart, t.CDSStart+2
		if !t.IsForwardStrand() {
			startCodonStart, startCodonEnd = t.CDSEnd-2, t.CDSEnd
		}
		if indelEnd >= startCodonStart && v.Pos <= startCodonEnd &&
			!startCodonPlaceable(v, t) {
			// See startCodonPlaceable: with no CDS coordinate for its own 5'
			// end VEP declines start_lost and names the regions the span
			// covers. TP53 17:7579507-7580166 reaches the ATG from inside
			// intron 1 and is splice_acceptor_variant,splice_donor_variant,
			// splice_donor_5th_base_variant,coding_sequence_variant,
			// 5_prime_UTR_variant,intron_variant.
			if joined := joinConsequenceTerms(spanSpliceTerms(v, t)); joined != "" {
				result.Consequence = joined
				result.Impact = GetImpact(firstTerm(joined))
				result.HGVSp = ""
				return result
			}
		}
		if indelEnd >= startCodonStart && v.Pos <= startCodonEnd {
			result.Consequence = ConsequenceStartLost
			// VEP names the residue the change actually lands on, which is Met1
			// only when the variant sits in the first codon. c.7del is p.V3?,
			// c.4del is p.A2?; the hardcoded 1 reported all of them as p.M1?.
			result.ProteinPosition = 1
			// The position is the one the variant lands on AFTER the 3' shift,
			// not where the MAF places it. ERCC5 deletes c.3, inside the start
			// codon, which is why this branch fires -- but the deletion shifts
			// to c.7 and VEP reports p.V3?, not p.M1?.
			firstDel := v.Pos
			if t.IsReverseStrand() {
				firstDel = v.Pos + int64(len(v.Ref)) - 1
			}
			if cdsPos := GenomicToCDS(firstDel, t); cdsPos >= 1 && len(t.CDSSequence) > 0 {
				delStart := int(cdsPos - 1)
				delEnd := delStart + len(v.Ref) - 1
				if delEnd >= delStart {
					sStart, _ := shiftDeletionThreePrime(delStart, delEnd,
						t.CDSSequence, cdsExonEndIdx(delEnd, t))
					if codon, _ := CDSToCodonPosition(int64(sStart) + 1); codon >= 1 {
						result.ProteinPosition = codon
						if cd := GetCodon(t.CDSSequence, codon); len(cd) == 3 {
							result.RefAA = TranslateCodon(cd)
						}
					}
				}
			}
			result.Impact = GetImpact(result.Consequence)
			result.HGVSp = FormatHGVSp(result)
			// A deletion that runs off the 5' end of the CDS into the UTR has
			// no well-defined protein description and VEP emits none, leaving
			// genome-nexus to synthesize the p.*N* placeholder from the protein
			// position. Emitting our own p.Met1? here would win the argument
			// and lose the comparison.
			if variantLeavesCDS(v, t, true) {
				result.Consequence = joinConsequenceTerms(
					append(spanSpliceTerms(v, t), ConsequenceStartLost))
				result.Impact = GetImpact(firstTerm(result.Consequence))
				result.HGVSp = ""
			}
			if result.Consequence == ConsequenceStartLost {
				// start_lost is a co-term for a deletion, not a verdict on its
				// own: VEP writes frameshift_variant,start_lost or
				// start_lost,inframe_deletion, and never start_lost alone for
				// a DEL - of the 103 deletions carrying it in a VEP111
				// MSK-IMPACT MAF, every one names the length change or the
				// 5'UTR alongside it.
				//
				// This runs after FormatHGVSp above, which keys off the
				// consequence and needs to see the start_lost form to render
				// p.Met1?.
				lengthTerm := ConsequenceInframeDeletion
				if (len(v.Ref)-len(v.Alt))%3 != 0 {
					lengthTerm = ConsequenceFrameshiftVariant
				}
				result.Consequence = joinConsequenceTerms(
					[]string{ConsequenceStartLost, lengthTerm})
				result.Impact = GetImpact(firstTerm(result.Consequence))
			}
			return result
		}
	}

	// A multi-base variant reaching into a splice site is reported as the splice
	// consequence, plus coding_sequence_variant because it also alters coding
	// sequence: VEP gives "splice_acceptor_variant,coding_sequence_variant".
	if spliceSite := indelSpliceSiteType(v, t); spliceSite != "" {
		// Name every region the variant reaches, not only the site it breaks:
		// a deletion running off the end of an exon is also an intron_variant,
		// and may cover the donor 5th base. spliceSite itself is kept because
		// it decides the impact and the two overlap-tests can disagree at the
		// very edge of the span.
		terms := append(spanSpliceTerms(v, t), spliceSite, ConsequenceCodingSequenceVariant)
		// With a CDS coordinate at each end VEP also says what the deletion
		// does to the protein, and that naming replaces coding_sequence_variant
		// (see joinConsequenceTerms). See spliceSpanProteinTerms.
		terms = append(terms, spliceSpanProteinTerms(v, t)...)
		result.Consequence = joinConsequenceTerms(terms)
		result.Impact = GetImpact(spliceSite)
		if t.IsProteinCoding() {
			if pos := nearestSpliceBoundaryProteinPos(v.Pos, t); pos > 0 {
				result.ProteinPosition = pos
				result.ProteinStart = pos
				result.HGVSp = FormatHGVSp(result)
			}
		}
	} else if spanTouchesSpliceRegion(v, t) {
		// Append splice_region_variant if any base the variant covers is near
		// an exon boundary.
		result.Consequence = appendSpliceRegion(result.Consequence)
	}

	return result
}

// predictCodingConsequence calculates the effect on coding sequence.
func predictCodingConsequence(v *vcf.Variant, t *cache.Transcript, exon *cache.Exon, result *ConsequenceResult) *ConsequenceResult {
	// Calculate CDS position
	cdsPos := GenomicToCDS(v.Pos, t)
	if cdsPos < 1 && len(v.Ref) == 0 {
		// Junction insertion (see PredictConsequence): anchor on the coding
		// flank. The +-1 restores the usual convention that the anchor names
		// the base the insertion follows in CODING order, which is the other
		// flank once we cross the junction.
		if p := GenomicToCDS(v.Pos+1, t); p > 0 {
			if t.IsForwardStrand() {
				cdsPos = p - 1
			} else {
				cdsPos = p + 1
			}
		}
	}
	if cdsPos < 1 {
		result.Consequence = ConsequenceIntronVariant
		result.Impact = GetImpact(result.Consequence)
		return result
	}
	result.CDSPosition = cdsPos

	// Calculate codon position
	codonNum, posInCodon := CDSToCodonPosition(cdsPos)
	result.ProteinPosition = codonNum
	result.ProteinStart = codonNum
	if pStart, pEnd := vepProteinSpan(v, t); pStart > 0 {
		result.ProteinStart = pStart
		if pEnd > pStart {
			result.ProteinEnd = pEnd
		}
	}

	// Handle indels
	if v.IsIndel() {
		return predictIndelConsequence(v, t, result)
	}

	// Multi-nucleotide variants (len(ref)==len(alt)>1) may span multiple
	// codons and need special handling.
	if len(v.Ref) > 1 {
		return predictMNVConsequence(v, t, result)
	}

	// SNV - calculate amino acid change
	if len(t.CDSSequence) == 0 || len(v.Alt) == 0 {
		// No CDS sequence or empty alt allele — can't determine AA change
		result.Consequence = ConsequenceCodingSequenceVariant
		result.Impact = GetImpact(result.Consequence)
		return result
	}

	// Get reference codon
	refCodon := GetCodon(t.CDSSequence, codonNum)
	if len(refCodon) != 3 {
		result.Consequence = ConsequenceCodingSequenceVariant
		result.Impact = GetImpact(result.Consequence)
		return result
	}
	result.RefCodon = refCodon

	// Determine the alternate base on the coding strand
	altBase := v.Alt[0]
	if t.IsReverseStrand() {
		altBase = Complement(altBase)
	}

	// Mutate the codon
	altCodon := MutateCodon(refCodon, posInCodon, altBase)
	result.AltCodon = altCodon

	// Translate codons
	result.RefAA = TranslateCodon(refCodon)
	result.AltAA = TranslateCodon(altCodon)

	// Format codon change (lowercase mutated base)
	result.CodonChange = formatCodonChange(refCodon, altCodon, posInCodon)

	// Determine consequence type
	if result.RefAA == result.AltAA {
		if result.RefAA == '*' {
			result.Consequence = ConsequenceStopRetained
		} else {
			result.Consequence = ConsequenceSynonymousVariant
		}
	} else if result.AltAA == '*' {
		result.Consequence = ConsequenceStopGained
		result.AminoAcidChange = formatAAChange(result.RefAA, codonNum, '*')
	} else if result.RefAA == '*' {
		result.Consequence = ConsequenceStopLost
		result.AminoAcidChange = formatAAChange('*', codonNum, result.AltAA)
		// Compute extension length by scanning 3'UTR for next in-frame stop
		result.StopLostExtDist = computeStopLostExtension(t, altCodon)
	} else if result.RefAA == 'M' && codonNum == 1 {
		result.Consequence = ConsequenceStartLost
		result.AminoAcidChange = formatAAChange('M', 1, result.AltAA)
	} else {
		result.Consequence = ConsequenceMissenseVariant
		result.AminoAcidChange = formatAAChange(result.RefAA, codonNum, result.AltAA)
	}

	result.Impact = GetImpact(result.Consequence)
	result.HGVSp = FormatHGVSp(result)
	return result
}

// predictMNVConsequence handles multi-nucleotide variants (len(ref)==len(alt)>1)
// that may affect multiple codons. Single-codon changes are classified as
// missense/stop/synonymous; multi-codon changes are emitted as delins.
func predictMNVConsequence(v *vcf.Variant, t *cache.Transcript, result *ConsequenceResult) *ConsequenceResult {
	if len(t.CDSSequence) == 0 || len(v.Alt) == 0 {
		result.Consequence = ConsequenceCodingSequenceVariant
		result.Impact = GetImpact(result.Consequence)
		return result
	}

	// Codons over the substituted CDS range. MNVs are not 3'-shifted, so the
	// unshifted range is the one VEP reports.
	mnvAlt := v.Alt
	mnvStart := result.CDSPosition
	if t.IsReverseStrand() {
		mnvAlt = ReverseComplement(v.Alt)
		mnvStart = result.CDSPosition - int64(len(v.Ref)) + 1
	}
	result.CodonChange = formatCodonChangeMNV(t.CDSSequence, mnvStart,
		mnvStart+int64(len(v.Ref))-1, mnvAlt)

	startPos, _, deletedAAs, insertedAAs := computeInframeProteinChange(v, t, result.CDSPosition)

	switch {
	case startPos == 0 || (len(deletedAAs) == 0 && len(insertedAAs) == 0):
		// No protein-level change detected: synonymous.
		result.Consequence = ConsequenceSynonymousVariant

	case len(deletedAAs) == 1 && len(insertedAAs) == 1:
		result.ProteinPosition = startPos
		result.RefAA = deletedAAs[0]
		result.AltAA = insertedAAs[0]
		switch {
		case result.AltAA == '*':
			result.Consequence = ConsequenceStopGained
			result.AminoAcidChange = formatAAChange(result.RefAA, startPos, '*')
		case result.RefAA == '*':
			result.Consequence = ConsequenceStopLost
		case result.RefAA == 'M' && startPos == 1:
			result.Consequence = ConsequenceStartLost
		default:
			result.Consequence = ConsequenceMissenseVariant
			result.AminoAcidChange = formatAAChange(result.RefAA, startPos, result.AltAA)
		}

	default:
		// Multi-codon change: emit as delins.
		result.IsDelIns = true
		result.ProteinPosition = startPos
		if len(deletedAAs) > 0 {
			result.RefAA = deletedAAs[0]
		}
		result.InsertedAAs = insertedAAs
		if len(deletedAAs) > 1 {
			result.ProteinEndPosition = startPos + int64(len(deletedAAs)) - 1
			result.EndAA = deletedAAs[len(deletedAAs)-1]
		}

		// Check for stop/start effects in multi-codon MNV. An MNV reaching
		// over several codons can do more than one of these at once, and VEP
		// names each: AR c.3_4delinsTT changes the last base of ATG and makes
		// the next codon a stop, and is reported stop_gained,start_lost. All
		// 6 multi-base substitutions touching a start codon in a VEP111
		// MSK-IMPACT MAF carry start_lost.
		hasInsStop := containsStop(insertedAAs)
		hasDelStop := containsStop(deletedAAs)
		var terms []string
		switch {
		case hasInsStop && !hasDelStop:
			terms = append(terms, ConsequenceStopGained)
		case hasDelStop && !hasInsStop:
			terms = append(terms, ConsequenceStopLost)
		}
		if startPos == 1 && (len(deletedAAs) == 0 || deletedAAs[0] == 'M') {
			terms = append(terms, ConsequenceStartLost)
		}
		if len(terms) == 0 {
			terms = append(terms, ConsequenceMissenseVariant)
		}
		result.Consequence = joinConsequenceTerms(terms)
	}

	result.Impact = GetImpact(result.Consequence)
	result.HGVSp = FormatHGVSp(result)
	return result
}

func predictIndelConsequence(v *vcf.Variant, t *cache.Transcript, result *ConsequenceResult) *ConsequenceResult {
	refLen := len(v.Ref)
	altLen := len(v.Alt)
	diff := altLen - refLen

	// For deletions, compute protein position from the first deleted base
	// instead of the VCF anchor (which is one position before the deletion).
	if refLen > altLen && len(t.CDSSequence) > 0 {
		var firstDelGenomic, lastDelGenomic int64
		if t.IsForwardStrand() {
			firstDelGenomic = v.Pos + 1
			lastDelGenomic = v.Pos + int64(refLen) - 1
		} else {
			firstDelGenomic = v.Pos + int64(refLen) - 1
			lastDelGenomic = v.Pos + 1
		}
		if firstDelCDS := GenomicToCDS(firstDelGenomic, t); firstDelCDS > 0 {
			delCodonNum, _ := CDSToCodonPosition(firstDelCDS)
			result.ProteinPosition = delCodonNum
		}
		// protein_start is reported at the first genuinely changed base. The
		// positions above assume a VCF anchor (v.Pos+1), but genome-nexus sends
		// MAF-style alleles where v.Pos IS the first deleted base, which put
		// protein_start one codon too high on every such deletion.
		shared := 0
		for shared < len(v.Ref) && shared < len(v.Alt) && v.Ref[shared] == v.Alt[shared] {
			shared++
		}
		firstChanged, lastChanged := v.Pos+int64(shared), v.Pos+int64(refLen)-1
		if !t.IsForwardStrand() {
			firstChanged, lastChanged = lastChanged, firstChanged
		}
		if cds := GenomicToCDS(firstChanged, t); cds > 0 {
			n, _ := CDSToCodonPosition(cds)
			result.ProteinStart = n
		}
		if cds := GenomicToCDS(lastChanged, t); cds > 0 {
			n, _ := CDSToCodonPosition(cds)
			if n > result.ProteinStart {
				result.ProteinEnd = n
			}
		}
		// Compute end position for multi-codon deletions
		if lastDelCDS := GenomicToCDS(lastDelGenomic, t); lastDelCDS > 0 {
			endCodonNum, _ := CDSToCodonPosition(lastDelCDS)
			if endCodonNum > result.ProteinPosition {
				result.ProteinEndPosition = endCodonNum
				result.ProteinEnd = endCodonNum
				endCodon := GetCodon(t.CDSSequence, endCodonNum)
				if len(endCodon) == 3 {
					result.EndAA = TranslateCodon(endCodon)
				}
			}
		}
	}

	// Look up the reference amino acid at the (possibly updated) protein position
	if result.ProteinPosition > 0 && len(t.CDSSequence) > 0 {
		refCodon := GetCodon(t.CDSSequence, result.ProteinPosition)
		if len(refCodon) == 3 {
			result.RefAA = TranslateCodon(refCodon)
		}
	}

	// Check if the indel overlaps the stop codon
	stopCodonCDSPos := int64(0)
	if len(t.CDSSequence) >= 3 {
		stopCodonCDSPos = int64(len(t.CDSSequence)) - 2 // 1-based start of last codon
	}

	// Insertions at or after the stop codon: if the stop codon remains intact
	// in the mutant CDS, classify as inframe_insertion,stop_retained_variant.
	if diff > 0 && stopCodonCDSPos > 0 && result.CDSPosition >= stopCodonCDSPos {
		if stopCodonPreserved(v, t, result.CDSPosition) {
			result.Consequence = ConsequenceInframeInsertion + "," + ConsequenceStopRetained
			result.Impact = GetImpact(ConsequenceInframeInsertion)
			result.HGVSp = FormatHGVSp(result)
			return result
		}
	}

	if diff%3 == 0 {
		// In-frame
		if diff > 0 {
			result.Consequence = ConsequenceInframeInsertion
			// Derive the protein change BEFORE reclassifying as stop_gained.
			// The block below is guarded on inframe_insertion, so running the
			// stop check first skipped it and left the delins fields empty --
			// c.4547_4548insTTA rendered as p.K1516* rather than VEP's
			// p.K1516delinsN*.
			// Compute inserted amino acids via protein comparison
			if result.CDSPosition > 0 && result.Consequence == ConsequenceInframeInsertion {
				startPos, endPos, delAAs, insAAs := computeInframeProteinChange(v, t, result.CDSPosition)
				if len(insAAs) > 0 {
					result.InsertedAAs = insAAs
					if len(delAAs) > 0 {
						// Anchor codon(s) modified — emit HGVS delins format.
						result.IsDelIns = true
						result.ProteinPosition = startPos
						refCodon := GetCodon(t.CDSSequence, result.ProteinPosition)
						if len(refCodon) == 3 {
							result.RefAA = TranslateCodon(refCodon)
						}
						if endPos > startPos {
							result.ProteinEndPosition = endPos
							endCodon := GetCodon(t.CDSSequence, result.ProteinEndPosition)
							if len(endCodon) == 3 {
								result.EndAA = TranslateCodon(endCodon)
							}
						}
					} else {
						// Pure insertion: anchor AA unchanged; flanking AAs bracket the insert.
						if startPos > 1 {
							result.ProteinPosition = startPos - 1
							refCodon := GetCodon(t.CDSSequence, result.ProteinPosition)
							if len(refCodon) == 3 {
								result.RefAA = TranslateCodon(refCodon)
							}
							// Downstream flanking AA for p.AA1pos_AA2(pos+1)ins format
							endCodon := GetCodon(t.CDSSequence, result.ProteinPosition+1)
							if len(endCodon) == 3 {
								result.EndAA = TranslateCodon(endCodon)
							}
						}
						// Check for protein-level duplication
						insLen := int64(len(insAAs))
						dupStart := result.ProteinPosition - insLen + 1
						if dupStart >= 1 {
							isDup := true
							for i := int64(0); i < insLen; i++ {
								codon := GetCodon(t.CDSSequence, dupStart+i)
								if len(codon) != 3 || TranslateCodon(codon) != insAAs[i] {
									isDup = false
									break
								}
							}
							if isDup {
								result.IsDup = true
								result.IsDelIns = false
								result.ProteinPosition = dupStart
								result.ProteinEndPosition = dupStart + insLen - 1
								refCodon := GetCodon(t.CDSSequence, result.ProteinPosition)
								if len(refCodon) == 3 {
									result.RefAA = TranslateCodon(refCodon)
								}
								if result.ProteinEndPosition > result.ProteinPosition {
									endCodon := GetCodon(t.CDSSequence, result.ProteinEndPosition)
									if len(endCodon) == 3 {
										result.EndAA = TranslateCodon(endCodon)
									}
								}
							}
						}
					}
				}
			}
			// Now that the protein change is known, reclassify if the
			// insertion introduced a stop. Skipped for a delins: the junction
			// scan misreads those, reporting p.H707* where VEP gives the full
			// p.H707_H712delinsQCF, and VEP puts stop_gained on only 30 of 776
			// in-frame delins variants (3.9%).
			if !isDelIns(v.Ref, v.Alt) && indelCreatesStop(v, t, result.CDSPosition) {
				// The length term stays, as it does for the deletion branch
				// below: VEP writes stop_gained,inframe_insertion. Whether it
				// is dropped again is decided from the codons once those are
				// known, by addStopGainedFromCodons.
				result.Consequence = ConsequenceStopGainedInframeIns
			}
		} else {
			result.Consequence = ConsequenceInframeDeletion
			// Check if in-frame deletion creates a stop codon at the junction.
			// Skipped for a delins, see the insertion branch above.
			if !isDelIns(v.Ref, v.Alt) && indelCreatesStop(v, t, result.CDSPosition) {
				result.Consequence = ConsequenceStopGainedInframeDel
			}
			// Check if in-frame deletion spans the stop codon (stop_lost).
			//
			// The span has to come from both genomic ends. CDSPosition names the
			// LAST deleted base on a reverse-strand transcript, so adding refLen
			// to it walks refLen past the real end and reaches the stop codon
			// from anywhere nearby -- ordinary in-frame deletions were reported
			// as stop_lost with a p.TerNNNXaaext*? protein change, where VEP has
			// a plain deletion (CDK6 c.967_978del is p.L323_A326del).
			if stopCodonCDSPos > 0 && result.CDSPosition > 0 {
				indelEndCDS := cdsSpanEnd(v, t, result.CDSPosition+int64(refLen)-1)
				if indelEndCDS >= stopCodonCDSPos {
					// 3_prime_UTR_variant is earned by leaving the CDS, not
					// by reaching the stop codon: a deletion that takes part
					// of the stop but stops there keeps its length term, and
					// VEP writes stop_lost,inframe_deletion. The split is
					// exact over the 94 in-frame stop-losing deletions in a
					// VEP111 MSK-IMPACT MAF -- the 90 running past the stop
					// carry the UTR term, the 4 ending inside it never do.
					if variantLeavesCDS(v, t, false) {
						result.Consequence = ConsequenceStopLost3PrimeUTR
						result.SuppressHGVSp = true
					} else {
						lengthTerm := ConsequenceInframeDeletion
						if isDelIns(v.Ref, v.Alt) {
							lengthTerm = ConsequenceProteinAltering
						}
						result.Consequence = joinConsequenceTerms(
							[]string{ConsequenceStopLost, lengthTerm})
					}
				}
			}
			// Compute protein-level change (may reveal delins if junction creates new AA)
			if result.CDSPosition > 0 {
				startPos, endPos, _, insAAs := computeInframeProteinChange(v, t, result.CDSPosition)
				if startPos > 0 {
					result.ProteinPosition = startPos
					result.ProteinEndPosition = endPos
					refCodon := GetCodon(t.CDSSequence, startPos)
					if len(refCodon) == 3 {
						result.RefAA = TranslateCodon(refCodon)
					}
					if endPos > startPos {
						endCodon := GetCodon(t.CDSSequence, endPos)
						if len(endCodon) == 3 {
							result.EndAA = TranslateCodon(endCodon)
						}
					}
					if len(insAAs) > 0 {
						// The junction produced new residues, so this is a
						// delins. stop_gained,inframe_deletion renders through
						// the stop_gained branch, which needs the flag to show
						// the whole span (p.L288_K289delins*, not p.L288*).
						result.InsertedAAs = insAAs
						result.IsDelIns = true
					}
				}
			}
		}
	} else {
		// Frameshift
		result.Consequence = ConsequenceFrameshiftVariant
		// An insertion landing between two bases of the start codon destroys
		// it, however far HGVS then shifts the description 3'. VEP decides
		// this on the unshifted position, so CDKN1A 6:36651880^36651881 --
		// written c.12_34dup -- is frameshift_variant,start_lost. Of the 11
		// insertions interior to a start codon in a VEP111 MSK-IMPACT MAF,
		// the 9 frameshifting ones all carry start_lost and the 2 in-frame
		// ones never do: whole codons slot in without displacing the ATG.
		if insertionInsideStartCodon(v, t) {
			result.Consequence = joinConsequenceTerms(
				[]string{ConsequenceFrameshiftVariant, ConsequenceStartLost})
		}
		// Check if frameshift overlaps stop codon (stop_lost)
		if stopCodonCDSPos > 0 && result.CDSPosition > 0 {
			indelEndCDS := cdsSpanEnd(v, t, result.CDSPosition+int64(refLen)-1)
			if indelEndCDS >= stopCodonCDSPos {
				result.Consequence = ConsequenceFrameshiftStopLost
				if variantLeavesCDS(v, t, false) {
					// Past the stop codon and into the 3'UTR there is no
					// well-defined protein product. VEP reports
					// stop_lost,3_prime_UTR_variant and emits no HGVSp,
					// leaving genome-nexus to synthesize the p.*N*
					// placeholder from the protein position.
					result.Consequence = ConsequenceStopLost3PrimeUTR
					result.SuppressHGVSp = true
				}
			}
		}
		// Compute the new amino acid and distance to first stop codon
		if result.CDSPosition > 0 {
			proteinPos, refAA, altAA, stopDist := computeFrameshiftDetails(v, t, result.CDSPosition)
			if proteinPos > 0 {
				result.ProteinPosition = proteinPos
				result.RefAA = refAA
			}
			if altAA != 0 {
				result.AltAA = altAA
			}
			result.FrameshiftStopDist = stopDist
			// A frameshift that immediately hits a stop stays a
			// frameshift_variant. VEP never reports stop_gained on its own for
			// a frameshift-length indel: across 62,967 such variants in a
			// VEP111 MSK-IMPACT MAF it reports frameshift_variant (58,253),
			// frameshift_variant with a splice term, or
			// "stop_gained,frameshift_variant" (517) -- but stop_gained alone
			// zero times. Reclassifying here cost 2,427 rows.
		}
	}

	// A deletion that 3'-shifts onto the essential donor dinucleotide destroys
	// the splice site, so the downstream reading frame is no longer something
	// the transcript predicts. VEP emits no protein description for these and
	// genome-nexus renders p.<pos>fs from the protein position instead.
	if _, _, delLen, k, ok := deletionShiftIntoIntron(v, t); ok && deletionHitsDonor(delLen, k) {
		result.SuppressHGVSp = true
		// The residue must be the one at the position genome-nexus reports,
		// which is ProteinStart. result.RefAA belongs to the frameshift's own
		// first changed codon and can be one further along: PTEN c.801+1del
		// reports position 267 (Lys) while RefAA is the Asp at 268.
		pos := result.ProteinStart
		if pos < 1 {
			pos = result.ProteinPosition
		}
		// A frameshift names the single residue at that position; an in-frame
		// deletion names every residue it removes, which is how VEP arrives at
		// p.DDDEEdel for the five codons of NPM1 c.511_524+1del.
		n := 1
		if delLen%3 == 0 {
			n = delLen / 3
		}
		var aas []byte
		for i := 0; i < n; i++ {
			codon := GetCodon(t.CDSSequence, pos+int64(i))
			if len(codon) != 3 {
				break
			}
			aa := TranslateCodon(codon)
			if aa == 0 {
				break
			}
			aas = append(aas, aa)
		}
		if len(aas) > 0 {
			result.AminoAcidsVEP = string(aas)
		}
	}

	result.Impact = GetImpact(result.Consequence)
	result.HGVSp = FormatHGVSp(result)
	if result.SuppressHGVSp {
		result.HGVSp = ""
	}

	// An in-frame change that both deletes and inserts is reported by VEP as
	// protein_altering_variant rather than a clean inframe insertion/deletion.
	// Relabel last: HGVSp and Impact are derived from the inframe term, and the
	// protein change itself ("p.Lys745_Ala750delinsAsnSer") is unaffected.
	if diff%3 == 0 && isDelIns(v.Ref, v.Alt) {
		result.Consequence = relabelInframeAsProteinAltering(result.Consequence)
	}
	return result
}

// splicedReader provides virtual access to a concatenation of up to 4 string
// parts without allocating the full concatenated string. Used to avoid the
// O(CDS) allocation when scanning mutant sequences for frameshifts and stops.
type splicedReader struct {
	parts   [4]string
	offsets [5]int // cumulative lengths: offsets[i] = sum of len(parts[0..i-1])
}

func newSplicedReader(a, b, c, d string) splicedReader {
	s := splicedReader{parts: [4]string{a, b, c, d}}
	s.offsets[1] = len(a)
	s.offsets[2] = s.offsets[1] + len(b)
	s.offsets[3] = s.offsets[2] + len(c)
	s.offsets[4] = s.offsets[3] + len(d)
	return s
}

func (s *splicedReader) Len() int { return s.offsets[4] }

func (s *splicedReader) At(i int) byte {
	for p := 0; p < 4; p++ {
		if i < s.offsets[p+1] {
			return s.parts[p][i-s.offsets[p]]
		}
	}
	return 0
}

// Codon returns the 3-byte codon at position i. Returns a substring (no alloc)
// when all 3 bytes are in the same part; allocates only at part boundaries.
func (s *splicedReader) Codon(i int) string {
	for p := 0; p < 4; p++ {
		if i < s.offsets[p+1] {
			localIdx := i - s.offsets[p]
			if localIdx+3 <= len(s.parts[p]) {
				return s.parts[p][localIdx : localIdx+3]
			}
			var buf [3]byte
			buf[0] = s.At(i)
			buf[1] = s.At(i + 1)
			buf[2] = s.At(i + 2)
			return string(buf[:])
		}
	}
	return ""
}

// computeFrameshiftDetails builds the mutant CDS for a frameshift variant and
// finds the first amino acid position that actually changes, along with the
// distance to the first new stop codon. If no stop is found in the CDS, it
// continues scanning into the 3'UTR sequence.
// Returns proteinPos=0 if no changed codon is found, stopDist=0 if no stop.
func computeFrameshiftDetails(v *vcf.Variant, t *cache.Transcript, cdsPos int64) (proteinPos int64, refAA byte, altAA byte, stopDist int) {
	if len(t.CDSSequence) == 0 || cdsPos < 1 {
		return 0, 0, 0, 0
	}

	ref := v.Ref
	alt := v.Alt
	if t.IsReverseStrand() {
		ref = ReverseComplement(ref)
		alt = ReverseComplement(alt)
	}

	// Splicing at the raw anchor inserts one base early, which in a repeat run
	// yields a mutant that is NOT the 3'-shifted one and makes the scan below
	// report a first-changed codon several residues too early: BRCA2
	// c.2090_2091dup came out as p.E696Kfs*35 where the protein first differs
	// at 698 (p.L698Nfs*33), residue 696 being Glu in both.
	startIdx := cdsEditIndex(cdsPos, ref, t)

	endIdx := startIdx + len(ref)
	if endIdx > len(t.CDSSequence) {
		endIdx = len(t.CDSSequence)
	}

	// Virtual mutant sequence: CDS prefix + alt + CDS suffix + 3'UTR.
	// Uses splicedReader to avoid allocating the entire concatenated string.
	mut := newSplicedReader(t.CDSSequence[:startIdx], alt, t.CDSSequence[endIdx:], t.UTR3Sequence)

	// Scan from the codon containing the variant start to find the first
	// codon that actually changes amino acid identity in the mutant.
	// This naturally aligns with the HGVS 3' rule: in homopolymer runs, the
	// anchor's codon may appear unchanged (repeat fills the deletion), so the
	// scan advances to the first codon that is truly altered — corresponding
	// to the HGVS-canonicalised (3'-shifted) variant position.
	codonStart := (startIdx / 3) * 3
	for i := codonStart; i+3 <= mut.Len(); i += 3 {
		mutAA := TranslateCodon(mut.Codon(i))
		if proteinPos == 0 {
			// Still looking for the first changed codon.
			var origAA byte
			if i+3 <= len(t.CDSSequence) {
				origAA = TranslateCodon(t.CDSSequence[i : i+3])
			}
			if mutAA != origAA {
				proteinPos = int64(i/3) + 1
				refAA = origAA
				altAA = mutAA
				stopDist = 1
			}
		} else {
			stopDist++
		}
		if proteinPos > 0 && mutAA == '*' {
			return
		}
	}

	// No stop found.
	stopDist = 0
	return
}

// computeStopLostExtension scans the 3'UTR for the next in-frame stop codon
// after a stop-lost variant. The altCodon is the mutated stop codon.
// Returns the number of codons to the new stop (including the new stop codon),
// or 0 if no stop is found within the available 3'UTR sequence.
func computeStopLostExtension(t *cache.Transcript, altCodon string) int {
	if len(t.UTR3Sequence) < 3 {
		return 0
	}

	// The extension starts after the mutated codon. We need to scan the
	// 3'UTR in-frame. The last 3 bases of CDSSequence are the (mutated) stop
	// codon, so the 3'UTR immediately follows in frame.
	dist := 1
	for i := 0; i+3 <= len(t.UTR3Sequence); i += 3 {
		codon := t.UTR3Sequence[i : i+3]
		aa := TranslateCodon(codon)
		if aa == '*' {
			return dist
		}
		dist++
	}

	return 0
}

// computeInframeProteinChange compares original and mutant protein sequences
// for an in-frame indel and returns the protein-level change.
// Returns start/end positions (1-based) of the affected protein range,
// the deleted AAs and inserted AAs (single-letter codes).
func computeInframeProteinChange(v *vcf.Variant, t *cache.Transcript, cdsPos int64) (startPos, endPos int64, deletedAAs, insertedAAs string) {
	if len(t.CDSSequence) == 0 || cdsPos < 1 {
		return
	}

	ref, alt := v.Ref, v.Alt
	if t.IsReverseStrand() {
		ref = ReverseComplement(ref)
		alt = ReverseComplement(alt)
	}

	// For reverse-strand variants, GenomicToCDS(v.Pos) maps the leftmost
	// genomic position to the HIGHEST CDS index (the last ref base on the
	// transcript).  For any multi-base ref we must shift cdsIdx left by
	// len(ref)-1 so that it points to the FIRST ref base.  This applies to
	// insertions, deletions, and complex substitutions uniformly.
	cdsIdx := cdsEditIndex(cdsPos, ref, t)

	refEndIdx := cdsIdx + len(ref)
	if refEndIdx > len(t.CDSSequence) {
		refEndIdx = len(t.CDSSequence)
	}

	codonStart := (cdsIdx / 3) * 3

	// Limit translation to a window around the variant for performance.
	// For in-frame indels, the protein change is local — we only need enough
	// codons to find the first/last differences plus suffix matching.
	indelCodons := (len(ref) + len(alt) + 5) / 3
	// The window must also cover the 3' shift applied below, not just the
	// change itself. At +10 the AR poly-Q deletion c.213_239del had exactly no
	// room left (first=11, L=9, origN=20) and could not shift at all, so it
	// reported p.Q68_Q76del where VEP has p.Q72_Q80del.
	windowCodons := indelCodons + 30
	if windowCodons < 40 {
		windowCodons = 40
	}

	// Translate original window
	origN := (len(t.CDSSequence) - codonStart) / 3
	if origN > windowCodons {
		origN = windowCodons
	}
	if origN <= 0 {
		return
	}

	// Build mutant CDS in the window only
	mutEnd := codonStart + windowCodons*3 + len(alt) - len(ref) + 3
	if refEndIdx > len(t.CDSSequence) {
		refEndIdx = len(t.CDSSequence)
	}
	var mutWindow string
	origEnd := codonStart + windowCodons*3
	if origEnd > len(t.CDSSequence) {
		origEnd = len(t.CDSSequence)
	}
	_ = mutEnd
	// Build local mutant: [codonStart..cdsIdx] + alt + [refEndIdx..origEnd]
	if refEndIdx <= origEnd {
		mutWindow = t.CDSSequence[codonStart:cdsIdx] + alt + t.CDSSequence[refEndIdx:origEnd]
	} else {
		mutWindow = t.CDSSequence[codonStart:cdsIdx] + alt
	}
	mutN := len(mutWindow) / 3
	if mutN <= 0 {
		return
	}

	origAAs := make([]byte, origN)
	mutAAs := make([]byte, mutN)
	for i := range origAAs {
		p := codonStart + i*3
		origAAs[i] = TranslateCodon(t.CDSSequence[p : p+3])
	}
	for i := range mutAAs {
		p := i * 3
		mutAAs[i] = TranslateCodon(mutWindow[p : p+3])
	}

	// Find first position that differs
	first := 0
	minLen := origN
	if mutN < minLen {
		minLen = mutN
	}
	for first < minLen && origAAs[first] == mutAAs[first] {
		first++
	}

	// Find matching suffix from end
	oi := origN - 1
	mi := mutN - 1
	for oi >= first && mi >= first && origAAs[oi] == mutAAs[mi] {
		oi--
		mi--
	}

	// No change found at protein level
	if first > oi && first > mi {
		return
	}

	basePos := int64(codonStart/3) + 1
	startPos = basePos + int64(first)
	if oi >= first {
		endPos = basePos + int64(oi)
		deletedAAs = string(origAAs[first : oi+1])
	}
	if mi >= first {
		insertedAAs = string(mutAAs[first : mi+1])
	}

	// The diff above is 5'-anchored: it reports the leftmost residue that
	// differs. HGVS assigns the most 3' position possible and VEP follows that
	// at the protein level, so inside a repeat the two disagree -- the 9-codon
	// deletion c.213_239del is p.Q72_Q80del to VEP but lands on codon 68 here,
	// with an identical HGVSc either way.
	//
	// A pure deletion slides 3' while the residue leaving the front of the
	// deleted block matches the one entering at the back. Bounded by the
	// translated window, so a shift that would run past it simply stops.
	if len(insertedAAs) == 0 && len(deletedAAs) > 0 {
		L := len(deletedAAs)
		shift := 0
		for first+shift+L < origN && origAAs[first+shift] == origAAs[first+shift+L] {
			shift++
		}
		if shift > 0 {
			startPos = basePos + int64(first+shift)
			endPos = startPos + int64(L) - 1
			deletedAAs = string(origAAs[first+shift : first+shift+L])
		}
	}
	return
}

// stopCodonPreserved checks whether an insertion at/after the stop codon
// preserves the stop codon in the mutant CDS. Uses splicedReader to avoid
// allocating the entire mutant CDS string.
func stopCodonPreserved(v *vcf.Variant, t *cache.Transcript, cdsPos int64) bool {
	if len(t.CDSSequence) < 3 || cdsPos < 1 {
		return false
	}

	ref, alt := v.Ref, v.Alt
	if t.IsReverseStrand() {
		ref = ReverseComplement(ref)
		alt = ReverseComplement(alt)
	}
	cdsIdx := cdsEditIndex(cdsPos, ref, t)

	endIdx := cdsIdx + len(ref)
	if endIdx > len(t.CDSSequence) {
		endIdx = len(t.CDSSequence)
	}

	origStop := len(t.CDSSequence) - 3
	if origStop < 0 {
		return false
	}
	if TranslateCodon(t.CDSSequence[origStop:origStop+3]) != '*' {
		return false
	}

	mut := newSplicedReader(t.CDSSequence[:cdsIdx], alt, t.CDSSequence[endIdx:], "")

	// Check at original position
	if origStop+3 <= mut.Len() && TranslateCodon(mut.Codon(origStop)) == '*' {
		return true
	}

	// Check at shifted position (insertion pushes stop codon right by diff
	// bases). Only an in-frame length change moves the stop codon to a new
	// codon boundary; after a frameshift the same three bases still sit
	// there, but no longer as a codon, and reading them as one finds a stop
	// that is not translated. TP53 c.1180_1181insT splits TGA into T|TGA and
	// is frameshift_variant,stop_lost, not stop_retained_variant.
	diff := len(alt) - len(ref)
	shifted := origStop + diff
	if diff%3 == 0 && shifted >= 0 && shifted+3 <= mut.Len() &&
		TranslateCodon(mut.Codon(shifted)) == '*' {
		return true
	}

	return false
}

// indelCreatesStop checks if an indel creates a stop codon near the variant site.
// Uses splicedReader to avoid allocating the entire mutant CDS string.
func indelCreatesStop(v *vcf.Variant, t *cache.Transcript, cdsPos int64) bool {
	if len(t.CDSSequence) == 0 || cdsPos < 1 {
		return false
	}

	ref := v.Ref
	alt := v.Alt
	if t.IsReverseStrand() {
		ref = ReverseComplement(ref)
		alt = ReverseComplement(alt)
	}

	// Without the anchor corrections the mutant CDS is built from the wrong
	// bases and reports a stop that does not exist: the in-frame deletion
	// c.4299_4301del came out as p.V1434* instead of p.V1434del.
	cdsIdx := cdsEditIndex(cdsPos, ref, t)

	endIdx := cdsIdx + len(ref)
	if endIdx > len(t.CDSSequence) {
		endIdx = len(t.CDSSequence)
	}

	mut := newSplicedReader(t.CDSSequence[:cdsIdx], alt, t.CDSSequence[endIdx:], "")

	// Find the codon-aligned start position for the variant
	codonStart := (cdsIdx / 3) * 3

	// Check codons at and after the variant position for new stop codons.
	nCodons := (len(ref) + 2) / 3
	if nCodons < 2 {
		nCodons = 2
	}
	for i := 0; i < nCodons; i++ {
		pos := codonStart + i*3
		if pos+3 > mut.Len() {
			break
		}
		if TranslateCodon(mut.Codon(pos)) == '*' {
			// Only report if this is a NEW stop. Past the edit the mutant and
			// the original are offset by the length difference, so comparing
			// t.CDSSequence at the same index reads an unrelated codon and a
			// pre-existing stop looks new -- that is how the codon-aligned
			// in-frame deletion c.847_882del came out as p.L283* instead of
			// p.L283_D294del.
			origPos := pos
			if pos >= cdsIdx {
				origPos = pos + len(ref) - len(alt)
			}
			if origPos >= 0 && origPos+3 <= len(t.CDSSequence) {
				if TranslateCodon(t.CDSSequence[origPos:origPos+3]) == '*' {
					return false
				}
			}
			return true
		}
	}
	return false
}

// vepProteinSpan returns the positions VEP reports as protein_start /
// protein_end: the codons holding the lowest and highest CDS coordinate the
// variant touches.
//
// It cannot be derived from v.Pos alone. v.Pos is the lowest *genomic*
// coordinate, which on a reverse-strand transcript maps to the *highest* CDS
// coordinate of the span, so a reverse-strand MNV or insertion would otherwise
// report the codon one past the one VEP names. Both ends are converted and
// ordered in CDS space instead.
//
// A pure insertion replaces nothing and sits between two bases; VEP reports it
// against the flanking pair, so the span is the anchor base and its neighbour.
// Returns 0, 0 when neither end lands in coding sequence.
func vepProteinSpan(v *vcf.Variant, t *cache.Transcript) (start, end int64) {
	lo, hi := v.Pos, v.Pos+1
	if len(v.Ref) > 0 {
		hi = v.Pos + int64(len(v.Ref)) - 1
	}
	a, b := GenomicToCDS(lo, t), GenomicToCDS(hi, t)
	switch {
	case a < 1 && b < 1:
		// An indel reaching out of the exon keeps whichever end is coding.
		return 0, 0
	case a < 1:
		a = b
	case b < 1:
		b = a
	}
	if b < a {
		a, b = b, a
	}
	start, _ = CDSToCodonPosition(a)
	end, _ = CDSToCodonPosition(b)
	return start, end
}

// cdsEditIndex returns the 0-based CDSSequence index at which a variant's edit
// begins, from the CDS position GenomicToCDS produced for v.Pos.
//
// Two corrections are needed, and getting either wrong mis-frames the codon and
// yields a mutant protein the variant does not produce. Both have been fixed
// separately in three call sites before this helper existed, so keep new
// mutant-CDS construction on this function rather than re-deriving it:
//
//   - On a reverse-strand transcript GenomicToCDS maps the leftmost genomic
//     base to the HIGHEST CDS index, so for a multi-base ref the position names
//     the LAST affected base and must be walked back len(ref)-1.
//   - A pure insertion replaces nothing; on the forward strand its bases belong
//     AFTER the anchor base, one position further along.
//
// ref must already be oriented to the coding strand.
func cdsEditIndex(cdsPos int64, ref string, t *cache.Transcript) int {
	idx := int(cdsPos - 1)
	if t.IsReverseStrand() && len(ref) > 1 {
		idx -= len(ref) - 1
	}
	if len(ref) == 0 && t.IsForwardStrand() {
		idx++
	}
	if idx < 0 {
		idx = 0
	}
	if idx > len(t.CDSSequence) {
		idx = len(t.CDSSequence)
	}
	return idx
}

// variantLeavesCDS reports whether the variant reaches outside the CDS on the
// given side, in coding orientation: fivePrime means past the start codon,
// otherwise past the stop codon.
func variantLeavesCDS(v *vcf.Variant, t *cache.Transcript, fivePrime bool) bool {
	lo := v.Pos
	hi := v.Pos + int64(len(v.Ref)) - 1
	if hi < lo {
		hi = lo
	}
	if t.IsForwardStrand() == fivePrime {
		return lo < t.CDSStart
	}
	return hi > t.CDSEnd
}

// cdsSpanEnd returns the highest CDS position the variant touches.
//
// It has to come from both genomic ends. CDSPosition names the LAST affected
// base on a reverse-strand transcript, so adding the ref length to it overshoots
// by that much and reaches whatever lies beyond -- typically the stop codon,
// which is how ordinary in-frame deletions came to be reported as stop_lost.
// Four separate sites derived this by hand before it was factored out; prefer
// this function to re-deriving it.
func cdsSpanEnd(v *vcf.Variant, t *cache.Transcript, fallback int64) int64 {
	lo := GenomicToCDS(v.Pos, t)
	hi := GenomicToCDS(v.Pos+int64(len(v.Ref))-1, t)
	if lo > 0 && hi > 0 {
		if hi > lo {
			return hi
		}
		return lo
	}
	return fallback
}

// GenomicToCDS converts a genomic position to CDS position within a transcript.
// Returns 0 if the position is not in the CDS.
//
// The result carries CDSStartOffset, so on a 5'-incomplete transcript it is the
// position VEP reports -- c.1 being the first base of the notional complete
// first codon -- and remains a valid index into the padded CDSSequence.
func GenomicToCDS(genomicPos int64, t *cache.Transcript) int64 {
	pos := genomicToCDSRaw(genomicPos, t)
	if pos < 1 {
		return 0
	}
	return pos + int64(t.CDSStartOffset)
}

// genomicToCDSRaw counts coding bases from the first annotated CDS base,
// ignoring any missing 5' partial codon.
func genomicToCDSRaw(genomicPos int64, t *cache.Transcript) int64 {
	if !t.IsProteinCoding() || !t.ContainsCDS(genomicPos) {
		return 0
	}

	// Fast path: binary search pre-computed CDS regions.
	if len(t.CDSRegions) > 0 {
		regions := t.CDSRegions
		lo, hi := 0, len(regions)-1
		for lo <= hi {
			mid := lo + (hi-lo)/2
			r := &regions[mid]
			if genomicPos >= r.GenomicStart && genomicPos <= r.GenomicEnd {
				if t.Strand == 1 {
					return r.CDSOffset + (genomicPos - r.GenomicStart) + 1
				}
				return r.CDSOffset + (r.GenomicEnd - genomicPos) + 1
			}
			if genomicPos < r.GenomicStart {
				hi = mid - 1
			} else {
				lo = mid + 1
			}
		}
		return 0
	}

	// Fallback: linear scan (for transcripts without pre-built index).
	var cdsPos int64
	if t.Strand == 1 {
		for _, exon := range t.Exons {
			if !exon.IsCoding() {
				continue
			}
			if genomicPos >= exon.CDSStart && genomicPos <= exon.CDSEnd {
				return cdsPos + genomicPos - exon.CDSStart + 1
			}
			if genomicPos > exon.CDSEnd {
				cdsPos += exon.CDSEnd - exon.CDSStart + 1
			}
		}
	} else {
		for i := len(t.Exons) - 1; i >= 0; i-- {
			exon := t.Exons[i]
			if !exon.IsCoding() {
				continue
			}
			if genomicPos >= exon.CDSStart && genomicPos <= exon.CDSEnd {
				return cdsPos + exon.CDSEnd - genomicPos + 1
			}
			if genomicPos < exon.CDSStart {
				cdsPos += exon.CDSEnd - exon.CDSStart + 1
			}
		}
	}
	return 0
}

// CDSToGenomic converts a CDS position (1-based) to a genomic position within a transcript.
// Returns 0 if the CDS position is out of range.
// This is the reverse of GenomicToCDS.
func CDSToGenomic(cdsPos int64, t *cache.Transcript) int64 {
	if !t.IsProteinCoding() || cdsPos < 1 {
		return 0
	}
	// Undo the 5'-incomplete shift applied by GenomicToCDS. Positions inside the
	// padded first codon have no genomic counterpart.
	cdsPos -= int64(t.CDSStartOffset)
	if cdsPos < 1 {
		return 0
	}

	// Fast path: binary search pre-computed CDS regions.
	if len(t.CDSRegions) > 0 {
		regions := t.CDSRegions
		n := len(regions)
		lo, hi := 0, n-1
		for lo <= hi {
			mid := lo + (hi-lo)/2
			r := &regions[mid]
			regionLen := r.GenomicEnd - r.GenomicStart + 1
			if cdsPos > r.CDSOffset && cdsPos <= r.CDSOffset+regionLen {
				offset := cdsPos - r.CDSOffset - 1
				if t.Strand == 1 {
					return r.GenomicStart + offset
				}
				return r.GenomicEnd - offset
			}
			if t.Strand == 1 {
				// Forward: CDSOffset ascending with array index.
				if cdsPos <= r.CDSOffset {
					hi = mid - 1
				} else {
					lo = mid + 1
				}
			} else {
				// Reverse: CDSOffset descending with array index.
				if cdsPos <= r.CDSOffset {
					lo = mid + 1
				} else {
					hi = mid - 1
				}
			}
		}
		return 0
	}

	// Fallback: linear scan.
	var cumulative int64
	if t.Strand == 1 {
		for _, exon := range t.Exons {
			if !exon.IsCoding() {
				continue
			}
			exonLen := exon.CDSEnd - exon.CDSStart + 1
			if cumulative+exonLen >= cdsPos {
				return exon.CDSStart + (cdsPos - cumulative - 1)
			}
			cumulative += exonLen
		}
	} else {
		for i := len(t.Exons) - 1; i >= 0; i-- {
			exon := t.Exons[i]
			if !exon.IsCoding() {
				continue
			}
			exonLen := exon.CDSEnd - exon.CDSStart + 1
			if cumulative+exonLen >= cdsPos {
				return exon.CDSEnd - (cdsPos - cumulative - 1)
			}
			cumulative += exonLen
		}
	}
	return 0
}

// GenomicToTranscriptPos converts a genomic position to a transcript-relative
// exonic position (1-based). This is used for non-coding transcripts where
// positions are counted from the transcript 5' end. Returns 0 if the position
// is not within an exon.
func GenomicToTranscriptPos(genomicPos int64, t *cache.Transcript) int64 {
	// Fast path: binary search + pre-computed cumulative bases.
	if len(t.ExonCumBases) == len(t.Exons) && len(t.Exons) > 0 {
		idx := t.FindExonIdx(genomicPos)
		if idx < 0 {
			return 0
		}
		if t.Strand == 1 {
			return t.ExonCumBases[idx] + (genomicPos - t.Exons[idx].Start) + 1
		}
		return t.ExonCumBases[idx] + (t.Exons[idx].End - genomicPos) + 1
	}

	// Fallback: linear scan.
	var pos int64
	if t.Strand == 1 {
		for _, exon := range t.Exons {
			if genomicPos >= exon.Start && genomicPos <= exon.End {
				return pos + genomicPos - exon.Start + 1
			}
			if genomicPos > exon.End {
				pos += exon.End - exon.Start + 1
			}
		}
	} else {
		for i := len(t.Exons) - 1; i >= 0; i-- {
			exon := t.Exons[i]
			if genomicPos >= exon.Start && genomicPos <= exon.End {
				return pos + exon.End - genomicPos + 1
			}
			if genomicPos < exon.Start {
				pos += exon.End - exon.Start + 1
			}
		}
	}
	return 0
}

// CDSToCodonPosition converts a CDS position to codon number and position within codon.
// CDS positions are 1-based. Returns codon number (1-based) and position in codon (0, 1, or 2).
func CDSToCodonPosition(cdsPos int64) (codonNumber int64, positionInCodon int) {
	if cdsPos < 1 {
		return 0, 0
	}
	codonNumber = (cdsPos-1)/3 + 1
	positionInCodon = int((cdsPos - 1) % 3)
	return
}

// firstTerm returns the leading (most severe) term of a comma-joined list.
func firstTerm(consequence string) string {
	for i := 0; i < len(consequence); i++ {
		if consequence[i] == ',' {
			return consequence[:i]
		}
	}
	return consequence
}

// spliceSiteType returns the splice site consequence (splice_donor_variant or
// splice_acceptor_variant) if the position is at ±1-2bp on the intron side of
// an exon boundary, or empty string if not at a splice site.
//
// Forward strand: exon.End+1/+2 = donor, exon.Start-1/-2 = acceptor
// Reverse strand: exon.Start-1/-2 = donor, exon.End+1/+2 = acceptor
func spliceSiteType(pos int64, t *cache.Transcript) string {
	idx := t.FindNearestExonIdx(pos)
	if idx < 0 {
		return ""
	}
	// Check the nearest exon and its immediate neighbors (pos is in the intron between exons).
	for _, i := range [3]int{idx - 1, idx, idx + 1} {
		if i < 0 || i >= len(t.Exons) {
			continue
		}
		exon := &t.Exons[i]
		if pos == exon.End+1 || pos == exon.End+2 {
			if t.IsForwardStrand() {
				return ConsequenceSpliceDonor
			}
			return ConsequenceSpliceAcceptor
		}
		if pos == exon.Start-1 || pos == exon.Start-2 {
			if t.IsForwardStrand() {
				return ConsequenceSpliceAcceptor
			}
			return ConsequenceSpliceDonor
		}
	}
	return ""
}

// indelSpliceSiteType checks if an indel's span overlaps a splice site.
// spannedExonNumber renders the Exon_Number VEP reports: the exons a variant's
// genomic span overlaps, as a range when it covers more than one.
//
// The span decides it, not the start position. On a reverse-strand transcript
// the leftmost genomic base of TP53 c.541_559+2del is the intronic one, so a
// start-position lookup finds no exon at all and leaves the field empty where
// VEP writes 5/11. A deletion running from one exon through an intron into the
// next is reported as a range: KMT2D c.5964_6183+90del is 28-29/54.
//
// A junction insertion still comes out empty, which is what VEP does, because
// a pure insertion spans no genomic base and v.Pos is the intronic flank.
func spannedExonNumber(v *vcf.Variant, t *cache.Transcript) string {
	if len(v.Ref) == 0 {
		// A pure insertion covers no base. VEP gives it a zero-length feature
		// whose start is one past its end, so it lies within an exon only
		// when both of those coordinates do - an insertion sitting on a
		// junction is in neither the exon nor the intron, and VEP names no
		// exon for it: CDH1 c.1565+1_1565+2dup and TP53 c.97-1_97insAG are
		// both empty.
		a, b := t.FindExon(v.Pos), t.FindExon(v.Pos+1)
		if a == nil || b == nil || a.Number != b.Number {
			return ""
		}
		return formatExonNumber(a.Number, len(t.Exons))
	}
	spanEnd := v.Pos + int64(len(v.Ref)) - 1
	if spanEnd < v.Pos {
		spanEnd = v.Pos
	}
	first, last := 0, 0
	for i := range t.Exons {
		e := &t.Exons[i]
		if e.End < v.Pos || e.Start > spanEnd {
			continue
		}
		// Exons are numbered in transcript order, so on a reverse-strand
		// transcript the lower number is the rightmost genomic one.
		if first == 0 || e.Number < first {
			first = e.Number
		}
		if e.Number > last {
			last = e.Number
		}
	}
	switch {
	case first == 0:
		return ""
	case first == last:
		return formatExonNumber(first, len(t.Exons))
	default:
		return formatExonRange(first, last, len(t.Exons))
	}
}

// formatExonRange builds the "N-K/M" string VEP uses for a span covering
// several exons.
func formatExonRange(first, last, total int) string {
	var buf [32]byte
	n := putInt64(buf[:], int64(first))
	buf[n] = '-'
	n++
	n += putInt64(buf[n:], int64(last))
	buf[n] = '/'
	n++
	n += putInt64(buf[n:], int64(total))
	return string(buf[:n])
}

// formatExonNumber builds "N/M" string using a stack-allocated buffer.
func formatExonNumber(num, total int) string {
	var buf [16]byte
	n := putInt64(buf[:], int64(num))
	buf[n] = '/'
	n++
	n += putInt64(buf[n:], int64(total))
	return string(buf[:n])
}

// formatAAChange builds an amino acid change string like "G12C" or "*100A"
// using a stack-allocated buffer to avoid multiple string concatenations.
func formatAAChange(refAA byte, pos int64, altAA byte) string {
	var buf [24]byte // refAA + up to 20 digit position + altAA
	n := 0
	buf[n] = refAA
	n++
	n += putInt64(buf[n:], pos)
	buf[n] = altAA
	n++
	return string(buf[:n])
}

// putInt64 writes a positive int64 as decimal digits into buf and returns
// the number of bytes written.
func putInt64(buf []byte, v int64) int {
	if v == 0 {
		buf[0] = '0'
		return 1
	}
	// Write digits in reverse, then reverse.
	var tmp [20]byte
	n := 0
	for v > 0 {
		tmp[n] = byte('0' + v%10)
		v /= 10
		n++
	}
	for i := 0; i < n; i++ {
		buf[i] = tmp[n-1-i]
	}
	return n
}

// spliceRegionCompound maps single consequences to their compound form with
// splice_region_variant appended. Pre-built to avoid runtime concatenation.
// The order within each entry is Ensembl's severity order, the same one
// consequenceRank encodes: splice_region_variant sits between the terms that
// outrank it and those it outranks. Writing it after every term, as this table
// used to, put synonymous_variant and coding_sequence_variant on the wrong
// side - VEP reports splice_region_variant,synonymous_variant.
var spliceRegionCompound = map[string]string{
	ConsequenceMissenseVariant:       ConsequenceMissenseVariant + "," + ConsequenceSpliceRegion,
	ConsequenceStopGained:            ConsequenceStopGained + "," + ConsequenceSpliceRegion,
	ConsequenceFrameshiftVariant:     ConsequenceFrameshiftVariant + "," + ConsequenceSpliceRegion,
	ConsequenceInframeDeletion:       ConsequenceInframeDeletion + "," + ConsequenceSpliceRegion,
	ConsequenceInframeInsertion:      ConsequenceInframeInsertion + "," + ConsequenceSpliceRegion,
	ConsequenceSynonymousVariant:     ConsequenceSpliceRegion + "," + ConsequenceSynonymousVariant,
	ConsequenceStopRetained:          ConsequenceSpliceRegion + "," + ConsequenceStopRetained,
	ConsequenceCodingSequenceVariant: ConsequenceSpliceRegion + "," + ConsequenceCodingSequenceVariant,
}

// appendSpliceRegion appends ",splice_region_variant" to a consequence string.
// Uses a pre-built lookup for common cases to avoid runtime concatenation.
func appendSpliceRegion(consequence string) string {
	// A variant already named by a more specific splice term does not also
	// get splice_region_variant (see joinConsequenceTerms).
	for _, t := range splitConsequenceTerms(consequence) {
		switch t {
		case ConsequenceSpliceDonor, ConsequenceSpliceAcceptor,
			ConsequenceSpliceDonor5thBase, ConsequenceSpliceDonorRegion:
			return consequence
		}
	}
	if compound, ok := spliceRegionCompound[consequence]; ok {
		return compound
	}
	// A consequence the table does not cover is ranked properly rather than
	// having the term stuck on the end: stop_lost,3_prime_UTR_variant takes it
	// in the middle, as stop_lost,splice_region_variant,3_prime_UTR_variant.
	return joinConsequenceTerms(append(splitConsequenceTerms(consequence), ConsequenceSpliceRegion))
}

// For deletions, the affected range is [pos, pos+len(ref)-1]. If any position
// in that range hits a ±1-2bp splice site, returns the splice consequence.
// Returns empty string for SNVs or if no splice site is hit.
func indelSpliceSiteType(v *vcf.Variant, t *cache.Transcript) string {
	// Any variant spanning more than one base can reach a splice site, not just
	// indels: a DNP/TNP/ONP straddling an exon boundary (c.421-1_421delinsAA)
	// hits the acceptor dinucleotide just as an indel would.
	if len(v.Ref) <= 1 {
		return ""
	}

	endPos := v.Pos + int64(len(v.Ref)) - 1
	for pos := v.Pos; pos <= endPos; pos++ {
		if site := spliceSiteType(pos, t); site != "" {
			return site
		}
	}
	return ""
}

// isSpliceRegion checks if a position is within a splice region of any exon.
// Per SO:0001630, splice_region_variant = within 3bp exon side or 3-8bp intron
// side of a splice site. The 1-2bp immediately into the intron are splice
// donor/acceptor territory, not splice region.
func isSpliceRegion(pos int64, t *cache.Transcript) bool {
	idx := t.FindNearestExonIdx(pos)
	if idx < 0 {
		return false
	}
	// Check the nearest exon and its immediate neighbors.
	for _, i := range [3]int{idx - 1, idx, idx + 1} {
		if i < 0 || i >= len(t.Exons) {
			continue
		}
		exon := &t.Exons[i]
		// A splice region needs a splice site, so each boundary only counts
		// when there is an intron on the far side of it. The outermost edges
		// of the transcript have none: VEP reports no splice_region_variant
		// for KMT2D c.2T>C, two bases into the first exon, nor anywhere in
		// single-exon genes like H3C4.
		if i > 0 {
			// Near exon Start boundary
			if pos >= exon.Start && pos <= exon.Start+2 {
				return true
			}
			if pos >= exon.Start-8 && pos <= exon.Start-3 {
				return true
			}
		}
		if i < len(t.Exons)-1 {
			// Near exon End boundary
			if pos >= exon.End-2 && pos <= exon.End {
				return true
			}
			if pos >= exon.End+3 && pos <= exon.End+8 {
				return true
			}
		}
	}
	return false
}

// spanTouchesSpliceRegion reports whether any base a variant covers lies in a
// splice region.
//
// VEP decides this on the whole span, not the start: MDM2 c.880_918del is
// written well inside the exon but its 39 deleted bases reach the boundary,
// and VEP reports frameshift_variant,splice_region_variant where vibe-vep,
// testing only v.Pos, reported the frameshift alone.
//
// A pure insertion covers no base, and VEP gives it a zero-length feature
// whose start is one past its end. Overlapping a region then needs both of
// those coordinates inside it, so an insertion between the third and fourth
// base of an exon is not in the splice region while one between the second
// and third is. Testing either flank on its own over-called it 237 times.
func spanTouchesSpliceRegion(v *vcf.Variant, t *cache.Transcript) bool {
	if len(v.Ref) == 0 {
		if isSpliceRegion(v.Pos, t) && isSpliceRegion(v.Pos+1, t) {
			return true
		}
		// An insertion landing on the junction itself is the exception: its
		// two bases sit either side of the exon edge, so neither of the two
		// windows contains both, yet VEP reports the splice region for
		// RAF1 c.990_990+1insT and 51 others like it. Requiring the exonic
		// flank to be in a splice region keeps the transcript's outer edges
		// out, since those border no intron.
		exonSide, intronSide := v.Pos, v.Pos+1
		if t.FindExon(exonSide) == nil {
			exonSide, intronSide = intronSide, exonSide
		}
		if t.FindExon(exonSide) != nil && t.FindExon(intronSide) == nil {
			return isSpliceRegion(exonSide, t)
		}
		return false
	}
	for pos := v.Pos; pos < v.Pos+int64(len(v.Ref)); pos++ {
		if isSpliceRegion(pos, t) {
			return true
		}
	}
	return false
}

// nearestSpliceBoundaryProteinPos finds the nearest exon boundary to an intronic
// splice-site position and returns the corresponding protein (codon) position.
// Returns 0 if the position cannot be mapped to CDS.
func nearestSpliceBoundaryProteinPos(pos int64, t *cache.Transcript) int64 {
	if !t.IsProteinCoding() {
		return 0
	}
	// Find the closest exon boundary on the coding side of the splice site.
	var boundaryGenomic int64
	minDist := int64(1<<62 - 1)
	idx := t.FindNearestExonIdx(pos)
	if idx < 0 {
		return 0
	}
	for _, i := range [3]int{idx - 1, idx, idx + 1} {
		if i < 0 || i >= len(t.Exons) {
			continue
		}
		exon := &t.Exons[i]
		if !exon.IsCoding() {
			continue
		}
		if d := abs64(pos - exon.End); d < minDist && d <= 2 {
			minDist = d
			boundaryGenomic = exon.CDSEnd
		}
		if d := abs64(pos - exon.Start); d < minDist && d <= 2 {
			minDist = d
			boundaryGenomic = exon.CDSStart
		}
	}
	if boundaryGenomic == 0 {
		return 0
	}
	cdsPos := GenomicToCDS(boundaryGenomic, t)
	if cdsPos < 1 {
		return 0
	}
	codonNum, _ := CDSToCodonPosition(cdsPos)
	return codonNum
}

// containsStop returns true if the amino acid string contains a stop codon ('*').
func containsStop(aas string) bool {
	for i := 0; i < len(aas); i++ {
		if aas[i] == '*' {
			return true
		}
	}
	return false
}

func abs64(x int64) int64 {
	if x < 0 {
		return -x
	}
	return x
}

// formatCodonChange formats the codon change the way Ensembl VEP does: both
// codons are lowercased except the changed base, which is uppercased in the
// reference codon as well as the alternate one (e.g. "gGt/gTt").
// Uses byte arithmetic for case conversion to avoid allocations.
func formatCodonChange(refCodon, altCodon string, posInCodon int) string {
	var buf [7]byte // 3 ref + '/' + 3 alt
	for i := 0; i < 3; i++ {
		if i == posInCodon {
			buf[i] = refCodon[i] &^ 0x20   // uppercase mutated ref
			buf[4+i] = altCodon[i] &^ 0x20 // uppercase mutated alt
		} else {
			buf[i] = refCodon[i] | 0x20   // lowercase unchanged ref
			buf[4+i] = altCodon[i] | 0x20 // lowercase unchanged alt
		}
	}
	buf[3] = '/'
	return string(buf[:])
}

// insertionInsideStartCodon reports whether a pure insertion sits strictly
// between two bases of the transcript's start codon. The insertion occupies no
// base of its own: it falls between v.Pos and v.Pos+1, so both flanks must be
// start-codon bases for the ATG to be broken apart.
func insertionInsideStartCodon(v *vcf.Variant, t *cache.Transcript) bool {
	if len(v.Ref) != 0 || len(v.Alt) == 0 || !t.IsProteinCoding() {
		return false
	}
	start, end := t.CDSStart, t.CDSStart+2
	if !t.IsForwardStrand() {
		start, end = t.CDSEnd-2, t.CDSEnd
	}
	return v.Pos >= start && v.Pos+1 <= end
}

// startCodonPlaceable reports whether the variant's 5'-most base, in transcript
// orientation, lies in an exon of the transcript. VEP derives translation_start
// from that base; an intronic or out-of-transcript 5' end leaves it undefined,
// and every consequence that needs it -- start_lost above all -- is then
// unavailable however much coding sequence the variant goes on to remove.
func startCodonPlaceable(v *vcf.Variant, t *cache.Transcript) bool {
	five := v.Pos
	if !t.IsForwardStrand() && len(v.Ref) > 0 {
		five = v.Pos + int64(len(v.Ref)) - 1
	}
	return t.FindExon(five) != nil
}
