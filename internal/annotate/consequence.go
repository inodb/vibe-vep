// Package annotate provides variant effect prediction functionality.
package annotate

import (
	"fmt"
	"strings"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// ConsequenceResult holds the result of consequence prediction for a variant.
type ConsequenceResult struct {
	Consequence     string
	Impact          string
	CDSPosition     int64
	ProteinPosition int64
	RefCodon        string
	AltCodon        string
	RefAA           byte
	AltAA           byte
	AminoAcidChange string
	CodonChange     string
	ExonNumber      string
	IntronNumber    string
	CDNAPosition    int64
}

// PredictConsequence determines the effect of a variant on a transcript.
func PredictConsequence(v *vcf.Variant, t *cache.Transcript) *ConsequenceResult {
	result := &ConsequenceResult{}

	// Check if variant is within transcript boundaries
	if !t.Contains(v.Pos) {
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
	if exon == nil {
		// Intronic - check splice sites (±1-2bp), then splice region (±3-8bp)
		if spliceSite := spliceSiteType(v.Pos, t); spliceSite != "" {
			result.Consequence = spliceSite
			result.Impact = GetImpact(spliceSite)
		} else if isSpliceRegion(v.Pos, t) {
			result.Consequence = ConsequenceSpliceRegion + "," + ConsequenceIntronVariant
			result.Impact = GetImpact(ConsequenceSpliceRegion)
		} else {
			result.Consequence = ConsequenceIntronVariant
			result.Impact = GetImpact(ConsequenceIntronVariant)
		}
		return result
	}

	// Set exon number
	result.ExonNumber = fmt.Sprintf("%d/%d", exon.Number, len(t.Exons))

	// Check if transcript is protein coding
	if !t.IsProteinCoding() {
		result.Consequence = ConsequenceNonCodingExon
		result.Impact = GetImpact(result.Consequence)
		return result
	}

	// Check UTR regions
	if t.IsForwardStrand() {
		if v.Pos < t.CDSStart {
			result.Consequence = Consequence5PrimeUTR
			result.Impact = GetImpact(result.Consequence)
			return result
		}
		if v.Pos > t.CDSEnd {
			result.Consequence = Consequence3PrimeUTR
			result.Impact = GetImpact(result.Consequence)
			return result
		}
	} else {
		// Reverse strand: CDSEnd is where the start codon is (higher genomic coord),
		// CDSStart is where the stop codon is (lower genomic coord)
		if v.Pos > t.CDSEnd {
			result.Consequence = Consequence5PrimeUTR
			result.Impact = GetImpact(result.Consequence)
			return result
		}
		if v.Pos < t.CDSStart {
			result.Consequence = Consequence3PrimeUTR
			result.Impact = GetImpact(result.Consequence)
			return result
		}
	}

	// Variant is in CDS - calculate coding effect
	result = predictCodingConsequence(v, t, exon, result)

	// Append splice_region_variant if near exon boundary
	if isSpliceRegion(v.Pos, t) {
		result.Consequence += "," + ConsequenceSpliceRegion
	}

	return result
}

// predictCodingConsequence calculates the effect on coding sequence.
func predictCodingConsequence(v *vcf.Variant, t *cache.Transcript, exon *cache.Exon, result *ConsequenceResult) *ConsequenceResult {
	// Calculate CDS position
	cdsPos := GenomicToCDS(v.Pos, t)
	if cdsPos < 1 {
		result.Consequence = ConsequenceIntronVariant
		result.Impact = GetImpact(result.Consequence)
		return result
	}
	result.CDSPosition = cdsPos

	// Calculate codon position
	codonNum, posInCodon := CDSToCodonPosition(cdsPos)
	result.ProteinPosition = codonNum

	// Handle indels
	if v.IsIndel() {
		return predictIndelConsequence(v, t, result)
	}

	// SNV - calculate amino acid change
	if len(t.CDSSequence) == 0 {
		// No CDS sequence available, can't determine AA change
		result.Consequence = ConsequenceMissenseVariant // Assume missense
		result.Impact = GetImpact(result.Consequence)
		return result
	}

	// Get reference codon
	refCodon := GetCodon(t.CDSSequence, codonNum)
	if len(refCodon) != 3 {
		result.Consequence = ConsequenceMissenseVariant
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
		result.Consequence = ConsequenceSynonymousVariant
	} else if result.AltAA == '*' {
		result.Consequence = ConsequenceStopGained
		result.AminoAcidChange = fmt.Sprintf("%c%d*", result.RefAA, codonNum)
	} else if result.RefAA == '*' {
		result.Consequence = ConsequenceStopLost
		result.AminoAcidChange = fmt.Sprintf("*%d%c", codonNum, result.AltAA)
	} else if result.RefAA == 'M' && codonNum == 1 {
		result.Consequence = ConsequenceStartLost
		result.AminoAcidChange = fmt.Sprintf("M1%c", result.AltAA)
	} else {
		result.Consequence = ConsequenceMissenseVariant
		result.AminoAcidChange = fmt.Sprintf("%c%d%c", result.RefAA, codonNum, result.AltAA)
	}

	result.Impact = GetImpact(result.Consequence)
	return result
}

// predictIndelConsequence determines the effect of an insertion or deletion.
func predictIndelConsequence(v *vcf.Variant, t *cache.Transcript, result *ConsequenceResult) *ConsequenceResult {
	refLen := len(v.Ref)
	altLen := len(v.Alt)
	diff := altLen - refLen

	if diff%3 == 0 {
		// In-frame
		if diff > 0 {
			result.Consequence = ConsequenceInframeInsertion
		} else {
			result.Consequence = ConsequenceInframeDeletion
		}
	} else {
		// Frameshift
		result.Consequence = ConsequenceFrameshiftVariant
	}

	result.Impact = GetImpact(result.Consequence)
	return result
}

// GenomicToCDS converts a genomic position to CDS position within a transcript.
// Returns 0 if the position is not in the CDS.
func GenomicToCDS(genomicPos int64, t *cache.Transcript) int64 {
	if !t.IsProteinCoding() {
		return 0
	}

	if !t.ContainsCDS(genomicPos) {
		return 0
	}

	var cdsPos int64 = 0

	if t.IsForwardStrand() {
		for _, exon := range t.Exons {
			// Skip non-coding exons
			if !exon.IsCoding() {
				continue
			}

			// Determine the CDS portion of this exon
			cdsStart := exon.CDSStart
			cdsEnd := exon.CDSEnd

			if genomicPos >= cdsStart && genomicPos <= cdsEnd {
				// Position is in this exon's CDS
				cdsPos += genomicPos - cdsStart + 1
				return cdsPos
			}

			if genomicPos > cdsEnd {
				// Position is after this exon, add full exon length
				cdsPos += cdsEnd - cdsStart + 1
			}
		}
	} else {
		// Reverse strand - iterate exons in reverse order
		for i := len(t.Exons) - 1; i >= 0; i-- {
			exon := t.Exons[i]
			if !exon.IsCoding() {
				continue
			}

			cdsStart := exon.CDSStart
			cdsEnd := exon.CDSEnd

			if genomicPos >= cdsStart && genomicPos <= cdsEnd {
				// Position is in this exon's CDS (count from end for reverse strand)
				cdsPos += cdsEnd - genomicPos + 1
				return cdsPos
			}

			if genomicPos < cdsStart {
				// Position is after this exon (in genomic terms, before in CDS terms)
				cdsPos += cdsEnd - cdsStart + 1
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

// spliceSiteType returns the splice site consequence (splice_donor_variant or
// splice_acceptor_variant) if the position is at ±1-2bp on the intron side of
// an exon boundary, or empty string if not at a splice site.
//
// Forward strand: exon.End+1/+2 = donor, exon.Start-1/-2 = acceptor
// Reverse strand: exon.Start-1/-2 = donor, exon.End+1/+2 = acceptor
func spliceSiteType(pos int64, t *cache.Transcript) string {
	for _, exon := range t.Exons {
		// Positions after exon.End (intron side): +1, +2
		if pos == exon.End+1 || pos == exon.End+2 {
			if t.IsForwardStrand() {
				return ConsequenceSpliceDonor
			}
			return ConsequenceSpliceAcceptor
		}
		// Positions before exon.Start (intron side): -1, -2
		if pos == exon.Start-1 || pos == exon.Start-2 {
			if t.IsForwardStrand() {
				return ConsequenceSpliceAcceptor
			}
			return ConsequenceSpliceDonor
		}
	}
	return ""
}

// isSpliceRegion checks if a position is within a splice region of any exon.
// Per SO:0001630, splice_region_variant = within 3bp exon side or 3-8bp intron
// side of a splice site. The 1-2bp immediately into the intron are splice
// donor/acceptor territory, not splice region.
func isSpliceRegion(pos int64, t *cache.Transcript) bool {
	for _, exon := range t.Exons {
		// Near exon Start boundary
		// Exon side: exon.Start, exon.Start+1, exon.Start+2
		if pos >= exon.Start && pos <= exon.Start+2 {
			return true
		}
		// Intron side: 3-8bp before exon start (skip ±1,±2 = splice donor/acceptor)
		if pos >= exon.Start-8 && pos <= exon.Start-3 {
			return true
		}

		// Near exon End boundary
		// Exon side: exon.End-2, exon.End-1, exon.End
		if pos >= exon.End-2 && pos <= exon.End {
			return true
		}
		// Intron side: 3-8bp after exon end
		if pos >= exon.End+3 && pos <= exon.End+8 {
			return true
		}
	}
	return false
}

// formatCodonChange formats the codon change string with lowercase mutated base.
func formatCodonChange(refCodon, altCodon string, posInCodon int) string {
	// Format: lowercase for unchanged bases, uppercase for mutated base in alt
	refBytes := []byte(refCodon)
	altBytes := []byte(altCodon)

	for i := 0; i < 3; i++ {
		if i == posInCodon {
			refBytes[i] = byte(strings.ToLower(string(refBytes[i]))[0])
			altBytes[i] = byte(strings.ToUpper(string(altBytes[i]))[0])
		} else {
			refBytes[i] = byte(strings.ToLower(string(refBytes[i]))[0])
			altBytes[i] = byte(strings.ToLower(string(altBytes[i]))[0])
		}
	}

	return string(refBytes) + "/" + string(altBytes)
}
