// Package annotate provides variant effect prediction functionality.
package annotate

import (
	"strconv"

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
	CodonChange        string
	ExonNumber         string
	IntronNumber       string
	CDNAPosition       int64
	HGVSp              string
	HGVSc              string
	FrameshiftStopDist int // Distance to new stop codon in frameshift (1=at variant pos, 0=unknown)
	StopLostExtDist    int // Distance to next stop codon in stop-lost (0=unknown)
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
		// For indels, check the entire span for splice site overlap
		if spliceSite := indelSpliceSiteType(v, t); spliceSite != "" {
			result.Consequence = spliceSite
			result.Impact = GetImpact(spliceSite)
		} else if spliceSite := spliceSiteType(v.Pos, t); spliceSite != "" {
			result.Consequence = spliceSite
			result.Impact = GetImpact(spliceSite)
		} else if isSpliceRegion(v.Pos, t) {
			result.Consequence = ConsequenceSpliceRegion + "," + ConsequenceIntronVariant
			result.Impact = GetImpact(ConsequenceSpliceRegion)
		} else {
			result.Consequence = ConsequenceIntronVariant
			result.Impact = GetImpact(ConsequenceIntronVariant)
		}
		// For splice donor/acceptor variants on protein-coding transcripts,
		// compute p.X###_splice notation using the nearest exon boundary.
		if t.IsProteinCoding() &&
			(result.Consequence == ConsequenceSpliceDonor || result.Consequence == ConsequenceSpliceAcceptor) {
			if pos := nearestSpliceBoundaryProteinPos(v.Pos, t); pos > 0 {
				result.ProteinPosition = pos
				result.HGVSp = FormatHGVSp(result)
			}
		}
		return result
	}

	// Set exon number
	result.ExonNumber = strconv.Itoa(exon.Number) + "/" + strconv.Itoa(len(t.Exons))

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
		// For large indels starting in UTR, check if deletion spans into
		// higher-impact regions (splice sites, start codon, CDS)
		if v.IsIndel() && len(v.Ref) > 1 {
			// Check splice site overlap first (highest priority)
			if spliceSite := indelSpliceSiteType(v, t); spliceSite != "" {
				result.Consequence = spliceSite
				result.Impact = GetImpact(spliceSite)
				return result
			}
			// Check if deletion spans into start codon
			indelEnd := v.Pos + int64(len(v.Ref)) - 1
			startCodonStart, startCodonEnd := t.CDSStart, t.CDSStart+2
			if !t.IsForwardStrand() {
				startCodonStart, startCodonEnd = t.CDSEnd-2, t.CDSEnd
			}
			if indelEnd >= startCodonStart && v.Pos <= startCodonEnd {
				result.Consequence = ConsequenceStartLost
				result.Impact = GetImpact(result.Consequence)
				return result
			}
		}
		result.Consequence = utrConsequence
		result.Impact = GetImpact(result.Consequence)
		return result
	}

	// Variant is in CDS - calculate coding effect
	result = predictCodingConsequence(v, t, exon, result)

	// For indels, check for higher-impact consequences
	if v.IsIndel() && len(v.Ref) > 1 {
		indelEnd := v.Pos + int64(len(v.Ref)) - 1
		// Check if indel spans the start codon → start_lost
		startCodonStart, startCodonEnd := t.CDSStart, t.CDSStart+2
		if !t.IsForwardStrand() {
			startCodonStart, startCodonEnd = t.CDSEnd-2, t.CDSEnd
		}
		if indelEnd >= startCodonStart && v.Pos <= startCodonEnd {
			result.Consequence = ConsequenceStartLost
			result.Impact = GetImpact(result.Consequence)
			return result
		}
	}

	// For indels spanning into a splice site, upgrade to splice donor/acceptor
	if spliceSite := indelSpliceSiteType(v, t); spliceSite != "" {
		result.Consequence = spliceSite
		result.Impact = GetImpact(spliceSite)
	} else if isSpliceRegion(v.Pos, t) {
		// Append splice_region_variant if near exon boundary
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
		if result.RefAA == '*' {
			result.Consequence = ConsequenceStopRetained
		} else {
			result.Consequence = ConsequenceSynonymousVariant
		}
	} else if result.AltAA == '*' {
		result.Consequence = ConsequenceStopGained
		result.AminoAcidChange = string(result.RefAA) + strconv.FormatInt(codonNum, 10) + "*"
	} else if result.RefAA == '*' {
		result.Consequence = ConsequenceStopLost
		result.AminoAcidChange = "*" + strconv.FormatInt(codonNum, 10) + string(result.AltAA)
		// Compute extension length by scanning 3'UTR for next in-frame stop
		result.StopLostExtDist = computeStopLostExtension(t, altCodon)
	} else if result.RefAA == 'M' && codonNum == 1 {
		result.Consequence = ConsequenceStartLost
		result.AminoAcidChange = "M1" + string(result.AltAA)
	} else {
		result.Consequence = ConsequenceMissenseVariant
		result.AminoAcidChange = string(result.RefAA) + strconv.FormatInt(codonNum, 10) + string(result.AltAA)
	}

	result.Impact = GetImpact(result.Consequence)
	result.HGVSp = FormatHGVSp(result)
	return result
}

// predictIndelConsequence determines the effect of an insertion or deletion.
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
		// Compute end position for multi-codon deletions
		if lastDelCDS := GenomicToCDS(lastDelGenomic, t); lastDelCDS > 0 {
			endCodonNum, _ := CDSToCodonPosition(lastDelCDS)
			if endCodonNum > result.ProteinPosition {
				result.ProteinEndPosition = endCodonNum
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

	if diff%3 == 0 {
		// In-frame
		if diff > 0 {
			result.Consequence = ConsequenceInframeInsertion
			// Check if in-frame insertion creates a stop codon
			if indelCreatesStop(v, t, result.CDSPosition) {
				result.Consequence = ConsequenceStopGained
			}
		} else {
			result.Consequence = ConsequenceInframeDeletion
		}
	} else {
		// Frameshift
		result.Consequence = ConsequenceFrameshiftVariant
		// Check if frameshift overlaps stop codon (stop_lost)
		if stopCodonCDSPos > 0 && result.CDSPosition > 0 {
			indelEndCDS := result.CDSPosition + int64(refLen) - 1
			if indelEndCDS >= stopCodonCDSPos {
				result.Consequence = ConsequenceFrameshiftVariant + "," + ConsequenceStopLost
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
		}
	}

	result.Impact = GetImpact(result.Consequence)
	result.HGVSp = FormatHGVSp(result)
	return result
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

	cdsIdx := int(cdsPos - 1) // 0-based index in CDS

	ref := v.Ref
	alt := v.Alt
	if t.IsReverseStrand() {
		ref = ReverseComplement(ref)
		alt = ReverseComplement(alt)
	}

	endIdx := cdsIdx + len(ref)
	if endIdx > len(t.CDSSequence) {
		endIdx = len(t.CDSSequence)
	}

	// Build mutant sequence: CDS with the variant applied, plus 3'UTR for scanning
	mutSeq := t.CDSSequence[:cdsIdx] + alt + t.CDSSequence[endIdx:] + t.UTR3Sequence

	// Scan from the codon containing the variant position forward,
	// comparing each mutant codon against the original to find the first change.
	codonStart := (cdsIdx / 3) * 3
	for i := codonStart; i+3 <= len(mutSeq); i += 3 {
		mutAA := TranslateCodon(mutSeq[i : i+3])
		if proteinPos == 0 {
			// Still looking for the first changed codon
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

	// No stop found
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

// indelCreatesStop checks if an indel creates a stop codon near the variant site.
// It reconstructs the local CDS around the indel and checks the first few codons.
func indelCreatesStop(v *vcf.Variant, t *cache.Transcript, cdsPos int64) bool {
	if len(t.CDSSequence) == 0 || cdsPos < 1 {
		return false
	}

	// Build the mutant CDS by splicing the alt allele into the CDS
	// VCF indels share a prefix base: REF=A ALT=ATAC means insert TAC after pos
	// For deletion: REF=ATAC ALT=A means delete TAC after pos
	cdsIdx := int(cdsPos - 1) // 0-based index in CDS

	// Determine the ref/alt on the coding strand
	ref := v.Ref
	alt := v.Alt
	if t.IsReverseStrand() {
		ref = ReverseComplement(ref)
		alt = ReverseComplement(alt)
	}

	// Reconstruct local mutant sequence: take CDS before variant, insert alt, then CDS after ref
	endIdx := cdsIdx + len(ref)
	if endIdx > len(t.CDSSequence) {
		endIdx = len(t.CDSSequence)
	}

	mutCDS := t.CDSSequence[:cdsIdx] + alt + t.CDSSequence[endIdx:]

	// Find the codon-aligned start position for the variant
	codonStart := (cdsIdx / 3) * 3

	// Check only the codon containing the variant for a new stop codon.
	// This catches cases like in-frame insertions that create a stop codon
	// at the exact variant position (e.g., TAC → TAG|TAC).
	if codonStart+3 <= len(mutCDS) {
		codon := mutCDS[codonStart : codonStart+3]
		if TranslateCodon(codon) == '*' {
			// Only report if this is a NEW stop (not present in original)
			if codonStart+3 <= len(t.CDSSequence) {
				origCodon := t.CDSSequence[codonStart : codonStart+3]
				if TranslateCodon(origCodon) == '*' {
					return false
				}
			}
			return true
		}
	}
	return false
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

// GenomicToTranscriptPos converts a genomic position to a transcript-relative
// exonic position (1-based). This is used for non-coding transcripts where
// positions are counted from the transcript 5' end. Returns 0 if the position
// is not within an exon.
func GenomicToTranscriptPos(genomicPos int64, t *cache.Transcript) int64 {
	var pos int64

	if t.IsForwardStrand() {
		for _, exon := range t.Exons {
			if genomicPos >= exon.Start && genomicPos <= exon.End {
				pos += genomicPos - exon.Start + 1
				return pos
			}
			if genomicPos > exon.End {
				pos += exon.End - exon.Start + 1
			}
		}
	} else {
		for i := len(t.Exons) - 1; i >= 0; i-- {
			exon := t.Exons[i]
			if genomicPos >= exon.Start && genomicPos <= exon.End {
				pos += exon.End - genomicPos + 1
				return pos
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

// indelSpliceSiteType checks if an indel's span overlaps a splice site.
// For deletions, the affected range is [pos, pos+len(ref)-1]. If any position
// in that range hits a ±1-2bp splice site, returns the splice consequence.
// Returns empty string for SNVs or if no splice site is hit.
func indelSpliceSiteType(v *vcf.Variant, t *cache.Transcript) string {
	if !v.IsIndel() || len(v.Ref) <= 1 {
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
	for _, exon := range t.Exons {
		if !exon.IsCoding() {
			continue
		}
		// Check exon.End boundary (splice site at End+1, End+2)
		if d := abs64(pos - exon.End); d < minDist && d <= 2 {
			minDist = d
			boundaryGenomic = exon.CDSEnd
		}
		// Check exon.Start boundary (splice site at Start-1, Start-2)
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

func abs64(x int64) int64 {
	if x < 0 {
		return -x
	}
	return x
}

// formatCodonChange formats the codon change string with lowercase mutated base.
// Uses byte arithmetic for case conversion to avoid allocations.
func formatCodonChange(refCodon, altCodon string, posInCodon int) string {
	var buf [7]byte // 3 ref + '/' + 3 alt
	for i := 0; i < 3; i++ {
		buf[i] = refCodon[i] | 0x20 // lowercase all ref
		if i == posInCodon {
			buf[4+i] = altCodon[i] &^ 0x20 // uppercase mutated alt
		} else {
			buf[4+i] = altCodon[i] | 0x20 // lowercase unchanged alt
		}
	}
	buf[3] = '/'
	return string(buf[:])
}
