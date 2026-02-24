package annotate

import (
	"fmt"
	"strings"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// FormatHGVSc formats the HGVS coding DNA notation for a variant on a transcript.
// Returns empty string for non-coding transcripts or upstream/downstream variants.
func FormatHGVSc(v *vcf.Variant, t *cache.Transcript, result *ConsequenceResult) string {
	if t == nil || v == nil || result == nil {
		return ""
	}

	// Skip upstream/downstream/intergenic
	switch result.Consequence {
	case ConsequenceUpstreamGene, ConsequenceDownstreamGene, ConsequenceIntergenicVariant:
		return ""
	}

	// Non-coding transcripts use n. notation — skip for now
	if !t.IsProteinCoding() {
		return ""
	}

	// Get ref/alt on coding strand
	ref := v.Ref
	alt := v.Alt
	if t.IsReverseStrand() {
		ref = ReverseComplement(ref)
		alt = ReverseComplement(alt)
	}

	// Determine position string for the variant start
	startPosStr := genomicToHGVScPos(v.Pos, t)
	if startPosStr == "" {
		return ""
	}

	// SNV
	if !v.IsIndel() {
		return fmt.Sprintf("c.%s%s>%s", startPosStr, ref, alt)
	}

	// Indel handling
	// VCF indels share a prefix base: REF=A ALT=ATAC means insert TAC after pos
	// DEL: REF=ATAC ALT=A means delete TAC after pos
	refLen := len(v.Ref)
	altLen := len(v.Alt)

	if altLen > refLen {
		// Insertion
		insertedSeq := alt[1:] // skip shared prefix base on coding strand

		// Check for duplication: if inserted seq matches the preceding reference
		if isDuplication(v, t, insertedSeq) {
			dupLen := len(insertedSeq)
			if dupLen == 1 {
				// Single base dup: position is the duplicated base
				dupPosStr := genomicToHGVScPos(dupGenomicPos(v, t), t)
				return fmt.Sprintf("c.%sdup", dupPosStr)
			}
			// Multi-base dup
			dupStart, dupEnd := dupGenomicRange(v, t)
			startStr := genomicToHGVScPos(dupStart, t)
			endStr := genomicToHGVScPos(dupEnd, t)
			return fmt.Sprintf("c.%s_%sdup", startStr, endStr)
		}

		// Plain insertion: between the prefix base position and next position
		// Position is between the two flanking bases
		insAfterPos := v.Pos
		insBeforePos := v.Pos + 1
		if t.IsReverseStrand() {
			// On reverse strand, the insertion is before v.Pos in transcript terms
			insAfterPos = v.Pos + 1
			insBeforePos = v.Pos
		}
		afterStr := genomicToHGVScPos(insAfterPos, t)
		beforeStr := genomicToHGVScPos(insBeforePos, t)
		return fmt.Sprintf("c.%s_%sins%s", afterStr, beforeStr, insertedSeq)
	}

	// Deletion
	if refLen > altLen {
		// Deleted bases start after the shared prefix
		delStartGenomic := v.Pos + int64(altLen)
		delEndGenomic := v.Pos + int64(refLen) - 1

		if t.IsReverseStrand() {
			// For reverse strand, swap start/end for transcript order
			delStartGenomic, delEndGenomic = delEndGenomic, delStartGenomic
		}

		delStartStr := genomicToHGVScPos(delStartGenomic, t)
		if delStartGenomic == delEndGenomic {
			return fmt.Sprintf("c.%sdel", delStartStr)
		}
		delEndStr := genomicToHGVScPos(delEndGenomic, t)
		return fmt.Sprintf("c.%s_%sdel", delStartStr, delEndStr)
	}

	// MNV (same length ref/alt, len > 1) — treat as delins
	if refLen > 1 {
		endPosStr := genomicToHGVScPos(v.Pos+int64(refLen)-1, t)
		if t.IsReverseStrand() {
			startPosStr, endPosStr = genomicToHGVScPos(v.Pos+int64(refLen)-1, t), startPosStr
		}
		return fmt.Sprintf("c.%s_%sdelins%s", startPosStr, endPosStr, alt)
	}

	return ""
}

// genomicToHGVScPos converts a genomic position to an HGVSc position string.
// Returns strings like "76", "88+1", "89-2", "-14", "*6".
func genomicToHGVScPos(pos int64, t *cache.Transcript) string {
	if !t.IsProteinCoding() {
		return ""
	}

	// Check if position is in an exon
	exon := t.FindExon(pos)
	if exon != nil {
		return exonicHGVScPos(pos, exon, t)
	}

	// Intronic position — find flanking exons
	return intronicHGVScPos(pos, t)
}

// exonicHGVScPos returns the HGVSc position for an exonic position.
func exonicHGVScPos(pos int64, exon *cache.Exon, t *cache.Transcript) string {
	// Check if in CDS
	if exon.IsCoding() && pos >= exon.CDSStart && pos <= exon.CDSEnd {
		cdsPos := GenomicToCDS(pos, t)
		if cdsPos > 0 {
			return fmt.Sprintf("%d", cdsPos)
		}
	}

	// 5'UTR or 3'UTR
	if t.IsForwardStrand() {
		if pos < t.CDSStart {
			return fiveprimeUTRPos(pos, t)
		}
		if pos > t.CDSEnd {
			return threeprimeUTRPos(pos, t)
		}
	} else {
		if pos > t.CDSEnd {
			return fiveprimeUTRPos(pos, t)
		}
		if pos < t.CDSStart {
			return threeprimeUTRPos(pos, t)
		}
	}

	return ""
}

// fiveprimeUTRPos returns the HGVSc position for a 5'UTR position (e.g., "-14").
// Counts exonic bases from the position to CDS start.
func fiveprimeUTRPos(pos int64, t *cache.Transcript) string {
	dist := exonicDistance(pos, t.CDSStart, t) // forward strand
	if t.IsReverseStrand() {
		dist = exonicDistance(t.CDSEnd, pos, t) // reverse strand: CDSEnd is 5' end
	}
	if dist < 0 {
		dist = -dist
	}
	return fmt.Sprintf("-%d", dist)
}

// threeprimeUTRPos returns the HGVSc position for a 3'UTR position (e.g., "*6").
// Counts exonic bases from CDS end to the position.
func threeprimeUTRPos(pos int64, t *cache.Transcript) string {
	dist := exonicDistance(t.CDSEnd, pos, t) // forward strand
	if t.IsReverseStrand() {
		dist = exonicDistance(pos, t.CDSStart, t)
	}
	if dist < 0 {
		dist = -dist
	}
	return fmt.Sprintf("*%d", dist)
}

// exonicDistance counts the number of exonic bases between two genomic positions.
// Both positions are inclusive. Only counts bases that fall within exons.
// from < to in genomic coordinates.
func exonicDistance(from, to int64, t *cache.Transcript) int64 {
	if from > to {
		from, to = to, from
	}
	var dist int64
	for _, exon := range t.Exons {
		// Overlap between [from, to] and [exon.Start, exon.End]
		overlapStart := from
		if exon.Start > overlapStart {
			overlapStart = exon.Start
		}
		overlapEnd := to
		if exon.End < overlapEnd {
			overlapEnd = exon.End
		}
		if overlapStart <= overlapEnd {
			dist += overlapEnd - overlapStart + 1
		}
	}
	// Subtract 1 because we want distance not span (from position is the anchor)
	if dist > 0 {
		dist--
	}
	return dist
}

// intronicHGVScPos computes the HGVSc position for an intronic position.
// Format: c.{boundary}+{offset} or c.{boundary}-{offset}
func intronicHGVScPos(pos int64, t *cache.Transcript) string {
	// Find flanking exons
	var upstreamExon, downstreamExon *cache.Exon
	for i := range t.Exons {
		exon := &t.Exons[i]
		if exon.End < pos {
			if upstreamExon == nil || exon.End > upstreamExon.End {
				upstreamExon = exon
			}
		}
		if exon.Start > pos {
			if downstreamExon == nil || exon.Start < downstreamExon.Start {
				downstreamExon = exon
			}
		}
	}

	if upstreamExon == nil && downstreamExon == nil {
		return ""
	}

	// Calculate distances to each flanking exon boundary
	var distToUpstream, distToDownstream int64
	if upstreamExon != nil {
		distToUpstream = pos - upstreamExon.End
	}
	if downstreamExon != nil {
		distToDownstream = downstreamExon.Start - pos
	}

	// Pick the closer exon boundary
	useUpstream := true
	if upstreamExon == nil {
		useUpstream = false
	} else if downstreamExon != nil && distToDownstream < distToUpstream {
		useUpstream = false
	}

	if t.IsForwardStrand() {
		if useUpstream {
			// After upstream exon end: c.{CDSpos}+{offset}
			boundaryPos := exonBoundaryHGVScPos(upstreamExon.End, upstreamExon, t)
			return fmt.Sprintf("%s+%d", boundaryPos, distToUpstream)
		}
		// Before downstream exon start: c.{CDSpos}-{offset}
		boundaryPos := exonBoundaryHGVScPos(downstreamExon.Start, downstreamExon, t)
		return fmt.Sprintf("%s-%d", boundaryPos, distToDownstream)
	}

	// Reverse strand: genomic upstream exon end = transcript 5' direction
	// On reverse strand, the exon with higher genomic coords is upstream in transcript
	if useUpstream {
		// pos is after upstreamExon.End genomically, meaning it's in the intron
		// on reverse strand, this is *before* the exon in transcript order
		// So: c.{CDSpos of exon.End}-{offset}  (approaching from 3' side)
		boundaryPos := exonBoundaryHGVScPos(upstreamExon.End, upstreamExon, t)
		return fmt.Sprintf("%s-%d", boundaryPos, distToUpstream)
	}
	// pos is before downstreamExon.Start genomically
	// on reverse strand, this is *after* the exon in transcript order
	// So: c.{CDSpos of exon.Start}+{offset}
	boundaryPos := exonBoundaryHGVScPos(downstreamExon.Start, downstreamExon, t)
	return fmt.Sprintf("%s+%d", boundaryPos, distToDownstream)
}

// exonBoundaryHGVScPos returns the HGVSc position string for an exon boundary.
// The boundary is the exon position closest to the intron.
func exonBoundaryHGVScPos(genomicPos int64, exon *cache.Exon, t *cache.Transcript) string {
	// If the boundary is in the CDS, use CDS position
	if exon.IsCoding() && genomicPos >= exon.CDSStart && genomicPos <= exon.CDSEnd {
		cdsPos := GenomicToCDS(genomicPos, t)
		if cdsPos > 0 {
			return fmt.Sprintf("%d", cdsPos)
		}
	}

	// If in 5'UTR region
	if t.IsForwardStrand() {
		if genomicPos < t.CDSStart {
			return fiveprimeUTRPos(genomicPos, t)
		}
		if genomicPos > t.CDSEnd {
			return threeprimeUTRPos(genomicPos, t)
		}
	} else {
		if genomicPos > t.CDSEnd {
			return fiveprimeUTRPos(genomicPos, t)
		}
		if genomicPos < t.CDSStart {
			return threeprimeUTRPos(genomicPos, t)
		}
	}

	return ""
}

// isDuplication checks if an insertion duplicates the preceding reference sequence.
func isDuplication(v *vcf.Variant, t *cache.Transcript, insertedSeq string) bool {
	if len(insertedSeq) == 0 {
		return false
	}

	// For forward strand, check if the bases before the insertion match
	// For reverse strand, we need to check in the coding direction
	// Simple heuristic: check if the genomic bases preceding the insertion
	// match the inserted sequence
	_ = t // May use transcript sequence for more accurate check in future

	// Check the CDS sequence if available
	if len(t.CDSSequence) == 0 {
		return false
	}

	cdsPos := GenomicToCDS(v.Pos, t)
	if cdsPos < 1 {
		return false
	}

	cdsIdx := int(cdsPos - 1)
	seqLen := len(insertedSeq)

	// Check if the bases at [cdsIdx-seqLen+1 .. cdsIdx] match insertedSeq
	startIdx := cdsIdx - seqLen + 1
	if startIdx < 0 {
		return false
	}
	endIdx := cdsIdx + 1
	if endIdx > len(t.CDSSequence) {
		return false
	}

	precedingSeq := strings.ToUpper(t.CDSSequence[startIdx:endIdx])
	return precedingSeq == strings.ToUpper(insertedSeq)
}

// dupGenomicPos returns the genomic position of the duplicated base for single-base dups.
func dupGenomicPos(v *vcf.Variant, t *cache.Transcript) int64 {
	// The duplicated base is at v.Pos (the shared prefix base)
	return v.Pos
}

// dupGenomicRange returns the genomic start and end of the duplicated region.
func dupGenomicRange(v *vcf.Variant, t *cache.Transcript) (int64, int64) {
	dupLen := int64(len(v.Alt) - len(v.Ref))
	// Duplicated region ends at v.Pos and extends back dupLen bases
	start := v.Pos - dupLen + 1
	end := v.Pos
	if t.IsReverseStrand() {
		return end, start
	}
	return start, end
}
