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
		// VCF indels share a prefix base on the genomic strand.
		// After reverse-complementing for reverse-strand transcripts, the
		// shared base becomes the suffix, so the inserted bases are at the
		// start of alt rather than after position 1.
		var insertedSeq string
		if t.IsReverseStrand() {
			insertedSeq = alt[:len(alt)-len(ref)] // shared base is at the end
		} else {
			insertedSeq = alt[len(ref):] // shared base is at the start
		}

		// 3' shift insertion in CDS space, then check dup at shifted position
		cdsPos := GenomicToCDS(v.Pos, t)
		if cdsPos > 0 && len(t.CDSSequence) > 0 {
			cdsIdx := int(cdsPos - 1) // 0-based anchor index
			shiftedSeq, shiftedIdx := shiftInsertionThreePrime(insertedSeq, cdsIdx, t.CDSSequence)
			seqLen := len(shiftedSeq)
			upperCDS := strings.ToUpper(t.CDSSequence)
			upperIns := strings.ToUpper(shiftedSeq)

			// Check dup: inserted bases match preceding bases at shifted position
			dupStart := shiftedIdx - seqLen + 1
			if dupStart >= 0 && shiftedIdx+1 <= len(upperCDS) &&
				upperCDS[dupStart:shiftedIdx+1] == upperIns {
				if seqLen == 1 {
					return fmt.Sprintf("c.%ddup", shiftedIdx+1)
				}
				return fmt.Sprintf("c.%d_%ddup", dupStart+1, shiftedIdx+1)
			}

			// Check dup: inserted bases match following bases at shifted position
			afterStart := shiftedIdx + 1
			afterEnd := afterStart + seqLen
			if afterStart >= 0 && afterEnd <= len(upperCDS) &&
				upperCDS[afterStart:afterEnd] == upperIns {
				if seqLen == 1 {
					return fmt.Sprintf("c.%ddup", afterStart+1)
				}
				return fmt.Sprintf("c.%d_%ddup", afterStart+1, afterEnd)
			}

			// Plain insertion at shifted CDS position
			if shiftedIdx+2 <= len(t.CDSSequence) {
				return fmt.Sprintf("c.%d_%dins%s", shiftedIdx+1, shiftedIdx+2, shiftedSeq)
			}
		}

		// Fall back to genomic-based positions for non-CDS insertions
		dup := checkDuplication(v, t, insertedSeq)
		if dup.isDup {
			if dup.cdsStart == dup.cdsEnd {
				return fmt.Sprintf("c.%ddup", dup.cdsStart)
			}
			return fmt.Sprintf("c.%d_%ddup", dup.cdsStart, dup.cdsEnd)
		}

		insAfterPos := v.Pos
		insBeforePos := v.Pos + 1
		if t.IsReverseStrand() {
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

		// Try 3' shifting in CDS space
		delStartCDS := GenomicToCDS(delStartGenomic, t)
		delEndCDS := GenomicToCDS(delEndGenomic, t)
		if delStartCDS > 0 && delEndCDS > 0 && len(t.CDSSequence) > 0 {
			sStart, sEnd := shiftDeletionThreePrime(int(delStartCDS-1), int(delEndCDS-1), t.CDSSequence)
			if sStart == sEnd {
				return fmt.Sprintf("c.%ddel", sStart+1)
			}
			return fmt.Sprintf("c.%d_%ddel", sStart+1, sEnd+1)
		}

		// Fall back to genomic-based positions for intronic/UTR deletions
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

// shiftInsertionThreePrime shifts an insertion rightward (3' direction) in CDS space.
// cdsAnchorIdx is the 0-based CDS index of the VCF anchor base (insertion is after this position).
// Returns the shifted inserted sequence and new anchor index.
func shiftInsertionThreePrime(insertedSeq string, cdsAnchorIdx int, cdsSeq string) (string, int) {
	seq := []byte(strings.ToUpper(insertedSeq))
	idx := cdsAnchorIdx
	upper := strings.ToUpper(cdsSeq)
	for idx+1 < len(upper) && upper[idx+1] == seq[0] {
		// Rotate: move first base to end, advance anchor
		first := seq[0]
		copy(seq, seq[1:])
		seq[len(seq)-1] = first
		idx++
	}
	return string(seq), idx
}

// shiftDeletionThreePrime shifts a deletion rightward (3' direction) in CDS space.
// delStart and delEnd are 0-based inclusive CDS indices of the deleted bases.
// Returns the shifted start and end indices.
func shiftDeletionThreePrime(delStart, delEnd int, cdsSeq string) (int, int) {
	upper := strings.ToUpper(cdsSeq)
	for delEnd+1 < len(upper) && upper[delStart] == upper[delEnd+1] {
		delStart++
		delEnd++
	}
	return delStart, delEnd
}

// dupResult holds the CDS position range of a detected duplication.
type dupResult struct {
	isDup    bool
	cdsStart int64 // 1-based CDS position of first duplicated base
	cdsEnd   int64 // 1-based CDS position of last duplicated base
}

// checkDuplication checks if an insertion duplicates adjacent reference bases.
// Per HGVS convention, the inserted sequence must match either the bases
// immediately before or immediately after the insertion point.
// Returns the CDS positions of the duplicated bases for direct HGVS formatting.
func checkDuplication(v *vcf.Variant, t *cache.Transcript, insertedSeq string) dupResult {
	if len(insertedSeq) == 0 || len(t.CDSSequence) == 0 {
		return dupResult{}
	}

	cdsPos := GenomicToCDS(v.Pos, t)
	if cdsPos < 1 {
		return dupResult{}
	}

	cdsIdx := int(cdsPos - 1) // 0-based index of the VCF anchor base
	seqLen := len(insertedSeq)
	upper := strings.ToUpper(insertedSeq)

	// Check if inserted bases match the bases at the anchor position
	// (the anchor and preceding bases in CDS).
	startIdx := cdsIdx - seqLen + 1
	if startIdx >= 0 && cdsIdx+1 <= len(t.CDSSequence) {
		if strings.ToUpper(t.CDSSequence[startIdx:cdsIdx+1]) == upper {
			return dupResult{
				isDup:    true,
				cdsStart: int64(startIdx + 1), // 1-based
				cdsEnd:   cdsPos,
			}
		}
	}

	// Check if inserted bases match the bases immediately after the anchor
	// (the bases right after the insertion point in CDS).
	afterStart := cdsIdx + 1
	afterEnd := afterStart + seqLen
	if afterStart >= 0 && afterEnd <= len(t.CDSSequence) {
		if strings.ToUpper(t.CDSSequence[afterStart:afterEnd]) == upper {
			return dupResult{
				isDup:    true,
				cdsStart: int64(afterStart + 1), // 1-based
				cdsEnd:   int64(afterEnd),
			}
		}
	}

	return dupResult{}
}
