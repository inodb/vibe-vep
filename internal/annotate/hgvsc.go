package annotate

import (
	"strconv"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// FormatHGVSc formats the HGVS coding DNA notation for a variant on a transcript.
// Returns empty string for non-coding transcripts or upstream/downstream variants.
//
// The function is optimized for minimal allocations on CDS-based paths (SNV, deletion,
// insertion) by using stack-allocated buffers and lazy computation of reverse complement.
func FormatHGVSc(v *vcf.Variant, t *cache.Transcript, result *ConsequenceResult) string {
	if t == nil || v == nil || result == nil {
		return ""
	}

	// Skip upstream/downstream/intergenic
	switch result.Consequence {
	case ConsequenceUpstreamGene, ConsequenceDownstreamGene, ConsequenceIntergenicVariant:
		return ""
	}

	// Determine prefix: c. for coding, n. for non-coding
	prefix := "c."
	if !t.IsProteinCoding() {
		prefix = "n."
	}

	// SNV (and MNV — same length ref/alt): compute RC and return immediately.
	// This avoids computing RC for the more common indel CDS paths.
	if !v.IsIndel() {
		ref := v.Ref
		alt := v.Alt
		if t.IsReverseStrand() {
			ref = ReverseComplement(ref)
			alt = ReverseComplement(alt)
		}
		startPosStr := genomicToHGVScPos(v.Pos, t)
		if startPosStr == "" {
			return ""
		}
		return prefix + startPosStr + ref + ">" + alt
	}

	// Indel handling — defer RC computation to branches that need it.
	refLen := len(v.Ref)
	altLen := len(v.Alt)

	if altLen > refLen {
		return formatHGVScInsertion(v, t, prefix, refLen, altLen)
	}

	if refLen > altLen {
		return formatHGVScDeletion(v, t, prefix, refLen, altLen)
	}

	// MNV (same length ref/alt, len > 1) — treat as delins.
	// Note: unreachable because !v.IsIndel() is true for same-length variants.
	if refLen > 1 {
		ref := v.Ref
		alt := v.Alt
		if t.IsReverseStrand() {
			ref = ReverseComplement(ref)
			alt = ReverseComplement(alt)
		}
		startPosStr := genomicToHGVScPos(v.Pos, t)
		endPosStr := genomicToHGVScPos(v.Pos+int64(refLen)-1, t)
		if t.IsReverseStrand() {
			startPosStr, endPosStr = genomicToHGVScPos(v.Pos+int64(refLen)-1, t), startPosStr
		}
		return prefix + startPosStr + "_" + endPosStr + "delins" + alt
	}

	return ""
}

// formatHGVScInsertion handles the insertion path of FormatHGVSc.
// It computes the inserted sequence on the coding strand using a stack buffer
// to avoid allocations from ReverseComplement.
func formatHGVScInsertion(v *vcf.Variant, t *cache.Transcript, prefix string, refLen, altLen int) string {
	insLen := altLen - refLen

	// Compute inserted sequence on coding strand into a stack buffer.
	// For reverse strand, we write the reverse complement directly to avoid
	// allocating a string from ReverseComplement().
	var seqStack [64]byte
	var seq []byte
	if insLen <= len(seqStack) {
		seq = seqStack[:insLen]
	} else {
		seq = make([]byte, insLen)
	}

	if t.IsReverseStrand() {
		// RC(v.Alt)[:altLen-refLen] == RC(v.Alt[refLen:])
		genomicInserted := v.Alt[refLen:]
		for i := 0; i < insLen; i++ {
			seq[i] = Complement(genomicInserted[insLen-1-i])
		}
	} else {
		copy(seq, v.Alt[refLen:])
	}

	// 3' shift insertion in CDS space, then check dup at shifted position
	cdsPos := GenomicToCDS(v.Pos, t)
	if cdsPos > 0 && len(t.CDSSequence) > 0 {
		cdsIdx := int(cdsPos - 1) // 0-based anchor index
		if t.IsReverseStrand() {
			cdsIdx-- // RC makes shared base a suffix; insertion is before CDS position
		}
		if cdsIdx < 0 {
			cdsIdx = 0
		}

		// Shift insertion in place (modifies seq)
		shiftedIdx := shiftInsertionBuf(seq, cdsIdx, t.CDSSequence)
		seqLen := len(seq)

		// Check dup: inserted bases match preceding bases at shifted position.
		// Go optimizes string([]byte) == string comparisons to avoid allocation.
		dupStart := shiftedIdx - seqLen + 1
		if dupStart >= 0 && shiftedIdx+1 <= len(t.CDSSequence) &&
			t.CDSSequence[dupStart:shiftedIdx+1] == string(seq) {
			return cdsPosRangeStr(dupStart+1, shiftedIdx+1, "dup")
		}

		// Check dup: inserted bases match following bases at shifted position
		afterStart := shiftedIdx + 1
		afterEnd := afterStart + seqLen
		if afterStart >= 0 && afterEnd <= len(t.CDSSequence) &&
			t.CDSSequence[afterStart:afterEnd] == string(seq) {
			return cdsPosRangeStr(afterStart+1, afterEnd, "dup")
		}

		// Plain insertion at shifted CDS position
		if shiftedIdx+2 <= len(t.CDSSequence) {
			var buf [128]byte
			n := copy(buf[:], "c.")
			n += putInt64(buf[n:], int64(shiftedIdx+1))
			buf[n] = '_'
			n++
			n += putInt64(buf[n:], int64(shiftedIdx+2))
			n += copy(buf[n:], "ins")
			n += copy(buf[n:], seq)
			return string(buf[:n])
		}
	}

	// Fall back to genomic-based positions for non-CDS insertions.
	// These paths use string concat (less common, not on the hot path).
	insertedSeq := string(seq)

	dup := checkDuplication(v, t, insertedSeq)
	if dup.isDup {
		if dup.cdsStart == dup.cdsEnd {
			return "c." + strconv.FormatInt(dup.cdsStart, 10) + "dup"
		}
		return "c." + strconv.FormatInt(dup.cdsStart, 10) + "_" + strconv.FormatInt(dup.cdsEnd, 10) + "dup"
	}

	// Check for duplication at splice junctions where one flanking position
	// is intronic (CDS=0) and the other is exonic. The CDS-based dup check
	// above bails when the anchor is intronic, but the inserted base may
	// still duplicate the exon boundary base.
	if sjDup := checkSpliceJunctionDup(v, t, insertedSeq); sjDup.isDup {
		if sjDup.cdsStart == sjDup.cdsEnd {
			return "c." + strconv.FormatInt(sjDup.cdsStart, 10) + "dup"
		}
		return "c." + strconv.FormatInt(sjDup.cdsStart, 10) + "_" + strconv.FormatInt(sjDup.cdsEnd, 10) + "dup"
	}

	insAfterPos := v.Pos
	insBeforePos := v.Pos + 1
	if t.IsReverseStrand() {
		insAfterPos = v.Pos + 1
		insBeforePos = v.Pos
	}
	afterStr := genomicToHGVScPos(insAfterPos, t)
	beforeStr := genomicToHGVScPos(insBeforePos, t)
	return prefix + afterStr + "_" + beforeStr + "ins" + insertedSeq
}

// formatHGVScDeletion handles the deletion (and delins) path of FormatHGVSc.
// It avoids calling ReverseComplement for the ref/alt at the top since the
// deletion CDS path only needs CDS positions and possibly RC of the extra alt bases.
func formatHGVScDeletion(v *vcf.Variant, t *cache.Transcript, prefix string, refLen, altLen int) string {
	// Compute actual shared prefix length on genomic strand
	sharedLen := 0
	for sharedLen < refLen && sharedLen < altLen && v.Ref[sharedLen] == v.Alt[sharedLen] {
		sharedLen++
	}
	if sharedLen == 0 {
		sharedLen = 1
	}
	if sharedLen > altLen {
		sharedLen = altLen
	}

	delStartGenomic := v.Pos + int64(sharedLen)
	delEndGenomic := v.Pos + int64(refLen) - 1
	extraAlt := ""
	if sharedLen < altLen {
		extraAlt = v.Alt[sharedLen:]
	}

	// Clip shared suffix between remaining ref and extra alt.
	if len(extraAlt) > 0 {
		remainRef := v.Ref[sharedLen:]
		sharedSuffix := 0
		for sharedSuffix < len(remainRef) && sharedSuffix < len(extraAlt) &&
			remainRef[len(remainRef)-1-sharedSuffix] == extraAlt[len(extraAlt)-1-sharedSuffix] {
			sharedSuffix++
		}
		if sharedSuffix > 0 {
			delEndGenomic -= int64(sharedSuffix)
			extraAlt = extraAlt[:len(extraAlt)-sharedSuffix]
		}
	}

	if t.IsReverseStrand() {
		delStartGenomic, delEndGenomic = delEndGenomic, delStartGenomic
	}

	// Try CDS-based positioning (hot path — optimized with stack buffers)
	delStartCDS := GenomicToCDS(delStartGenomic, t)
	delEndCDS := GenomicToCDS(delEndGenomic, t)
	if delStartCDS > 0 && delEndCDS > 0 && len(t.CDSSequence) > 0 {
		if len(extraAlt) > 0 {
			// Delins: no 3' shift per HGVS convention.
			// Write RC of extra alt directly into the output buffer for reverse strand.
			var buf [128]byte
			n := copy(buf[:], "c.")
			n += putInt64(buf[n:], delStartCDS)
			if delStartCDS != delEndCDS {
				buf[n] = '_'
				n++
				n += putInt64(buf[n:], delEndCDS)
			}
			n += copy(buf[n:], "delins")
			if t.IsReverseStrand() {
				// Write reverse complement of extraAlt directly into buffer
				eLen := len(extraAlt)
				for i := 0; i < eLen; i++ {
					buf[n+i] = Complement(extraAlt[eLen-1-i])
				}
				n += eLen
			} else {
				n += copy(buf[n:], extraAlt)
			}
			return string(buf[:n])
		}
		// Pure deletion: apply 3' shift
		sStart, sEnd := shiftDeletionThreePrime(int(delStartCDS-1), int(delEndCDS-1), t.CDSSequence)
		return cdsPosRangeStr(sStart+1, sEnd+1, "del")
	}

	// Fall back to genomic-based positions for intronic/UTR/non-coding deletions.
	// Convert extra alt bases to coding strand (only needed for fallback path).
	codingExtraAlt := extraAlt
	if t.IsReverseStrand() && len(extraAlt) > 0 {
		codingExtraAlt = ReverseComplement(extraAlt)
	}

	delStartStr := genomicToHGVScPos(delStartGenomic, t)
	if len(codingExtraAlt) > 0 {
		if delStartGenomic == delEndGenomic {
			return prefix + delStartStr + "delins" + codingExtraAlt
		}
		delEndStr := genomicToHGVScPos(delEndGenomic, t)
		return prefix + delStartStr + "_" + delEndStr + "delins" + codingExtraAlt
	}
	if delStartGenomic == delEndGenomic {
		return prefix + delStartStr + "del"
	}
	delEndStr := genomicToHGVScPos(delEndGenomic, t)
	return prefix + delStartStr + "_" + delEndStr + "del"
}

// cdsPosRangeStr builds a CDS position string like "c.34del", "c.34_36del",
// "c.35dup", etc. using a stack-allocated buffer for zero intermediate allocations.
func cdsPosRangeStr(start, end int, suffix string) string {
	var buf [64]byte
	n := copy(buf[:], "c.")
	n += putInt64(buf[n:], int64(start))
	if start != end {
		buf[n] = '_'
		n++
		n += putInt64(buf[n:], int64(end))
	}
	n += copy(buf[n:], suffix)
	return string(buf[:n])
}

// shiftInsertionBuf shifts an insertion rightward (3' direction) in CDS space.
// It operates on the provided mutable byte slice, avoiding allocations.
// Returns the new anchor index.
func shiftInsertionBuf(seq []byte, cdsAnchorIdx int, cdsSeq string) int {
	idx := cdsAnchorIdx
	for idx+1 < len(cdsSeq) && cdsSeq[idx+1] == seq[0] {
		first := seq[0]
		copy(seq, seq[1:])
		seq[len(seq)-1] = first
		idx++
	}
	return idx
}

// genomicToHGVScPos converts a genomic position to an HGVSc position string.
// For coding transcripts, returns strings like "76", "88+1", "89-2", "-14", "*6".
// For non-coding transcripts, returns transcript-relative positions like "127", "42+5".
func genomicToHGVScPos(pos int64, t *cache.Transcript) string {
	// Check if position is in an exon
	exon := t.FindExon(pos)
	if exon != nil {
		if t.IsProteinCoding() {
			return exonicHGVScPos(pos, exon, t)
		}
		return nonCodingExonicPos(pos, t)
	}

	// Intronic position — find flanking exons
	return intronicHGVScPos(pos, t)
}

// nonCodingExonicPos returns the transcript position for an exonic position in a
// non-coding transcript. Position 1 is the first exonic base at the 5' end.
func nonCodingExonicPos(pos int64, t *cache.Transcript) string {
	txPos := GenomicToTranscriptPos(pos, t)
	if txPos > 0 {
		return strconv.FormatInt(txPos, 10)
	}
	return ""
}

// exonicHGVScPos returns the HGVSc position for an exonic position.
func exonicHGVScPos(pos int64, exon *cache.Exon, t *cache.Transcript) string {
	// Check if in CDS
	if exon.IsCoding() && pos >= exon.CDSStart && pos <= exon.CDSEnd {
		cdsPos := GenomicToCDS(pos, t)
		if cdsPos > 0 {
			return strconv.FormatInt(cdsPos, 10)
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
	return "-" + strconv.FormatInt(dist, 10)
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
	return "*" + strconv.FormatInt(dist, 10)
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
// Uses FindNearestExonIdx for O(log n) flanking exon lookup.
func intronicHGVScPos(pos int64, t *cache.Transcript) string {
	idx := t.FindNearestExonIdx(pos)
	n := len(t.Exons)
	if n == 0 {
		return ""
	}

	// Find flanking exons by checking idx and its neighbors.
	// Exons are sorted ascending by Start. An intronic position has
	// an upstream exon (End < pos) and a downstream exon (Start > pos).
	var upstreamExon, downstreamExon *cache.Exon
	for _, i := range [3]int{idx - 1, idx, idx + 1} {
		if i < 0 || i >= n {
			continue
		}
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
			boundaryPos := exonBoundaryHGVScPos(upstreamExon.End, upstreamExon, t)
			return boundaryPos + "+" + strconv.FormatInt(distToUpstream, 10)
		}
		boundaryPos := exonBoundaryHGVScPos(downstreamExon.Start, downstreamExon, t)
		return boundaryPos + "-" + strconv.FormatInt(distToDownstream, 10)
	}

	// Reverse strand
	if useUpstream {
		boundaryPos := exonBoundaryHGVScPos(upstreamExon.End, upstreamExon, t)
		return boundaryPos + "-" + strconv.FormatInt(distToUpstream, 10)
	}
	boundaryPos := exonBoundaryHGVScPos(downstreamExon.Start, downstreamExon, t)
	return boundaryPos + "+" + strconv.FormatInt(distToDownstream, 10)
}

// exonBoundaryHGVScPos returns the HGVSc position string for an exon boundary.
// The boundary is the exon position closest to the intron.
func exonBoundaryHGVScPos(genomicPos int64, exon *cache.Exon, t *cache.Transcript) string {
	// Non-coding transcripts use transcript-relative position
	if !t.IsProteinCoding() {
		return nonCodingExonicPos(genomicPos, t)
	}

	// If the boundary is in the CDS, use CDS position
	if exon.IsCoding() && genomicPos >= exon.CDSStart && genomicPos <= exon.CDSEnd {
		cdsPos := GenomicToCDS(genomicPos, t)
		if cdsPos > 0 {
			return strconv.FormatInt(cdsPos, 10)
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
	seq := []byte(insertedSeq)
	idx := cdsAnchorIdx
	for idx+1 < len(cdsSeq) && cdsSeq[idx+1] == seq[0] {
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
	for delEnd+1 < len(cdsSeq) && cdsSeq[delStart] == cdsSeq[delEnd+1] {
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

// checkSpliceJunctionDup detects duplications at exon-intron boundaries where
// the VCF anchor is intronic (CDS=0) but the adjacent exonic base matches the
// inserted sequence. This handles two patterns:
//   - Forward strand c.X-1_XinsN: anchor (v.Pos) intronic, v.Pos+1 exonic
//   - Reverse strand c.X_X+1insN: anchor (v.Pos) exonic, v.Pos+1 intronic
//
// In both cases, the CDS-based dup detection in formatHGVScInsertion skips
// these because GenomicToCDS returns 0 for the intronic position.
func checkSpliceJunctionDup(v *vcf.Variant, t *cache.Transcript, insertedSeq string) dupResult {
	if len(insertedSeq) == 0 || len(t.CDSSequence) == 0 {
		return dupResult{}
	}

	// Try both flanking positions to find which one is exonic
	cdsLeft := GenomicToCDS(v.Pos, t)
	cdsRight := GenomicToCDS(v.Pos+1, t)

	// We only handle the case where exactly one side is exonic
	if (cdsLeft > 0) == (cdsRight > 0) {
		return dupResult{}
	}

	seqLen := len(insertedSeq)

	if cdsRight > 0 {
		// Forward strand pattern: anchor is intronic, next position is exonic.
		// The exonic boundary is v.Pos+1, mapping to cdsRight.
		// Check if CDS bases starting at cdsRight match the inserted sequence.
		startIdx := int(cdsRight - 1) // 0-based
		endIdx := startIdx + seqLen
		if endIdx <= len(t.CDSSequence) && t.CDSSequence[startIdx:endIdx] == insertedSeq {
			return dupResult{
				isDup:    true,
				cdsStart: cdsRight,
				cdsEnd:   cdsRight + int64(seqLen) - 1,
			}
		}
	} else {
		// Reverse strand pattern: anchor is exonic, next position is intronic.
		// The exonic boundary is v.Pos, mapping to cdsLeft.
		// Check if CDS bases ending at cdsLeft match the inserted sequence.
		endIdx := int(cdsLeft) // exclusive, since cdsLeft is 1-based and we want [start, cdsLeft]
		startIdx := endIdx - seqLen
		if startIdx >= 0 && t.CDSSequence[startIdx:endIdx] == insertedSeq {
			return dupResult{
				isDup:    true,
				cdsStart: int64(startIdx + 1),
				cdsEnd:   cdsLeft,
			}
		}
	}

	return dupResult{}
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

	// Check if inserted bases match the bases at the anchor position
	// (the anchor and preceding bases in CDS).
	startIdx := cdsIdx - seqLen + 1
	if startIdx >= 0 && cdsIdx+1 <= len(t.CDSSequence) {
		if t.CDSSequence[startIdx:cdsIdx+1] == insertedSeq {
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
		if t.CDSSequence[afterStart:afterEnd] == insertedSeq {
			return dupResult{
				isDup:    true,
				cdsStart: int64(afterStart + 1), // 1-based
				cdsEnd:   int64(afterEnd),
			}
		}
	}

	return dupResult{}
}
