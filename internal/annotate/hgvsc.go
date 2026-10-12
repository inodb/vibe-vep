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

	// A variant whose first base lies outside the transcript has no c.
	// coordinate to start from, and VEP emits no HGVSc for one: RAD51C
	// 17:56769924 deletes 91 bases ending inside exon 1 and VEP leaves the
	// field empty while still reporting coding_sequence_variant,
	// 5_prime_UTR_variant. Extrapolating a coordinate from the nearest exon,
	// as genomicToHGVScPos will happily do, invents a position the transcript
	// does not have.
	if !t.Contains(v.Pos) {
		return ""
	}

	// Determine prefix: c. for coding, n. for non-coding
	prefix := "c."
	if !t.IsProteinCoding() {
		prefix = "n."
	}

	// SNV and MNV (same length ref/alt): compute RC and format.
	if !v.IsIndel() {
		// No change (ref==alt): HGVS uses "=" notation
		if v.Ref == v.Alt {
			startPosStr := genomicToHGVScPos(v.Pos, t)
			if startPosStr == "" {
				return ""
			}
			return prefix + startPosStr + "="
		}
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
		// SNV: single nucleotide substitution (c.123A>G)
		if len(v.Ref) == 1 {
			return prefix + startPosStr + ref + ">" + alt
		}
		// MNV: multi-nucleotide changes must use delins per HGVS
		// (substitution is defined as exactly one nucleotide replaced by one)
		endPosStr := genomicToHGVScPos(v.Pos+int64(len(v.Ref))-1, t)
		if t.IsReverseStrand() {
			startPosStr, endPosStr = endPosStr, startPosStr
		}
		if endPosStr == "" {
			return ""
		}
		// An MNV whose alternate allele is the reverse complement of the
		// reference is an inversion, and HGVS gives it its own operator:
		// ATRX X:76938207 AA>TT is c.2540_2541inv, not c.2540_2541delinsAA.
		// This holds on either strand, since reverse-complementing both
		// alleles preserves the relation. It accounts for every one of the
		// 161 inv calls in a VEP111 MSK-IMPACT MAF and nothing else: the
		// other 7,263 same-length multi-base changes there are delins.
		if isInversion(ref, alt) {
			return prefix + startPosStr + "_" + endPosStr + "inv"
		}
		return prefix + startPosStr + "_" + endPosStr + "delins" + alt
	}

	// Indel handling — defer RC computation to branches that need it.
	refLen := len(v.Ref)
	altLen := len(v.Alt)

	// A longer alternate allele is only an insertion when the reference
	// survives inside it. MAF rows typed INS are not always pure insertions -
	// NF1 17:29686003 is ref C alt GG - and the insertion path assumes the
	// alt begins with the ref, so it read that as inserting a single G and
	// dropped the C>G substitution entirely. VEP calls it c.8067delinsGG.
	if altLen > refLen && !isDelIns(v.Ref, v.Alt) {
		return formatHGVScInsertion(v, t, prefix, refLen, altLen, result)
	}

	return formatHGVScDeletion(v, t, prefix, refLen, altLen, result)
}

// formatHGVScInsertion handles the insertion path of FormatHGVSc.
// It computes the inserted sequence on the coding strand using a stack buffer
// to avoid allocations from ReverseComplement.
func formatHGVScInsertion(v *vcf.Variant, t *cache.Transcript, prefix string, refLen, altLen int, result *ConsequenceResult) string {
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

		// Codons are reported at the UNSHIFTED position, with the inserted
		// bases in the order the variant gives them. VEP writes the HGVS at
		// the 3'-shifted position but expands codons around where the change
		// actually sits: SMARCA2 c.705_707dup is -/CAG, the allele as it
		// appears in the MAF, not the -/GCA the same bases become once they
		// have been walked 41 bases down the repeat.
		//
		// This runs before the shift because shiftInsertionBuf rotates seq in
		// place. Unshifting the anchor afterwards is not enough on its own -
		// that was the bug, and it left the position right but the bases
		// rotated. The anchor here is the same value that subtracting the
		// offset used to produce.
		result.CodonChange = formatCodonChangeCDS(t.CDSSequence,
			int64(cdsIdx)+2, int64(cdsIdx)+1, string(seq))

		// Shift insertion in place (modifies seq), stops at exon boundary
		maxIdx := cdsExonEndIdx(cdsIdx, t)
		shiftedIdx := shiftInsertionBuf(seq, cdsIdx, t.CDSSequence, maxIdx)
		result.HGVSOffset = shiftedIdx - cdsIdx
		seqLen := len(seq)

		// The shift stops at the exon edge; HGVS does not (see
		// insertionShiftIntoIntron).
		if k, rotated, eEnd, okShift := insertionShiftIntoIntron(seq, shiftedIdx, t); okShift {
			if g := CDSToGenomic(int64(shiftedIdx+1), t); g != 0 {
				if exon := t.FindExon(g); exon != nil {
					result.HGVSOffset += k
					// A frameshifting insertion that ends up inside the intron
					// has no protein product to describe: VEP emits no HGVSp
					// and genome-nexus renders p.<residue><pos>fs from the
					// protein position (KDM6A c.225+1dup is p.K75fs, not
					// p.A76Gfs*6).
					//
					// An in-frame one still gets a description -- EGFL7
					// c.313+1_313+3dup is p.G103dup to VEP -- so the withdrawal
					// keys on the frame, not on landing in the intron.
					if insLen%3 != 0 {
						result.SuppressHGVSp = true
					}
					if result.SuppressHGVSp {
						pos := result.ProteinStart
						if pos < 1 {
							pos = result.ProteinPosition
						}
						if cd := GetCodon(t.CDSSequence, pos); len(cd) == 3 {
							if aa := TranslateCodon(cd); aa != 0 {
								result.AminoAcidsVEP = string(aa)
							}
						}
					}
					return formatInsertionInIntron(rotated, exon, t, eEnd, k)
				}
			}
		}

		// HGVS defines insertion and duplication as mutually exclusive: an
		// insertion is a change "where the insertion is not a copy of a sequence
		// immediately 5'", and a duplication is a copy inserted "directly 3' of
		// the original copy of that sequence". So whenever the inserted bases are
		// a tandem copy of what precedes them, "ins" is not a valid description.
		//
		// This deliberately diverges from Ensembl VEP, which only writes "dup"
		// when the 3' shift moved the variant. Across a VEP111 MSK-IMPACT MAF
		// 6,769 of VEP's 8,124 "ins" descriptions are tandem copies by the
		// reference sequence -- e.g. c.1663_1664insC inserted into a run of seven
		// C's, which HGVS requires be written c.1663dup. VEP also contradicts
		// itself on 541 rows, writing "ins" at the CDS level and "dup" at the
		// protein level for the same variant (ASXL2 c.353_354insCAG / p.S117dup).
		// Matching VEP here would cost nothing in accuracy but would emit
		// descriptions that fail to match ClinVar and LOVD, which normalize to
		// "dup". Concordance with vep111 on HGVSc is ~1.4pp lower as a result.

		// Check dup: inserted bases match preceding bases at shifted position.
		// Go optimizes string([]byte) == string comparisons to avoid allocation.
		dupStart := shiftedIdx - seqLen + 1
		if dupStart >= 0 && shiftedIdx+1 <= len(t.CDSSequence) &&
			t.CDSSequence[dupStart:shiftedIdx+1] == string(seq) {
			return cdsPosRangeStr(t, dupStart+1, shiftedIdx+1, "dup")
		}

		// Check dup: inserted bases match following bases at shifted position
		afterStart := shiftedIdx + 1
		afterEnd := afterStart + seqLen
		if afterStart >= 0 && afterEnd <= len(t.CDSSequence) &&
			t.CDSSequence[afterStart:afterEnd] == string(seq) {
			return cdsPosRangeStr(t, afterStart+1, afterEnd, "dup")
		}

		// Plain insertion at shifted CDS position
		if shiftedIdx+2 <= len(t.CDSSequence) {
			var buf [128]byte
			n := copy(buf[:], "c.")
			n += putInt64(buf[n:], cdsDisplayPos(int64(shiftedIdx+1), t))
			buf[n] = '_'
			n++
			n += putInt64(buf[n:], cdsDisplayPos(int64(shiftedIdx+2), t))
			n += copy(buf[n:], "ins")
			n += copy(buf[n:], seq)
			return string(buf[:n])
		}
		// Insertion at end of CDS: second flanking position is *1 (3'UTR)
		if shiftedIdx+1 == len(t.CDSSequence) {
			return prefix + strconv.FormatInt(cdsDisplayPos(int64(shiftedIdx+1), t), 10) + "_*1ins" + string(seq)
		}
	}

	// Fall back to genomic-based positions for non-CDS insertions.
	// These paths use string concat (less common, not on the hot path).
	insertedSeq := string(seq)

	// Only the CDS-based block above filled in the Codons string, so an
	// insertion anchored on the intronic side of a junction came out with
	// none at all - 59 rows, and with it the stop_gained test that reads the
	// altered codons. The codons are still well defined: the insertion sits
	// between two bases, and whichever of them is coding gives the anchor.
	//
	// seq is unrotated here because shiftInsertionBuf runs only inside that
	// block, so these are the bases as the variant gives them - which is what
	// the Codons field shows.
	if result.CodonChange == "" && len(t.CDSSequence) > 0 {
		insAfter, insBefore := v.Pos, v.Pos+1
		if t.IsReverseStrand() {
			insAfter, insBefore = v.Pos+1, v.Pos
		}
		anchor := GenomicToCDS(insAfter, t)
		if anchor < 1 {
			if next := GenomicToCDS(insBefore, t); next >= 1 {
				anchor = next - 1
			}
		}
		if anchor >= 1 {
			result.CodonChange = formatCodonChangeCDS(t.CDSSequence, anchor+1, anchor, insertedSeq)
		}
	}

	dup := checkDuplication(v, t, insertedSeq)
	if dup.isDup {
		if dup.cdsStart == dup.cdsEnd {
			return "c." + strconv.FormatInt(cdsDisplayPos(dup.cdsStart, t), 10) + "dup"
		}
		return "c." + strconv.FormatInt(cdsDisplayPos(dup.cdsStart, t), 10) + "_" + strconv.FormatInt(cdsDisplayPos(dup.cdsEnd, t), 10) + "dup"
	}

	// Check for duplication at splice junctions where one flanking position
	// is intronic (CDS=0) and the other is exonic. The CDS-based dup check
	// above bails when the anchor is intronic, but the inserted base may
	// still duplicate the exon boundary base.
	// A pure insertion inside the intron may shift 3' into the exon.
	if shifted, k, ok := insertionShiftIntronToExon(v, t, insertedSeq); ok {
		result.HGVSOffset = k
		return shifted
	}

	if sjDup := checkSpliceJunctionDup(v, t, insertedSeq); sjDup.isDup {
		if sjDup.cdsStart == sjDup.cdsEnd {
			return "c." + strconv.FormatInt(cdsDisplayPos(sjDup.cdsStart, t), 10) + "dup"
		}
		return "c." + strconv.FormatInt(cdsDisplayPos(sjDup.cdsStart, t), 10) + "_" + strconv.FormatInt(cdsDisplayPos(sjDup.cdsEnd, t), 10) + "dup"
	}

	// Both flanks intronic: shift and test for a duplication against the
	// cached intron flank, which the CDS-based checks above cannot see.
	if dupStr, ok := intronicInsertionDup(v, t, insertedSeq, result); ok {
		return dupStr
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

// isInversion reports whether alt is the reverse complement of ref, compared
// in place rather than against an allocated reverse complement.
func isInversion(ref, alt string) bool {
	if len(ref) != len(alt) || len(ref) < 2 {
		return false
	}
	for i := 0; i < len(ref); i++ {
		if alt[i] != Complement(ref[len(ref)-1-i]) {
			return false
		}
	}
	return true
}

// formatHGVScDeletion handles the deletion (and delins) path of FormatHGVSc.
// It avoids calling ReverseComplement for the ref/alt at the top since the
// deletion CDS path only needs CDS positions and possibly RC of the extra alt bases.
func formatHGVScDeletion(v *vcf.Variant, t *cache.Transcript, prefix string, refLen, altLen int, result *ConsequenceResult) string {
	// Compute actual shared prefix length on genomic strand
	sharedLen := 0
	for sharedLen < refLen && sharedLen < altLen && v.Ref[sharedLen] == v.Alt[sharedLen] {
		sharedLen++
	}
	// A VCF-style deletion carries an anchor base the ref and alt share, which
	// the loop above finds. When they share nothing the variant replaces every
	// base it covers, so clipping one anyway would drop the inserted allele and
	// turn a delins into a plain deletion: KIT ref=AGTGGAA alt=C came out as
	// c.1669_1674del where VEP has c.1667_1673delinsC. A pure MAF deletion has
	// no alt at all and is unaffected either way.
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
			codingExtra := extraAlt
			if t.IsReverseStrand() {
				codingExtra = ReverseComplement(extraAlt)
			}
			result.CodonChange = formatCodonChangeCDS(t.CDSSequence,
				delStartCDS, delEndCDS, codingExtra)
			// Write RC of extra alt directly into the output buffer for reverse strand.
			var buf [128]byte
			n := copy(buf[:], "c.")
			n += putInt64(buf[n:], cdsDisplayPos(delStartCDS, t))
			if delStartCDS != delEndCDS {
				buf[n] = '_'
				n++
				n += putInt64(buf[n:], cdsDisplayPos(delEndCDS, t))
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
		// A deletion whose 3' end reaches the exon edge may keep shifting into
		// the intron (see deletionShiftIntoIntron).
		if sc, ee, _, k, ok := deletionShiftIntoIntron(v, t); ok {
			result.HGVSOffset = int(sc - delStartCDS)
			if result.HGVSOffset < 0 {
				result.HGVSOffset = 0
			}
			// Codons are still reported, and at the UNSHIFTED position: VEP
			// gives aaG/aa for PTEN c.801+1del, the codon around c.801. Leaving
			// this out blanked the column for every shifted row.
			result.CodonChange = formatCodonChangeCDS(t.CDSSequence, delStartCDS, delEndCDS, "")
			if g := genomicSpanCodons(v, t, delStartCDS, delEndCDS); g != "" {
				result.CodonChange = g
			}
			return formatDeletionAcrossDonor(t, sc, ee, k)
		}

		// Pure deletion: apply 3' shift (stops at exon boundary)
		maxIdx := cdsExonEndIdx(int(delEndCDS-1), t)
		sStart, sEnd := shiftDeletionThreePrime(int(delStartCDS-1), int(delEndCDS-1), t.CDSSequence, maxIdx)
		result.HGVSOffset = sStart - int(delStartCDS-1)
		// Codons are expanded around the unshifted deleted range (see above).
		result.CodonChange = formatCodonChangeCDS(t.CDSSequence,
			int64(sStart+1-result.HGVSOffset), int64(sEnd+1-result.HGVSOffset), "")
		if g := genomicSpanCodons(v, t, int64(sStart+1-result.HGVSOffset),
			int64(sEnd+1-result.HGVSOffset)); g != "" {
			result.CodonChange = g
		}
		return cdsPosRangeStr(t, sStart+1, sEnd+1, "del")
	}

	// Before falling back, a purely intronic pure deletion may still 3'-shift
	// across the junction into the exon (see shiftedIntronicDelCDS).
	if len(extraAlt) == 0 {
		if sStart, sEnd, k := shiftedIntronicDelCDS(v, t); sStart > 0 {
			result.HGVSOffset = k
			return cdsPosRangeStr(t, int(sStart), int(sEnd), "del")
		}
	}

	// Fall back to genomic-based positions for intronic/UTR/non-coding deletions.
	// Convert extra alt bases to coding strand (only needed for fallback path).
	codingExtraAlt := extraAlt
	if t.IsReverseStrand() && len(extraAlt) > 0 {
		codingExtraAlt = ReverseComplement(extraAlt)
	}

	// A deletion running from an exon into the following intron can still be
	// 3'-shifted, using the stored intron flank (see straddlingDeletionShift).
	if len(extraAlt) == 0 {
		// The genomic shift applies the HGVS rule directly, reading each base
		// through codingBaseAt, so it is asked first. The three below each
		// describe one arrangement of exon and intron and are approximations
		// of the same rule; they answer only where the genomic one has no
		// sequence to compare and returns nothing.
		k := genomicDeletionShift(delStartGenomic, delEndGenomic, t)
		if k == 0 {
			k, _ = straddlingDeletionShift(v, t)
		}
		if k == 0 {
			// The donor flank did not reach; the next exon may be nearer.
			k, _ = acceptorSideDeletionShift(v, t)
		}
		if k == 0 {
			// Or the deletion starts upstream of the start codon.
			k = utr5DeletionShift(v, t)
		}
		if k > 0 {
			d := int64(k)
			if t.IsReverseStrand() {
				d = -d
			}
			delStartGenomic += d
			delEndGenomic += d
			result.HGVSOffset = k
		}
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

// insertionShiftIntronToExon walks a pure insertion sitting inside an intron 3'
// into the following exon, using the acceptor flank to cross the gap.
//
// checkSpliceJunctionDup below handles the single-base case where one flank is
// already exonic. This is the general one: POLD1 inserts a G between the last
// two bases of an intron whose exon then opens GGGGGG, and VEP reports the
// result as c.2959dup, six bases inside the exon.
//
// Positions are addressed on one axis where the exon's first coding base is E
// and E-k is the kth intron base before it, so the walk does not care which
// side of the junction it is on.
func insertionShiftIntronToExon(v *vcf.Variant, t *cache.Transcript, insertedSeq string) (string, int, bool) {
	if insertedSeq == "" || len(t.CDSSequence) == 0 || !t.IsProteinCoding() {
		return "", 0, false
	}
	insAfter := v.Pos // base the insertion follows, in coding order
	if t.IsReverseStrand() {
		insAfter = v.Pos + 1
	}
	if GenomicToCDS(insAfter, t) >= 1 {
		return "", 0, false // already exonic; the CDS path handles it
	}

	// The exon the insertion is moving toward.
	idx := t.FindNearestExonIdx(insAfter)
	if idx < 0 {
		return "", 0, false
	}
	var exon *cache.Exon
	for _, i := range [3]int{idx - 1, idx, idx + 1} {
		if i < 0 || i >= len(t.Exons) {
			continue
		}
		e := &t.Exons[i]
		if !e.IsCoding() {
			continue
		}
		if t.IsReverseStrand() && e.End < insAfter {
			exon = e
			break
		}
		if !t.IsReverseStrand() && e.Start > insAfter {
			exon = e
			break
		}
	}
	if exon == nil || exon.IntronBeforeLen == 0 {
		return "", 0, false
	}
	firstCoding := exon.Start
	if t.IsReverseStrand() {
		firstCoding = exon.End
	}
	eCDS := GenomicToCDS(firstCoding, t)
	if eCDS < 1 {
		return "", 0, false
	}
	// How far before the exon the insertion currently sits.
	gap := firstCoding - insAfter
	if t.IsReverseStrand() {
		gap = insAfter - firstCoding
	}
	if gap < 1 || gap > int64(exon.IntronBeforeLen) {
		return "", 0, false
	}

	f := int64(exon.IntronBeforeLen)
	cds := t.CDSSequence
	baseAt := func(p int64) (byte, bool) {
		if p >= eCDS {
			if int(p) > len(cds) {
				return 0, false
			}
			return cds[p-1], true
		}
		return exon.IntronBeforeBase(int(f - (eCDS - p)))
	}

	L := len(insertedSeq)
	a := eCDS - gap // axis position the insertion follows
	k := 0
	for {
		entering, ok := baseAt(a + int64(k) + 1)
		if !ok || entering != insertedSeq[k%L] {
			break
		}
		k++
	}
	if k == 0 || a+int64(k) < eCDS {
		return "", 0, false // never made it into the exon
	}
	end := a + int64(k)

	// Duplication when the L bases it now follows repeat the insert.
	rotated := insertedSeq[k%L:] + insertedSeq[:k%L]
	start := end - int64(L) + 1
	match := start >= 1
	for i := 0; match && i < L; i++ {
		b, ok := baseAt(start + int64(i))
		if !ok || b != rotated[i] {
			match = false
		}
	}
	if match && start >= eCDS {
		if L == 1 {
			return "c." + strconv.FormatInt(cdsDisplayPos(end, t), 10) + "dup", k, true
		}
		return "c." + strconv.FormatInt(cdsDisplayPos(start, t), 10) + "_" + strconv.FormatInt(cdsDisplayPos(end, t), 10) + "dup", k, true
	}
	return "c." + strconv.FormatInt(cdsDisplayPos(end, t), 10) + "_" + strconv.FormatInt(cdsDisplayPos(end+1, t), 10) + "ins" + rotated, k, true
}

// insertionShiftIntoIntron continues an insertion's 3' shift past the end of an
// exon using the stored intron flank.
//
// shiftInsertionBuf stops at the exon edge, but HGVS does not: an insertion
// sitting on the last base of an exon keeps moving while the intron repeats it.
// KDM6A c.225dup is c.225+1dup to VEP, the duplicated base being the first base
// of the intron rather than the last of the exon.
//
// seq is the inserted sequence on the coding strand and shiftedIdx the 0-based
// CDS index of the base the insertion currently follows. Returns how far it
// moved, the rotated sequence, and the exon's last CDS position.
func insertionShiftIntoIntron(seq []byte, shiftedIdx int, t *cache.Transcript) (k int, rotated []byte, exonEndCDS int64, ok bool) {
	if len(seq) == 0 || shiftedIdx < 0 {
		return 0, nil, 0, false
	}
	g := CDSToGenomic(int64(shiftedIdx+1), t)
	if g == 0 {
		return 0, nil, 0, false
	}
	exon := t.FindExon(g)
	if exon == nil {
		return 0, nil, 0, false
	}
	exonCodingEnd := exon.End
	if t.IsReverseStrand() {
		exonCodingEnd = exon.Start
	}
	eEnd := GenomicToCDS(exonCodingEnd, t)
	if eEnd < 1 || int64(shiftedIdx+1) != eEnd {
		return 0, nil, 0, false // not sitting on the exon's last base
	}

	L := len(seq)
	n := 0
	for {
		entering, okB := exon.IntronAfterBase(n)
		if !okB || entering != seq[n%L] {
			break
		}
		n++
	}
	if n == 0 {
		return 0, nil, 0, false
	}
	r := make([]byte, L)
	for i := 0; i < L; i++ {
		r[i] = seq[(i+n%L)%L]
	}
	return n, r, eEnd, true
}

// formatInsertionInIntron renders an insertion that has moved k bases past the
// end of an exon, as a duplication when the bases it now follows repeat it,
// otherwise as a plain insertion between two intronic positions.
//
// The duplicated stretch can straddle the boundary: RPS6KB2 is c.777_798+12dup,
// starting inside the exon and ending twelve bases into the intron. Positions
// are therefore addressed on a single axis where exonEndCDS+j is the jth intron
// base, and rendered back into c. notation at the end.
func formatInsertionInIntron(rotated []byte, exon *cache.Exon, t *cache.Transcript, exonEndCDS int64, k int) string {
	L := len(rotated)
	baseAt := func(p int64) (byte, bool) {
		if p <= exonEndCDS {
			if p < 1 || int(p) > len(t.CDSSequence) {
				return 0, false
			}
			return t.CDSSequence[p-1], true
		}
		return exon.IntronAfterBase(int(p - exonEndCDS - 1))
	}
	pos := func(p int64) string {
		if p <= exonEndCDS {
			return strconv.FormatInt(cdsDisplayPos(p, t), 10)
		}
		return strconv.FormatInt(cdsDisplayPos(exonEndCDS, t), 10) + "+" + strconv.FormatInt(p-exonEndCDS, 10)
	}

	end := exonEndCDS + int64(k)
	start := end - int64(L) + 1
	match := start >= 1
	for i := 0; match && i < L; i++ {
		b, ok := baseAt(start + int64(i))
		if !ok || b != rotated[i] {
			match = false
		}
	}
	if match {
		if L == 1 {
			return "c." + pos(end) + "dup"
		}
		return "c." + pos(start) + "_" + pos(end) + "dup"
	}
	return "c." + pos(end) + "_" + pos(end+1) + "ins" + string(rotated)
}

// utr5DeletionShift handles a pure deletion that starts upstream of the start
// codon and ends inside the CDS.
//
// The base leaving the front sits in the 5'UTR, which the cache now keeps a
// tail of, so the shift is computable where before it was not:
//
//	SMAD4  c.-28_12del -> c.-27_13del
//	MYCN   c.-35_33del -> c.-34_34del
//
// c. numbering has no zero, so positions step from -1 straight to 1.
// Returns how far the deletion may move 3'.
func utr5DeletionShift(v *vcf.Variant, t *cache.Transcript) int {
	refLen := len(v.Ref)
	if refLen == 0 || len(v.Alt) != 0 || len(t.CDSSequence) == 0 || t.UTR5Len == 0 {
		return 0
	}
	lo, hi := v.Pos, v.Pos+int64(refLen)-1
	codingStart, codingEnd := lo, hi
	if t.IsReverseStrand() {
		codingStart, codingEnd = hi, lo
	}
	if GenomicToCDS(codingStart, t) >= 1 {
		return 0 // starts inside the CDS; a different shape
	}
	endCDS := GenomicToCDS(codingEnd, t)
	if endCDS < 1 {
		return 0 // never reaches the CDS
	}

	// The 5'UTR tail is contiguous with the CDS, so this only applies when no
	// intron sits between the deletion's start and the start codon -- i.e. the
	// start and the first coding base share an exon. Without this check a
	// deletion beginning in an INTRON walks the same way, counts intronic bases
	// as 5'UTR and shifts against sequence that is not there, which regressed
	// 83 rows of HGVSc.
	startExon := t.FindExon(codingStart)
	if startExon == nil {
		return 0
	}
	firstCoding := codingStart
	for i := int64(0); i <= int64(t.UTR5Len); i++ {
		q := firstCoding
		if t.IsReverseStrand() {
			q = codingStart - i
		} else {
			q = codingStart + i
		}
		if GenomicToCDS(q, t) >= 1 {
			firstCoding = q
			break
		}
	}
	if e := t.FindExon(firstCoding); e == nil || e != startExon {
		return 0
	}

	// How far upstream of the CDS the deletion begins, by walking in.
	var upstream int64
	p := codingStart
	for upstream <= int64(t.UTR5Len) {
		if GenomicToCDS(p, t) >= 1 {
			break
		}
		upstream++
		if t.IsReverseStrand() {
			p--
		} else {
			p++
		}
	}
	if upstream < 1 || upstream > int64(t.UTR5Len) {
		return 0
	}

	cds := t.CDSSequence
	baseAt := func(c int64) (byte, bool) {
		if c > 0 {
			if int(c) > len(cds) {
				return 0, false
			}
			return cds[c-1], true
		}
		return t.UTR5BaseAt(int(-c))
	}
	next := func(c int64) int64 {
		if c == -1 {
			return 1
		}
		return c + 1
	}

	front, back := -upstream, endCDS
	k := 0
	for {
		leaving, ok1 := baseAt(front)
		entering, ok2 := baseAt(next(back))
		if !ok1 || !ok2 || leaving != entering {
			break
		}
		front, back = next(front), next(back)
		k++
	}
	return k
}

// acceptorSideDeletionShift handles the same shape as straddlingDeletionShift --
// a deletion running from an exon into the following intron -- but reads the
// flank of the NEXT exon rather than the previous one.
//
// HGVS numbers an intronic position from whichever exon is nearer, so a
// deletion written c.715_863-229del ends 229 bases before the next exon while
// sitting far further than that beyond the previous one. Walking the donor
// flank would mean crossing most of the intron; the acceptor flank reaches it
// in 229. Both flanks are stored, and this is the one that was going unread.
//
// Returns the shift distance and the CDS position the deletion starts at after
// shifting, or zeros when the shape does not apply or the flank is too short.
func acceptorSideDeletionShift(v *vcf.Variant, t *cache.Transcript) (shift int, newStartCDS int64) {
	refLen := len(v.Ref)
	if refLen == 0 || len(v.Alt) != 0 || len(t.CDSSequence) == 0 || !t.IsProteinCoding() {
		return 0, 0
	}
	lo, hi := v.Pos, v.Pos+int64(refLen)-1
	codingStart, codingEnd := lo, hi
	if t.IsReverseStrand() {
		codingStart, codingEnd = hi, lo
	}
	startCDS := GenomicToCDS(codingStart, t)
	if startCDS < 1 || GenomicToCDS(codingEnd, t) > 0 {
		return 0, 0
	}
	idx := t.FindExonIdx(codingStart)
	if idx < 0 {
		return 0, 0
	}
	// Exons are sorted ascending by genomic start, so transcript order runs
	// forward on the plus strand and backward on the minus.
	next := idx + 1
	if t.IsReverseStrand() {
		next = idx - 1
	}
	if next < 0 || next >= len(t.Exons) {
		return 0, 0
	}
	nextExon := &t.Exons[next]

	// How many intronic bases lie between the deletion's 3' end and that exon.
	var gap int64
	if t.IsReverseStrand() {
		gap = codingEnd - nextExon.End
	} else {
		gap = nextExon.Start - codingEnd
	}
	if gap < 1 {
		return 0, 0
	}
	f := int64(nextExon.IntronBeforeLen)
	if f == 0 || gap > f {
		return 0, 0 // the flank does not reach back that far
	}

	cds := t.CDSSequence
	k := 0
	for {
		// Acceptor-relative position of the base entering at the back, counted
		// as a positive distance before the exon.
		rel := gap - 1 - int64(k)
		if rel < 1 {
			break // reached the exon boundary itself
		}
		entering, ok := nextExon.IntronBeforeBase(int(f - rel))
		ci := int(startCDS) - 1 + k
		if !ok || ci >= len(cds) || entering != cds[ci] {
			break
		}
		k++
	}
	if k == 0 {
		return 0, 0
	}
	return k, startCDS + int64(k)
}

// deletionShiftIntoIntron handles a pure deletion whose 3' end reaches the end
// of an exon and can keep shifting 3' into the following intron.
//
// HGVS places a deletion as far 3' as possible and does not stop at a splice
// boundary. CYLD c.2099del sits on the last base of an exon whose intron opens
// with the same base, so VEP reports c.2099+1del: the deletion has moved onto
// the essential donor base. vibe-vep stopped at the exon end and reported the
// unshifted c.2099del.
//
// Returns the post-shift exonic start (1-based CDS), the exon's last CDS
// position, the deletion length, and how many bases it moved into the intron.
// ok is false when the shape does not apply or no flank is stored.
func deletionShiftIntoIntron(v *vcf.Variant, t *cache.Transcript) (startCDS, exonEndCDS int64, delLen, k int, ok bool) {
	refLen := len(v.Ref)
	if refLen == 0 || len(v.Alt) != 0 || len(t.CDSSequence) == 0 || !t.IsProteinCoding() {
		return 0, 0, 0, 0, false
	}
	lo, hi := v.Pos, v.Pos+int64(refLen)-1
	codingStart, codingEnd := lo, hi
	if t.IsReverseStrand() {
		codingStart, codingEnd = hi, lo
	}
	sCDS := GenomicToCDS(codingStart, t)
	eCDS := GenomicToCDS(codingEnd, t)
	if sCDS < 1 || eCDS < 1 {
		return 0, 0, 0, 0, false // not fully exonic
	}
	exon := t.FindExon(codingEnd)
	if exon == nil {
		return 0, 0, 0, 0, false
	}
	exonCodingEnd := exon.End
	if t.IsReverseStrand() {
		exonCodingEnd = exon.Start
	}
	eEnd := GenomicToCDS(exonCodingEnd, t)
	if eEnd < 1 {
		return 0, 0, 0, 0, false
	}

	// Shift within the exon first, exactly as the ordinary path does.
	sStart, sEnd := shiftDeletionThreePrime(int(sCDS-1), int(eCDS-1), t.CDSSequence, int(eEnd-1))
	if int64(sEnd) != eEnd-1 {
		return 0, 0, 0, 0, false // 3' end never reaches the exon edge; nothing to cross
	}
	s := int64(sStart) + 1
	L := sEnd - sStart + 1

	// Keep going into the intron. Position p is CDS while p <= eEnd, and the
	// (p-eEnd)th intron base beyond that, so one accessor covers both sides.
	cds := t.CDSSequence
	baseAt := func(p int64) (byte, bool) {
		if p <= eEnd {
			if p < 1 || int(p) > len(cds) {
				return 0, false
			}
			return cds[p-1], true
		}
		return exon.IntronAfterBase(int(p - eEnd - 1))
	}
	n := 0
	for {
		leaving, okL := baseAt(s + int64(n))
		entering, okE := baseAt(eEnd + int64(n) + 1)
		if !okL || !okE || leaving != entering {
			break
		}
		n++
	}
	if n == 0 {
		return 0, 0, 0, 0, false
	}
	return s + int64(n), eEnd, L, n, true
}

// formatDeletionAcrossDonor renders a deletion that has shifted k bases past the
// end of an exon, keeping whatever exonic part remains.
func formatDeletionAcrossDonor(t *cache.Transcript, startCDS, exonEndCDS int64, k int) string {
	if startCDS <= exonEndCDS {
		// Part of the deletion is still exonic: c.START_EXONEND+k del
		return "c." + strconv.FormatInt(cdsDisplayPos(startCDS, t), 10) + "_" +
			strconv.FormatInt(cdsDisplayPos(exonEndCDS, t), 10) + "+" + strconv.Itoa(k) + "del"
	}
	// Entirely intronic now: c.EXONEND+a_EXONEND+k del
	a := int(startCDS - exonEndCDS)
	base := "c." + strconv.FormatInt(cdsDisplayPos(exonEndCDS, t), 10) + "+"
	if a == k {
		return base + strconv.Itoa(k) + "del"
	}
	return base + strconv.Itoa(a) + "_" + base + strconv.Itoa(k) + "del"
}

// deletionHitsDonor reports whether the shifted deletion covers the first two
// intron bases, the essential donor dinucleotide.
func deletionHitsDonor(delLen, k int) bool {
	first := k - delLen + 1
	if first < 1 {
		first = 1
	}
	return first <= 2 && k >= 1
}

// straddlingDeletionShift computes how far a pure deletion running from an exon
// into the following intron may be shifted 3'.
//
// HGVS places a deletion as far 3' as possible. Each step is allowed when the
// base leaving the front of the deleted span equals the base entering at the
// back; here the front is exonic and the back is intronic, so the test needs
// intron sequence. That is what Exon.IntronAfterPacked carries.
//
//	MET  vibe c.3015_3028+6del -> VEP c.3016_3028+7del  (offset 1)
//	B2M  vibe c.59_67+237del   -> VEP c.61_67+239del    (offset 2)
//
// Returns the shift distance and the CDS position the deletion starts at after
// shifting. Returns 0 when the shape does not apply or no flank is stored, in
// which case the caller keeps the unshifted description.
func straddlingDeletionShift(v *vcf.Variant, t *cache.Transcript) (shift int, newStartCDS int64) {
	refLen := len(v.Ref)
	if refLen == 0 || len(v.Alt) != 0 || len(t.CDSSequence) == 0 || !t.IsProteinCoding() {
		return 0, 0
	}
	lo, hi := v.Pos, v.Pos+int64(refLen)-1

	// Orient the span: codingStart is the 5' end of the deletion on the coding
	// strand, codingEnd the 3' end.
	codingStart, codingEnd := lo, hi
	if t.IsReverseStrand() {
		codingStart, codingEnd = hi, lo
	}
	startCDS := GenomicToCDS(codingStart, t)
	if startCDS < 1 || GenomicToCDS(codingEnd, t) > 0 {
		return 0, 0 // not a straddle: either the 5' end is intronic or the 3' end is still coding
	}
	exon := t.FindExon(codingStart)
	if exon == nil {
		return 0, 0
	}

	// How many intronic bases the deletion already covers.
	exonCodingEnd := exon.End
	if t.IsReverseStrand() {
		exonCodingEnd = exon.Start
	}
	tail := codingEnd - exonCodingEnd
	if t.IsReverseStrand() {
		tail = exonCodingEnd - codingEnd
	}
	if tail < 1 {
		return 0, 0
	}

	cds := t.CDSSequence
	k := 0
	for {
		fi := int(tail) + k // 0-based index of the base just past the deletion
		ci := int(startCDS) - 1 + k
		entering, ok := exon.IntronAfterBase(fi)
		if !ok || ci >= len(cds) || entering != cds[ci] {
			break
		}
		k++
	}
	if k == 0 {
		return 0, 0
	}
	return k, startCDS + int64(k)
}

// shiftedIntronicDelCDS handles a pure deletion lying entirely within an intron
// immediately 5' (on the coding strand) of a coding exon, where the deleted
// bases keep repeating into that exon. HGVS assigns the most 3' position
// possible, so such a deletion slides across the splice junction: POLD1 deletes
// a G from the last intron base of an exon that opens GGGGGG, and VEP reports
// c.2959del with hgvs_offset 6, not c.2954-1del.
//
// No intron sequence is needed, which matters because vibe-vep has none. The
// deleted bases come from the variant itself and the run continues in
// CDSSequence, so sliding k bases is allowed while the exon base at j matches
// deleted[j%L]. Once k reaches L the deletion lies wholly inside the exon;
// below that it still straddles the junction and HGVS keeps the intronic form.
//
// Returns the 1-based CDS range it lands on and the shift distance, or zeros.
func shiftedIntronicDelCDS(v *vcf.Variant, t *cache.Transcript) (cdsStart, cdsEnd int64, shift int) {
	refLen := len(v.Ref)
	if refLen == 0 || len(v.Alt) != 0 || len(t.CDSSequence) == 0 || !t.IsProteinCoding() {
		return 0, 0, 0
	}
	delStartGenomic := v.Pos
	delEndGenomic := v.Pos + int64(refLen) - 1
	if GenomicToCDS(delStartGenomic, t) > 0 || GenomicToCDS(delEndGenomic, t) > 0 {
		return 0, 0, 0 // touches coding sequence; the exonic path handles it
	}

	delSeq := v.Ref
	next := delEndGenomic + 1
	if t.IsReverseStrand() {
		delSeq = ReverseComplement(delSeq)
		next = delStartGenomic - 1
	}
	exonStart := GenomicToCDS(next, t)
	if exonStart < 1 {
		return 0, 0, 0
	}

	k := 0
	for {
		idx := int(exonStart) - 1 + k
		if idx >= len(t.CDSSequence) || t.CDSSequence[idx] != delSeq[k%refLen] {
			break
		}
		k++
	}
	if k < refLen {
		return 0, 0, 0
	}
	cdsStart = exonStart + int64(k-refLen)
	return cdsStart, cdsStart + int64(refLen) - 1, k
}

// cdsDisplayPos turns a CDS position into the number an HGVSc string shows.
// Internally a CDS position counts from the notional complete first codon, so
// on a 5'-incomplete (cds_start_NF) transcript it includes the bases GENCODE
// never annotated. VEP numbers c. from the first base that is actually there,
// so those padded bases come back off before the number is rendered. For every
// transcript with a complete start -- all but ~5% of GENCODE -- this is a no-op.
func cdsDisplayPos(cdsPos int64, t *cache.Transcript) int64 {
	return cdsPos - int64(t.CDSStartOffset)
}

// cdsPosRangeStr builds a CDS position string like "c.34del", "c.34_36del",
// "c.35dup", etc. using a stack-allocated buffer for zero intermediate allocations.
func cdsPosRangeStr(t *cache.Transcript, start, end int, suffix string) string {
	var buf [64]byte
	n := copy(buf[:], "c.")
	n += putInt64(buf[n:], cdsDisplayPos(int64(start), t))
	if start != end {
		buf[n] = '_'
		n++
		n += putInt64(buf[n:], cdsDisplayPos(int64(end), t))
	}
	n += copy(buf[n:], suffix)
	return string(buf[:n])
}

// shiftInsertionBuf shifts an insertion rightward (3' direction) in CDS space.
// It operates on the provided mutable byte slice, avoiding allocations.
// The shift stops at maxIdx (0-based) to avoid crossing exon-exon junctions.
// Returns the new anchor index.
func shiftInsertionBuf(seq []byte, cdsAnchorIdx int, cdsSeq string, maxIdx int) int {
	idx := cdsAnchorIdx
	for idx < maxIdx && idx+1 < len(cdsSeq) && cdsSeq[idx+1] == seq[0] {
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
			return strconv.FormatInt(cdsDisplayPos(cdsPos, t), 10)
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
			return strconv.FormatInt(cdsDisplayPos(cdsPos, t), 10)
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

// cdsExonEndIdx returns the 0-based CDS index of the last base in the same exon
// as the given 0-based CDS index. This is used to prevent 3' shifting from crossing
// exon-exon junctions per HGVS rules. Returns len(cdsSeq)-1 if no exon boundary
// information is available.
func cdsExonEndIdx(cdsIdx0 int, t *cache.Transcript) int {
	if len(t.CDSRegions) == 0 {
		return len(t.CDSSequence) - 1
	}
	for _, r := range t.CDSRegions {
		regionLen := int(r.GenomicEnd - r.GenomicStart + 1)
		regionStart := int(r.CDSOffset)
		regionEnd := regionStart + regionLen - 1
		if cdsIdx0 >= regionStart && cdsIdx0 <= regionEnd {
			return regionEnd
		}
	}
	return len(t.CDSSequence) - 1
}

// shiftInsertionThreePrime shifts an insertion rightward (3' direction) in CDS space.
// cdsAnchorIdx is the 0-based CDS index of the VCF anchor base (insertion is after this position).
// maxIdx limits the shift to avoid crossing exon-exon junctions.
// Returns the shifted inserted sequence and new anchor index.
func shiftInsertionThreePrime(insertedSeq string, cdsAnchorIdx int, cdsSeq string, maxIdx int) (string, int) {
	seq := []byte(insertedSeq)
	idx := cdsAnchorIdx
	for idx < maxIdx && idx+1 < len(cdsSeq) && cdsSeq[idx+1] == seq[0] {
		first := seq[0]
		copy(seq, seq[1:])
		seq[len(seq)-1] = first
		idx++
	}
	return string(seq), idx
}

// shiftDeletionThreePrime shifts a deletion rightward (3' direction) in CDS space.
// delStart and delEnd are 0-based inclusive CDS indices of the deleted bases.
// The shift stops at maxIdx (0-based) to avoid crossing exon-exon junctions.
// Returns the shifted start and end indices.
func shiftDeletionThreePrime(delStart, delEnd int, cdsSeq string, maxIdx int) (int, int) {
	for delEnd < maxIdx && delEnd+1 < len(cdsSeq) && cdsSeq[delStart] == cdsSeq[delEnd+1] {
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
