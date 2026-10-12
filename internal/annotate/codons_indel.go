package annotate

import (
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// formatCodonChangeCDS renders the Codons string Ensembl VEP reports for an
// indel or multi-nucleotide variant, working in CDS coordinates.
//
// VEP expands the affected range out to whole-codon boundaries, then writes
// ref/alt over that window with the bases that actually changed in uppercase
// and the untouched flanking bases in lowercase. Unlike a substitution the two
// sides differ in length:
//
//	c.5946del          agT/ag           (one base deleted from the codon end)
//	c.1856del          cCt/ct           (deleted base mid-codon)
//	c.4049_4051del     gAGGcc/gcc       (three bases spanning two codons)
//	c.29_30insCAG      cag/caCAGg       (insertion inside a codon)
//	c.627delinsAT      atG/atAT         (one base replaced by two)
//	c.1132_1133delinsT AGt/Tt
//
// delStart/delEnd are 1-based inclusive CDS positions of the deleted bases; an
// empty range (delEnd < delStart) means a pure insertion after delEnd. inserted
// holds the inserted bases on the coding strand, empty for a pure deletion.
//
// An insertion that lands exactly on a codon boundary disturbs no existing
// codon, and VEP reports it as "-/INSERTED".
func formatCodonChangeCDS(cdsSeq string, delStart, delEnd int64, inserted string) string {
	n := int64(len(cdsSeq))
	if n == 0 {
		return ""
	}
	pureInsertion := delEnd < delStart
	if pureInsertion && delEnd%3 == 0 {
		// Codon-boundary insertion: no existing codon is disturbed.
		var buf [256]byte
		if len(inserted)+2 > len(buf) {
			return "-/" + upper(inserted)
		}
		k := copy(buf[:], "-/")
		k += copyUpper(buf[k:], inserted)
		return string(buf[:k])
	}

	first, last := delStart, delEnd
	if pureInsertion {
		first, last = delEnd, delEnd
	}
	if first < 1 || last > n {
		return ""
	}

	codonStart := ((first-1)/3)*3 + 1
	codonEnd := ((last-1)/3 + 1) * 3
	if codonEnd > n {
		codonEnd = n
	}
	if codonStart < 1 || codonStart > codonEnd {
		return ""
	}

	// Worst case: the whole window twice, plus the insertion twice, plus '/'.
	need := int(codonEnd-codonStart+1)*2 + len(inserted) + 1
	var stack [256]byte
	buf := stack[:]
	if need > len(buf) {
		buf = make([]byte, need)
	}
	k := 0

	// ref side: deleted bases uppercase, flanks lowercase
	for pos := codonStart; pos <= codonEnd; pos++ {
		b := cdsSeq[pos-1]
		if !pureInsertion && pos >= delStart && pos <= delEnd {
			buf[k] = upperByte(b)
		} else {
			buf[k] = lowerByte(b)
		}
		k++
	}
	buf[k] = '/'
	k++
	altStart := k

	// alt side: inserted bases uppercase, retained flanks lowercase
	for pos := codonStart; pos <= codonEnd; pos++ {
		if !pureInsertion && pos == delStart && inserted != "" {
			k += copyUpper(buf[k:], inserted)
		}
		if !pureInsertion && pos >= delStart && pos <= delEnd {
			continue
		}
		buf[k] = lowerByte(cdsSeq[pos-1])
		k++
		if pureInsertion && pos == delEnd {
			k += copyUpper(buf[k:], inserted)
		}
	}
	// A window that is wholly deleted leaves the alt side empty; VEP writes
	// that as "-", e.g. "GGC/-".
	if k == altStart {
		buf[k] = '-'
		k++
	}
	return string(buf[:k])
}

func formatCodonChangeMNV(cdsSeq string, start, end int64, altBases string) string {
	n := int64(len(cdsSeq))
	if n == 0 || start < 1 || end > n || end < start {
		return ""
	}
	if int64(len(altBases)) != end-start+1 {
		return ""
	}
	codonStart := ((start-1)/3)*3 + 1
	codonEnd := ((end-1)/3 + 1) * 3
	if codonEnd > n {
		codonEnd = n
	}
	width := int(codonEnd - codonStart + 1)
	var stack [256]byte
	buf := stack[:]
	if width*2+1 > len(buf) {
		buf = make([]byte, width*2+1)
	}
	k := 0
	for pos := codonStart; pos <= codonEnd; pos++ {
		b := cdsSeq[pos-1]
		if pos >= start && pos <= end {
			buf[k] = upperByte(b)
		} else {
			buf[k] = lowerByte(b)
		}
		k++
	}
	buf[k] = '/'
	k++
	for pos := codonStart; pos <= codonEnd; pos++ {
		if pos >= start && pos <= end {
			buf[k] = upperByte(altBases[pos-start])
		} else {
			buf[k] = lowerByte(cdsSeq[pos-1])
		}
		k++
	}
	return string(buf[:k])
}

// copyUpper writes s into dst uppercased and returns the number of bytes written.
func copyUpper(dst []byte, s string) int {
	for i := 0; i < len(s); i++ {
		dst[i] = upperByte(s[i])
	}
	return len(s)
}

func upperByte(b byte) byte { return b &^ 0x20 }
func lowerByte(b byte) byte { return b | 0x20 }

func lower(s string) string {
	b := []byte(s)
	for i := range b {
		b[i] = lowerByte(b[i])
	}
	return string(b)
}

func upper(s string) string {
	b := []byte(s)
	for i := range b {
		b[i] = upperByte(b[i])
	}
	return string(b)
}

// codonsSpellStop reports whether reading the codon string in frame from its
// start hits a stop codon. The Codons field mixes case to mark the changed
// bases and writes "-" for an empty side, neither of which affects the frame.
func codonsSpellStop(codons string) bool {
	seq := make([]byte, 0, len(codons))
	for i := 0; i < len(codons); i++ {
		if c := codons[i]; c != '-' {
			seq = append(seq, upperBase(c))
		}
	}
	for i := 0; i+3 <= len(seq); i += 3 {
		switch string(seq[i : i+3]) {
		case "TAA", "TAG", "TGA":
			return true
		}
	}
	return false
}

func upperBase(c byte) byte {
	if c >= 'a' && c <= 'z' {
		return c - 'a' + 'A'
	}
	return c
}

// lengthChangeTerms are the consequences a coding indel or delins carries.
// The stop_gained test below applies to all of them; substitutions are left
// alone because the single-base path already classifies those.
var lengthChangeTerms = [...]string{
	ConsequenceFrameshiftVariant,
	ConsequenceInframeInsertion,
	ConsequenceInframeDeletion,
	ConsequenceProteinAltering,
}

// addStopGainedFromCodons adds stop_gained to a coding indel whose altered
// codons spell out a stop.
//
// VEP's stop_gained test is on the peptide the affected codons translate to,
// not on where the shifted reading frame eventually terminates, so the two
// come apart in both directions: GATA3 caGTCA/ca is p.S70* yet plain
// frameshift_variant, while gct/gcCTCTGGCCCTGGCTGCTGATGGCt is
// stop_gained,frameshift_variant despite running on to p.T20Sfs*6.
//
// Reading the Codons alt side in frame and requiring the ref side not to
// already hold a stop reproduces VEP exactly, and on every class of variant:
// 530 of 60,028 frameshifts, 978 of 17,207 in-frame indels and MNVs, and
// 44,171 of 377,879 substitutions in a VEP111 MSK-IMPACT MAF, with no false
// positives and none missed anywhere.
//
// The length term stays. VEP writes stop_gained,inframe_insertion and
// stop_gained,frameshift_variant, and even where it drops the term the
// leading stop_gained is what decides Variant_Classification, so keeping it
// costs nothing there.
func addStopGainedFromCodons(result *ConsequenceResult) {
	if result.CodonChange == "" {
		return
	}
	coding := false
	for _, term := range lengthChangeTerms {
		if hasTerm(result.Consequence, term) {
			coding = true
			break
		}
	}
	if !coding {
		return
	}
	// A consequence that already names stop_gained still needs the length-term
	// decision below, so this only skips re-adding the term.
	alreadyStopGained := hasTerm(result.Consequence, ConsequenceStopGained)
	slash := -1
	for i := 0; i < len(result.CodonChange); i++ {
		if result.CodonChange[i] == '/' {
			slash = i
			break
		}
	}
	if slash < 0 {
		return
	}
	ref, alt := result.CodonChange[:slash], result.CodonChange[slash+1:]
	if !alreadyStopGained && (!codonsSpellStop(alt) || codonsSpellStop(ref)) {
		return
	}
	terms := splitConsequenceTerms(result.Consequence)
	if dropsLengthTerm(ref, alt, result.Consequence) {
		kept := terms[:0]
		for _, t := range terms {
			if !isLengthChangeTerm(t) {
				kept = append(kept, t)
			}
		}
		terms = kept
	}
	result.Consequence = joinConsequenceTerms(append([]string{ConsequenceStopGained}, terms...))
	result.Impact = GetImpact(firstTerm(result.Consequence))
}

func isLengthChangeTerm(term string) bool {
	for _, t := range lengthChangeTerms {
		if t == term {
			return true
		}
	}
	return false
}

// dropsLengthTerm reports whether VEP reports stop_gained on its own rather
// than alongside the length change.
//
// It does so when the variant replaces a reference codon and the codon that
// takes its place is itself a stop: nothing survives at that position, so
// there is no insertion or replacement left to describe. CCGGga/Tga is
// stop_gained alone, while cgCTCCAACGag/cgTTag keeps
// protein_altering_variant because its first altered codon still spells Arg.
//
// Two shapes always keep the term. A pure deletion does - all 48 in a VEP111
// MSK-IMPACT MAF, even where the junction it creates spells a stop - and so
// does an insertion at a codon boundary, which replaces no reference codon at
// all. Frameshifts are never affected: VEP writes stop_gained,
// frameshift_variant and never stop_gained alone for one.
//
// Together these reproduce VEP's choice on all 122 non-frameshift coding
// indels that gain a stop, with no exceptions.
func dropsLengthTerm(ref, alt, consequence string) bool {
	if hasTerm(consequence, ConsequenceFrameshiftVariant) {
		return false
	}
	if ref == "" || ref == "-" {
		return false // whole codons inserted at a boundary
	}
	if !hasUpper(alt) {
		return false // pure deletion
	}
	return firstCodonIsStop(alt)
}

func hasUpper(s string) bool {
	for i := 0; i < len(s); i++ {
		if s[i] >= 'A' && s[i] <= 'Z' {
			return true
		}
	}
	return false
}

// firstCodonIsStop reports whether the first complete codon of a Codons side
// is a stop.
func firstCodonIsStop(codons string) bool {
	var c [3]byte
	n := 0
	for i := 0; i < len(codons) && n < 3; i++ {
		if codons[i] != '-' {
			c[n] = upperBase(codons[i])
			n++
		}
	}
	if n < 3 {
		return false
	}
	switch string(c[:]) {
	case "TAA", "TAG", "TGA":
		return true
	}
	return false
}

// splitConsequenceTerms splits a rendered consequence back into its terms.
func splitConsequenceTerms(consequence string) []string {
	var terms []string
	start := 0
	for i := 0; i <= len(consequence); i++ {
		if i == len(consequence) || consequence[i] == ',' {
			if i > start {
				terms = append(terms, consequence[start:i])
			}
			start = i + 1
		}
	}
	return terms
}

// hasTerm reports whether a rendered consequence contains the given term.
func hasTerm(consequence, term string) bool {
	for _, t := range splitConsequenceTerms(consequence) {
		if t == term {
			return true
		}
	}
	return false
}

// translateCodonsInto writes the peptide a Codons side spells into buf and
// returns its length, ignoring the mixed case that marks the changed bases and
// the "-" that marks an empty side. A trailing partial codon is dropped.
func translateCodonsInto(buf []byte, codons string) int {
	var seq [96]byte
	n := 0
	for i := 0; i < len(codons) && n < len(seq); i++ {
		if c := codons[i]; c != '-' {
			seq[n] = upperBase(c)
			n++
		}
	}
	out := 0
	for i := 0; i+3 <= n && out < len(buf); i += 3 {
		buf[out] = TranslateCodon(string(seq[i : i+3]))
		out++
	}
	return out
}

// relabelInsertionSplittingACodon swaps inframe_insertion for
// protein_altering_variant when the inserted bases break a codon apart rather
// than adding whole ones.
//
// EGFR c.2308_2309insGTT inserts three bases one base into the codon, so the
// aspartate it lands in is replaced by a glycine and a valine - a delins at
// the protein level, which VEP reports as protein_altering_variant.
//
// Two shapes keep the plain insertion term: bases added at a codon boundary,
// where the Codons reference side is empty, and bases that split a codon but
// leave the original residue standing at one end of the result. Together with
// this relabelling they account for all 1,698 in-frame insertions in a VEP111
// MSK-IMPACT MAF, with no exceptions either way.
func relabelInsertionSplittingACodon(result *ConsequenceResult) {
	if !hasTerm(result.Consequence, ConsequenceInframeInsertion) {
		return
	}
	slash := -1
	for i := 0; i < len(result.CodonChange); i++ {
		if result.CodonChange[i] == '/' {
			slash = i
			break
		}
	}
	if slash < 0 {
		return
	}
	ref, alt := result.CodonChange[:slash], result.CodonChange[slash+1:]
	var refPep, altPep [32]byte
	rn := translateCodonsInto(refPep[:], ref)
	an := translateCodonsInto(altPep[:], alt)
	if rn == 0 || an < rn {
		// Nothing was replaced: whole codons were added at a boundary.
		return
	}
	r := string(refPep[:rn])
	if string(altPep[:rn]) == r || string(altPep[an-rn:an]) == r {
		// The original residues survive at one end, so this reads as an
		// insertion into the protein rather than a replacement.
		return
	}
	result.Consequence = relabelInframeAsProteinAltering(result.Consequence)
	result.Impact = GetImpact(firstTerm(result.Consequence))
}

// clearCodonsForSpliceSite drops the Codons string for a variant that destroys
// a splice site.
//
// Once the donor or acceptor is gone the transcript's reading frame is no
// longer something the annotation predicts, and VEP reports no codons for it:
// of the 16,264 variants carrying splice_donor_variant or
// splice_acceptor_variant without stop_lost in a VEP111 MSK-IMPACT MAF, 16,262
// have an empty Codons field.
//
// A variant that also loses the stop codon keeps its codons - 61 of 63 - and
// those are deletions taking a whole terminal exon, where the junction still
// spells something readable.
func clearCodonsForSpliceSite(result *ConsequenceResult) {
	if result.CodonChange == "" {
		return
	}
	if !hasTerm(result.Consequence, ConsequenceSpliceDonor) &&
		!hasTerm(result.Consequence, ConsequenceSpliceAcceptor) {
		return
	}
	if hasTerm(result.Consequence, ConsequenceStopLost) {
		return
	}
	result.CodonChange = ""
}

// relabelStopLostWhenItShiftsClear turns stop_lost into stop_retained_variant
// for a variant whose 3'-shifted representation no longer reaches the stop
// codon.
//
// The stop_lost tests run on the unshifted genomic span, which can overlap the
// stop codon even when the change describes itself entirely past it. ERBB3
// deletes 510 bases starting at the stop codon, but the bases that follow
// spell the same stop, so HGVS writes it c.*1_*510del and the protein still
// terminates where it did. vibe-vep called that stop_lost; VEP reports
// stop_retained_variant,3_prime_UTR_variant.
//
// A stop codon the variant does not touch cannot be lost, so this only ever
// corrects an over-call.
func relabelStopLostWhenItShiftsClear(v *vcf.Variant, t *cache.Transcript, result *ConsequenceResult) {
	if !hasTerm(result.Consequence, ConsequenceStopLost) || len(v.Ref) == 0 {
		return
	}
	// A deletion that breaks a splice site gets stop_lost for a different
	// reason -- VEP reads it off the deleted genomic sequence, not off the
	// real stop codon (see spliceSpanProteinTerms) -- and relabelling that
	// would answer a question this function is not asking. VEP never pairs
	// stop_retained_variant with a splice site term.
	if hasTerm(result.Consequence, ConsequenceSpliceDonor) ||
		hasTerm(result.Consequence, ConsequenceSpliceAcceptor) {
		return
	}
	shift := int64(result.HGVSOffset)
	if t.IsReverseStrand() {
		shift = -shift
	}
	start := v.Pos + shift
	end := start + int64(len(v.Ref)) - 1
	stopStart, stopEnd := t.CDSEnd-2, t.CDSEnd
	if t.IsReverseStrand() {
		stopStart, stopEnd = t.CDSStart, t.CDSStart+2
	}
	if end >= stopStart && start <= stopEnd {
		return // still covers the stop codon
	}
	terms := splitConsequenceTerms(result.Consequence)
	for i, term := range terms {
		if term == ConsequenceStopLost {
			terms[i] = ConsequenceStopRetained
		}
	}
	result.Consequence = joinConsequenceTerms(terms)
	result.Impact = GetImpact(firstTerm(result.Consequence))
}
