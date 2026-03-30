package annotate

import (
	"fmt"
	"regexp"
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// ReverseMapProteinChange maps a protein change (e.g. KRAS G12C) back to
// genomic variant(s) using the canonical transcript's CDS sequence.
func ReverseMapProteinChange(c *cache.Cache, geneName string, refAA byte, protPos int64, altAA byte) ([]*vcf.Variant, error) {
	transcripts := c.FindTranscriptsByGene(geneName)
	if len(transcripts) == 0 {
		return nil, fmt.Errorf("gene %q not found in transcript cache", geneName)
	}

	// Find canonical protein-coding transcript
	var canonical *cache.Transcript
	for _, t := range transcripts {
		if t.IsCanonicalMSK && t.IsProteinCoding() {
			canonical = t
			break
		}
	}
	if canonical == nil {
		// Fall back to first protein-coding transcript
		for _, t := range transcripts {
			if t.IsProteinCoding() {
				canonical = t
				break
			}
		}
	}
	if canonical == nil {
		return nil, fmt.Errorf("no protein-coding transcript found for gene %q", geneName)
	}

	return reverseMapProtein(canonical, refAA, protPos, altAA)
}

// reverseMapProtein maps a protein change on a specific transcript to genomic variant(s).
func reverseMapProtein(t *cache.Transcript, refAA byte, protPos int64, altAA byte) ([]*vcf.Variant, error) {
	if len(t.CDSSequence) == 0 {
		return nil, fmt.Errorf("transcript %s has no CDS sequence", t.ID)
	}

	// Get reference codon
	refCodon := GetCodon(t.CDSSequence, protPos)
	if len(refCodon) != 3 {
		return nil, fmt.Errorf("protein position %d out of range for transcript %s", protPos, t.ID)
	}

	// Verify reference AA matches
	actualRefAA := TranslateCodon(refCodon)
	if actualRefAA != refAA {
		return nil, fmt.Errorf("reference amino acid mismatch at position %d: expected %c, got %c in transcript %s",
			protPos, refAA, actualRefAA, t.ID)
	}

	// CDS start position for this codon (1-based)
	cdsStart := (protPos-1)*3 + 1

	// Enumerate all single-base mutations of the codon that produce altAA
	var variants []*vcf.Variant
	for posInCodon := 0; posInCodon < 3; posInCodon++ {
		for _, base := range "ACGT" {
			if byte(base) == refCodon[posInCodon] {
				continue // skip ref base
			}
			mutCodon := MutateCodon(refCodon, posInCodon, byte(base))
			if TranslateCodon(mutCodon) != altAA {
				continue
			}

			cdsPos := cdsStart + int64(posInCodon)
			genomicPos := CDSToGenomic(cdsPos, t)
			if genomicPos == 0 {
				continue
			}

			refBase := string(refCodon[posInCodon])
			altBase := string(base)

			// Account for reverse strand
			if t.IsReverseStrand() {
				refBase = string(Complement(refCodon[posInCodon]))
				altBase = string(Complement(byte(base)))
			}

			variants = append(variants, &vcf.Variant{
				Chrom: t.Chrom,
				Pos:   genomicPos,
				Ref:   refBase,
				Alt:   altBase,
			})
		}
	}

	if len(variants) == 0 {
		return nil, fmt.Errorf("no single-base mutation of codon %s at position %d produces %c",
			refCodon, protPos, altAA)
	}

	return variants, nil
}

// reCDSChange parses a simple CDS substitution like "35G>T".
var reCDSChange = regexp.MustCompile(`^(\d+)([ACGT])>([ACGT])$`)

// reCDSDeletion parses a CDS deletion like "923del" or "100_102del".
var reCDSDeletion = regexp.MustCompile(`^(\d+)(?:_(\d+))?del$`)

// ReverseMapHGVSc maps an HGVSc notation (e.g. KRAS c.35G>T or KRAS c.34del)
// back to a genomic variant using the transcript's CDS-to-genomic mapping.
func ReverseMapHGVSc(c *cache.Cache, geneOrTranscript string, cdsChange string) ([]*vcf.Variant, error) {
	variants, _, err := ReverseMapHGVScWithWarning(c, geneOrTranscript, cdsChange)
	return variants, err
}

// ReverseMapHGVScWithWarning is like ReverseMapHGVSc but also returns a warning
// string if the exact transcript version was requested but not found.
func ReverseMapHGVScWithWarning(c *cache.Cache, geneOrTranscript string, cdsChange string) ([]*vcf.Variant, string, error) {
	// Check for version mismatch warning.
	result := FindHGVScTranscriptWithWarning(c, geneOrTranscript)
	warning := result.Warning

	// Try substitution first
	if m := reCDSChange.FindStringSubmatch(cdsChange); m != nil {
		variants, err := reverseMapHGVScSubstitution(c, geneOrTranscript, m)
		return variants, warning, err
	}

	// Try deletion
	if m := reCDSDeletion.FindStringSubmatch(cdsChange); m != nil {
		variants, err := reverseMapHGVScDeletion(c, geneOrTranscript, m)
		return variants, warning, err
	}

	return nil, warning, fmt.Errorf("unsupported CDS change notation %q (supported: substitutions like 35G>T, deletions like 923del or 100_102del)", cdsChange)
}

// findHGVScTranscript resolves a gene name or transcript ID to a transcript.
// TranscriptLookupResult holds the found transcript and any warnings.
type TranscriptLookupResult struct {
	Transcript *cache.Transcript
	Warning    string // non-empty if the exact version was not found
}

func findHGVScTranscript(c *cache.Cache, geneOrTranscript string) (*cache.Transcript, error) {
	result := FindHGVScTranscriptWithWarning(c, geneOrTranscript)
	if result.Transcript == nil {
		return nil, fmt.Errorf("transcript %q not found", geneOrTranscript)
	}
	return result.Transcript, nil
}

// FindHGVScTranscriptWithWarning looks up a transcript by ID (with optional version)
// or gene name. Returns a warning if the exact version was requested but not found.
func FindHGVScTranscriptWithWarning(c *cache.Cache, geneOrTranscript string) TranscriptLookupResult {
	if strings.HasPrefix(geneOrTranscript, "ENST") {
		// Try exact match first (includes version).
		transcript := c.GetTranscript(geneOrTranscript)
		if transcript != nil {
			return TranscriptLookupResult{Transcript: transcript}
		}

		// Fallback: match by base ID (strip version).
		transcripts := findTranscriptByPrefix(c, geneOrTranscript)
		if len(transcripts) > 0 {
			transcript = transcripts[0]
		}
		if transcript == nil {
			return TranscriptLookupResult{}
		}

		// Check if a specific version was requested but didn't match.
		var warning string
		if strings.Contains(geneOrTranscript, ".") {
			warning = fmt.Sprintf("requested transcript version %s not found, using %s instead", geneOrTranscript, transcript.ID)
		}
		return TranscriptLookupResult{Transcript: transcript, Warning: warning}
	}

	transcripts := c.FindTranscriptsByGene(geneOrTranscript)
	if len(transcripts) == 0 {
		return TranscriptLookupResult{}
	}
	for _, t := range transcripts {
		if t.IsCanonicalMSK && t.IsProteinCoding() {
			return TranscriptLookupResult{Transcript: t}
		}
	}
	for _, t := range transcripts {
		if t.IsProteinCoding() {
			return TranscriptLookupResult{Transcript: t}
		}
	}
	return TranscriptLookupResult{}
}

// reverseMapHGVScSubstitution handles simple CDS substitutions like "35G>T".
func reverseMapHGVScSubstitution(c *cache.Cache, geneOrTranscript string, m []string) ([]*vcf.Variant, error) {
	cdsPos, _ := strconv.ParseInt(m[1], 10, 64)
	cdsRef := m[2]
	cdsAlt := m[3]

	transcript, err := findHGVScTranscript(c, geneOrTranscript)
	if err != nil {
		return nil, err
	}

	// Verify CDS ref base
	if len(transcript.CDSSequence) > 0 && cdsPos <= int64(len(transcript.CDSSequence)) {
		actualRef := string(transcript.CDSSequence[cdsPos-1])
		if actualRef != cdsRef {
			return nil, fmt.Errorf("CDS reference mismatch at position %d: expected %s, got %s in transcript %s",
				cdsPos, cdsRef, actualRef, transcript.ID)
		}
	}

	genomicPos := CDSToGenomic(cdsPos, transcript)
	if genomicPos == 0 {
		return nil, fmt.Errorf("CDS position %d could not be mapped to genomic coordinates in transcript %s",
			cdsPos, transcript.ID)
	}

	ref := cdsRef
	alt := cdsAlt
	if transcript.IsReverseStrand() {
		ref = string(Complement(cdsRef[0]))
		alt = string(Complement(cdsAlt[0]))
	}

	return []*vcf.Variant{{
		Chrom: transcript.Chrom,
		Pos:   genomicPos,
		Ref:   ref,
		Alt:   alt,
	}}, nil
}

// reverseMapHGVScDeletion handles CDS deletions like "923del" or "100_102del".
// It returns all equivalent genomic variants — in repeat regions, multiple
// genomic positions produce the same CDS deletion after normalization.
func reverseMapHGVScDeletion(c *cache.Cache, geneOrTranscript string, m []string) ([]*vcf.Variant, error) {
	cdsStart, _ := strconv.ParseInt(m[1], 10, 64)
	cdsEnd := cdsStart
	if m[2] != "" {
		cdsEnd, _ = strconv.ParseInt(m[2], 10, 64)
	}
	if cdsEnd < cdsStart {
		return nil, fmt.Errorf("invalid CDS deletion range: %d_%d", cdsStart, cdsEnd)
	}

	transcript, err := findHGVScTranscript(c, geneOrTranscript)
	if err != nil {
		return nil, err
	}

	if cdsEnd > int64(len(transcript.CDSSequence)) {
		return nil, fmt.Errorf("CDS position %d out of range for transcript %s (CDS length %d)",
			cdsEnd, transcript.ID, len(transcript.CDSSequence))
	}

	// Find all equivalent CDS positions (repeat-aware)
	positions := equivalentCDSDeletionPositions(transcript.CDSSequence, cdsStart, cdsEnd)

	var variants []*vcf.Variant
	for _, pos := range positions {
		end := pos + (cdsEnd - cdsStart)
		v, err := buildCDSDeletionVariant(transcript, pos, end)
		if err != nil {
			continue // skip positions where padding is out of range
		}
		variants = append(variants, v)
	}

	if len(variants) == 0 {
		return nil, fmt.Errorf("CDS positions %d-%d could not be mapped to genomic coordinates in transcript %s",
			cdsStart, cdsEnd, transcript.ID)
	}
	return variants, nil
}

// equivalentCDSDeletionPositions finds all CDS start positions where deleting
// the same number of bases produces an identical CDS sequence. In a repeat
// region, a deletion can be placed at any position within the repeat.
// Returns 1-based CDS start positions in ascending order.
func equivalentCDSDeletionPositions(cdsSeq string, cdsStart, cdsEnd int64) []int64 {
	L := int(cdsEnd - cdsStart + 1)
	i := int(cdsStart - 1) // 0-indexed

	// Shift left: deleting [i-1, i-1+L) ≡ [i, i+L) when cdsSeq[i-1] == cdsSeq[i+L-1]
	left := i
	for left > 0 && cdsSeq[left-1] == cdsSeq[left+L-1] {
		left--
	}

	// Shift right: deleting [i+1, i+1+L) ≡ [i, i+L) when cdsSeq[i] == cdsSeq[i+L]
	right := i
	for right+L < len(cdsSeq) && cdsSeq[right] == cdsSeq[right+L] {
		right++
	}

	positions := make([]int64, 0, right-left+1)
	for pos := left; pos <= right; pos++ {
		positions = append(positions, int64(pos+1)) // 1-based
	}
	return positions
}

// buildCDSDeletionVariant builds a single VCF-convention deletion variant
// from a CDS deletion at the given positions.
func buildCDSDeletionVariant(t *cache.Transcript, cdsStart, cdsEnd int64) (*vcf.Variant, error) {
	deletedBases := t.CDSSequence[cdsStart-1 : cdsEnd]

	genomicStart := CDSToGenomic(cdsStart, t)
	genomicEnd := CDSToGenomic(cdsEnd, t)
	if genomicStart == 0 || genomicEnd == 0 {
		return nil, fmt.Errorf("CDS positions %d-%d unmapped", cdsStart, cdsEnd)
	}

	if t.IsReverseStrand() {
		padPos := genomicEnd - 1
		padCDSPos := cdsEnd + 1
		if padPos < 1 || padCDSPos > int64(len(t.CDSSequence)) {
			return nil, fmt.Errorf("padding base out of range")
		}
		padBase := string(Complement(t.CDSSequence[padCDSPos-1]))
		genomicDeleted := ReverseComplement(string(deletedBases))
		return &vcf.Variant{
			Chrom: t.Chrom,
			Pos:   padPos,
			Ref:   padBase + genomicDeleted,
			Alt:   padBase,
		}, nil
	}

	padPos := genomicStart - 1
	padCDSPos := cdsStart - 1
	if padPos < 1 || padCDSPos < 1 {
		return nil, fmt.Errorf("padding base out of range")
	}
	padBase := string(t.CDSSequence[padCDSPos-1])
	return &vcf.Variant{
		Chrom: t.Chrom,
		Pos:   padPos,
		Ref:   padBase + string(deletedBases),
		Alt:   padBase,
	}, nil
}

// reGenomicSubstitution parses a genomic substitution like "1293968C>T".
var reGenomicSubstitution = regexp.MustCompile(`^(\d+)([ACGT])>([ACGT])$`)

// reGenomicDeletion parses a genomic deletion like "1293968del" or "1293968_1293970del".
var reGenomicDeletion = regexp.MustCompile(`^(\d+)(?:_(\d+))?del$`)

// reGenomicInsertion parses a genomic insertion like "41242962_41242963insGA".
var reGenomicInsertion = regexp.MustCompile(`^(\d+)_(\d+)ins([ACGTacgt]+)$`)

// reGenomicDelIns parses a genomic deletion-insertion like "7_8delinsAA" or "7del1insAA".
var reGenomicDelIns = regexp.MustCompile(`^(\d+)(?:_(\d+))?delins([ACGTacgt]+)$`)

// reGenomicDup parses a genomic duplication like "41242962dup" or "41242962_41242964dup".
var reGenomicDup = regexp.MustCompile(`^(\d+)(?:_(\d+))?dup$`)

// ResolveHGVSg resolves an HGVSg notation (e.g. "1293968del") on a given chromosome
// to a VCF-convention variant. For deletions, insertions, and other notations, it
// looks up reference bases from a transcript's CDS sequence.
func ResolveHGVSg(c *cache.Cache, chrom string, genomicChange string) ([]*vcf.Variant, error) {
	// Try substitution first (bases are given, no lookup needed)
	if m := reGenomicSubstitution.FindStringSubmatch(genomicChange); m != nil {
		pos, _ := strconv.ParseInt(m[1], 10, 64)
		return []*vcf.Variant{{
			Chrom: chrom,
			Pos:   pos,
			Ref:   m[2],
			Alt:   m[3],
		}}, nil
	}

	// Try insertion (e.g. "41242962_41242963insGA")
	if m := reGenomicInsertion.FindStringSubmatch(genomicChange); m != nil {
		return resolveHGVSgInsertion(c, chrom, m)
	}

	// Try deletion-insertion (e.g. "7_8delinsAA")
	if m := reGenomicDelIns.FindStringSubmatch(genomicChange); m != nil {
		return resolveHGVSgDelIns(c, chrom, m)
	}

	// Try duplication (e.g. "41242962dup" or "41242962_41242964dup")
	if m := reGenomicDup.FindStringSubmatch(genomicChange); m != nil {
		return resolveHGVSgDup(c, chrom, m)
	}

	// Try deletion
	if m := reGenomicDeletion.FindStringSubmatch(genomicChange); m != nil {
		return resolveHGVSgDeletion(c, chrom, m)
	}

	return nil, fmt.Errorf("unsupported genomic change notation %q", genomicChange)
}

// resolveHGVSgInsertion resolves "pos1_pos2insSeq" to a VCF-convention variant.
// VCF convention: ref = base at pos1, alt = base at pos1 + inserted sequence.
func resolveHGVSgInsertion(c *cache.Cache, chrom string, m []string) ([]*vcf.Variant, error) {
	pos1, _ := strconv.ParseInt(m[1], 10, 64)
	insertedSeq := strings.ToUpper(m[3])

	// Look up the reference base at pos1 from a transcript
	refBase, err := lookupRefBase(c, chrom, pos1)
	if err != nil {
		return nil, fmt.Errorf("insertion at %s:%s: %w", chrom, m[1], err)
	}

	return []*vcf.Variant{{
		Chrom: chrom,
		Pos:   pos1,
		Ref:   refBase,
		Alt:   refBase + insertedSeq,
	}}, nil
}

// resolveHGVSgDelIns resolves "pos1_pos2delinsSeq" or "pos1delinsSeq".
func resolveHGVSgDelIns(c *cache.Cache, chrom string, m []string) ([]*vcf.Variant, error) {
	start, _ := strconv.ParseInt(m[1], 10, 64)
	end := start
	if m[2] != "" {
		end, _ = strconv.ParseInt(m[2], 10, 64)
	}
	insertedSeq := strings.ToUpper(m[3])

	// Look up padding base and deleted reference bases
	padPos := start - 1
	padBase, err := lookupRefBase(c, chrom, padPos)
	if err != nil {
		return nil, fmt.Errorf("delins at %s:%d: %w", chrom, start, err)
	}

	refBases := padBase
	for pos := start; pos <= end; pos++ {
		b, err := lookupRefBase(c, chrom, pos)
		if err != nil {
			return nil, fmt.Errorf("delins at %s:%d: %w", chrom, pos, err)
		}
		refBases += b
	}

	return []*vcf.Variant{{
		Chrom: chrom,
		Pos:   padPos,
		Ref:   refBases,
		Alt:   padBase + insertedSeq,
	}}, nil
}

// resolveHGVSgDup resolves "posdup" or "pos1_pos2dup".
func resolveHGVSgDup(c *cache.Cache, chrom string, m []string) ([]*vcf.Variant, error) {
	start, _ := strconv.ParseInt(m[1], 10, 64)
	end := start
	if m[2] != "" {
		end, _ = strconv.ParseInt(m[2], 10, 64)
	}

	// Look up the duplicated bases
	dupBases := ""
	for pos := start; pos <= end; pos++ {
		b, err := lookupRefBase(c, chrom, pos)
		if err != nil {
			return nil, fmt.Errorf("dup at %s:%d: %w", chrom, pos, err)
		}
		dupBases += b
	}

	// VCF: position before the dup, ref = pad base, alt = pad base + dup bases
	padPos := start - 1
	padBase, err := lookupRefBase(c, chrom, padPos)
	if err != nil {
		return nil, fmt.Errorf("dup at %s:%d: %w", chrom, padPos, err)
	}

	return []*vcf.Variant{{
		Chrom: chrom,
		Pos:   padPos,
		Ref:   padBase,
		Alt:   padBase + dupBases,
	}}, nil
}

// lookupRefBase looks up a single reference base at a genomic position
// by finding a transcript covering that position.
func lookupRefBase(c *cache.Cache, chrom string, pos int64) (string, error) {
	transcripts := c.FindTranscripts(chrom, pos)
	for _, t := range transcripts {
		if len(t.CDSSequence) == 0 {
			continue
		}
		cdsPos := GenomicToCDS(pos, t)
		if cdsPos == 0 {
			continue
		}
		idx := int(cdsPos) - 1
		if idx >= 0 && idx < len(t.CDSSequence) {
			if t.IsReverseStrand() {
				return string(Complement(t.CDSSequence[idx])), nil
			}
			return string(t.CDSSequence[idx]), nil
		}
	}
	return "", fmt.Errorf("could not resolve reference base at %s:%d (no CDS coverage)", chrom, pos)
}

// resolveHGVSgDeletion resolves a genomic deletion by looking up reference bases
// from a transcript's CDS sequence.
func resolveHGVSgDeletion(c *cache.Cache, chrom string, m []string) ([]*vcf.Variant, error) {
	genomicStart, _ := strconv.ParseInt(m[1], 10, 64)
	genomicEnd := genomicStart
	if m[2] != "" {
		genomicEnd, _ = strconv.ParseInt(m[2], 10, 64)
	}
	if genomicEnd < genomicStart {
		return nil, fmt.Errorf("invalid genomic deletion range: %d_%d", genomicStart, genomicEnd)
	}

	// Find a transcript covering the deletion and padding position
	padPos := genomicStart - 1
	transcripts := c.FindTranscripts(chrom, genomicStart)
	if len(transcripts) == 0 {
		return nil, fmt.Errorf("no transcript found covering %s:%d", chrom, genomicStart)
	}

	// Try each transcript until we find one that can resolve all bases
	for _, t := range transcripts {
		if len(t.CDSSequence) == 0 {
			continue
		}

		// Check that padding position and all deleted positions map to CDS
		padCDS := GenomicToCDS(padPos, t)
		if padCDS == 0 {
			continue
		}

		allMapped := true
		for pos := genomicStart; pos <= genomicEnd; pos++ {
			if GenomicToCDS(pos, t) == 0 {
				allMapped = false
				break
			}
		}
		if !allMapped {
			continue
		}

		// Get padding base (always on genomic/forward strand)
		var padBase string
		if t.IsReverseStrand() {
			padBase = string(Complement(t.CDSSequence[padCDS-1]))
		} else {
			padBase = string(t.CDSSequence[padCDS-1])
		}

		// Get deleted bases on genomic strand
		var deletedBases []byte
		for pos := genomicStart; pos <= genomicEnd; pos++ {
			cdsPos := GenomicToCDS(pos, t)
			base := t.CDSSequence[cdsPos-1]
			if t.IsReverseStrand() {
				base = Complement(base)
			}
			deletedBases = append(deletedBases, base)
		}

		return []*vcf.Variant{{
			Chrom: chrom,
			Pos:   padPos,
			Ref:   padBase + string(deletedBases),
			Alt:   padBase,
		}}, nil
	}

	return nil, fmt.Errorf("no transcript with CDS sequence covers positions %s:%d-%d (including padding base)", chrom, padPos, genomicEnd)
}

// findTranscriptByPrefix finds transcripts matching an ID prefix (without version).
func findTranscriptByPrefix(c *cache.Cache, prefix string) []*cache.Transcript {
	// Strip version suffix from prefix if present for matching
	base := prefix
	if idx := strings.IndexByte(prefix, '.'); idx >= 0 {
		base = prefix[:idx]
	}

	var matches []*cache.Transcript
	for _, chrom := range c.Chromosomes() {
		for _, t := range c.FindTranscriptsByChrom(chrom) {
			tid := t.ID
			if idx := strings.IndexByte(tid, '.'); idx >= 0 {
				tid = tid[:idx]
			}
			if tid == base {
				matches = append(matches, t)
			}
		}
	}
	return matches
}
