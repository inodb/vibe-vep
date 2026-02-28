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
		if t.IsCanonical && t.IsProteinCoding() {
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

// reCDSChange parses a CDS change like "35G>T" or "100del" or "100_102del".
// For now we support simple substitutions: position + ref > alt.
var reCDSChange = regexp.MustCompile(`^(\d+)([ACGT])>([ACGT])$`)

// ReverseMapHGVSc maps an HGVSc notation (e.g. KRAS c.35G>T) back to a
// genomic variant using the transcript's CDS-to-genomic mapping.
func ReverseMapHGVSc(c *cache.Cache, geneOrTranscript string, cdsChange string) ([]*vcf.Variant, error) {
	// Parse the CDS change
	m := reCDSChange.FindStringSubmatch(cdsChange)
	if m == nil {
		return nil, fmt.Errorf("unsupported CDS change notation %q (only simple substitutions like 35G>T are supported)", cdsChange)
	}

	cdsPos, _ := strconv.ParseInt(m[1], 10, 64)
	cdsRef := m[2]
	cdsAlt := m[3]

	// Find the transcript
	var transcript *cache.Transcript

	// Check if it looks like a transcript ID (starts with ENST)
	if strings.HasPrefix(geneOrTranscript, "ENST") {
		// Try exact match first, then try with version stripped
		transcript = c.GetTranscript(geneOrTranscript)
		if transcript == nil {
			// Try matching without version suffix
			transcripts := findTranscriptByPrefix(c, geneOrTranscript)
			if len(transcripts) > 0 {
				transcript = transcripts[0]
			}
		}
		if transcript == nil {
			return nil, fmt.Errorf("transcript %q not found", geneOrTranscript)
		}
	} else {
		// Treat as gene name, find canonical
		transcripts := c.FindTranscriptsByGene(geneOrTranscript)
		if len(transcripts) == 0 {
			return nil, fmt.Errorf("gene %q not found in transcript cache", geneOrTranscript)
		}
		for _, t := range transcripts {
			if t.IsCanonical && t.IsProteinCoding() {
				transcript = t
				break
			}
		}
		if transcript == nil {
			for _, t := range transcripts {
				if t.IsProteinCoding() {
					transcript = t
					break
				}
			}
		}
		if transcript == nil {
			return nil, fmt.Errorf("no protein-coding transcript found for gene %q", geneOrTranscript)
		}
	}

	// Verify CDS ref base
	if len(transcript.CDSSequence) > 0 && cdsPos <= int64(len(transcript.CDSSequence)) {
		actualRef := string(transcript.CDSSequence[cdsPos-1])
		if actualRef != cdsRef {
			return nil, fmt.Errorf("CDS reference mismatch at position %d: expected %s, got %s in transcript %s",
				cdsPos, cdsRef, actualRef, transcript.ID)
		}
	}

	// Map CDS position to genomic
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
