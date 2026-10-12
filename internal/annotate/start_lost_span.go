package annotate

import (
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/cache"
)

// nameStartLostSpan rewrites the HGVSp of a start-lost variant that covers
// more than one codon into the span form VEP uses for them.
//
// VEP writes p.Met1? while the change stays inside one codon -- 495 of 495
// substitutions, and the indels confined to a single codon with them -- and
// switches to p.<residues><first>_?<last> as soon as it covers two or more.
// All 14 such variants in a VEP111 MSK-IMPACT MAF take the span form and
// nothing else does.
//
// How much of the reference it names depends on the frame. An in-frame change
// gets every residue it covers, so PTEN c.3_26del is
// p.MetThrAlaIleIleLysGluIleVal1_?9; a frameshift gets only the first, so
// HLA-A c.3_54del is p.Met1_?18. A duplication is named from the codon after
// the duplicated block, which is how PIK3R2 c.6_9dup comes out p.Pro4_?3 with
// the second number lower than the first.
//
// The one case this does not reproduce is CDKN1A c.12_34del, which VEP calls
// p.Ala5_?12 where the deletion's first base is in codon 4.
func nameStartLostSpan(t *cache.Transcript, result *ConsequenceResult) {
	if !hasTerm(result.Consequence, ConsequenceStartLost) || result.SuppressHGVSp {
		return
	}
	first, last, op, insLen, ok := cdsSpanOfHGVSc(result.HGVSc)
	if !ok {
		return
	}

	// A duplication or insertion names the codon after the block it follows.
	start, end := codonOf(first), codonOf(last)
	netLen := insLen - (last - first + 1)
	switch op {
	case "dup":
		start, end = codonOf(last+1), codonOf(last)
		netLen = last - first + 1
	case "ins":
		start, end = codonOf(last), codonOf(first)
		netLen = insLen
	}
	if start == end {
		return // one codon: VEP keeps the plain p.Met1? form
	}

	peptide := ""
	if netLen%3 == 0 && start <= end {
		for p := start; p <= end; p++ {
			aa := residueAt(t, p)
			if aa == 0 {
				return
			}
			peptide += aaThree(aa)
		}
	} else {
		aa := residueAt(t, start)
		if aa == 0 {
			return
		}
		peptide = aaThree(aa)
	}
	result.HGVSp = "p." + peptide + strconv.FormatInt(start, 10) +
		"_?" + strconv.FormatInt(end, 10)
}

// residueAt translates the transcript's codon at a protein position.
func residueAt(t *cache.Transcript, pos int64) byte {
	codon := GetCodon(t.CDSSequence, pos)
	if len(codon) != 3 {
		return 0
	}
	return TranslateCodon(codon)
}

// codonOf gives the protein position a CDS coordinate falls in.
func codonOf(cds int64) int64 { return (cds + 2) / 3 }

// cdsSpanOfHGVSc reads the coding coordinates out of an HGVSc, along with the
// operator and how many bases it inserts. It reports false unless both
// coordinates are plain CDS numbers, since anything else has no codon.
func cdsSpanOfHGVSc(hgvsc string) (first, last int64, op string, insLen int64, ok bool) {
	if i := strings.Index(hgvsc, ":c."); i >= 0 {
		hgvsc = hgvsc[i+3:]
	} else if strings.HasPrefix(hgvsc, "c.") {
		hgvsc = hgvsc[2:]
	} else {
		return 0, 0, "", 0, false
	}
	k := strings.IndexFunc(hgvsc, func(r rune) bool { return r >= 'a' && r <= 'z' })
	if k < 0 {
		return 0, 0, "", 0, false
	}
	rangeStr, rest := hgvsc[:k], hgvsc[k:]
	for _, cand := range []string{"delins", "del", "dup", "ins", "inv"} {
		if strings.HasPrefix(rest, cand) {
			op = cand
			insLen = int64(len(rest) - len(cand))
			break
		}
	}
	if op == "" {
		return 0, 0, "", 0, false
	}
	if op == "del" || op == "dup" || op == "inv" {
		insLen = 0
	}
	a, b := rangeStr, rangeStr
	if i := strings.IndexByte(rangeStr, '_'); i >= 0 {
		a, b = rangeStr[:i], rangeStr[i+1:]
	}
	first, err1 := strconv.ParseInt(a, 10, 64)
	last, err2 := strconv.ParseInt(b, 10, 64)
	if err1 != nil || err2 != nil || first < 1 || last < first {
		return 0, 0, "", 0, false
	}
	return first, last, op, insLen, true
}
