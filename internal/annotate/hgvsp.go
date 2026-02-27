package annotate

import (
	"strconv"
	"strings"
)

// aaThree converts a single-letter amino acid code to its three-letter code.
// Returns "Xaa" for unknown amino acids.
func aaThree(aa byte) string {
	if three, ok := AminoAcidSingleToThree[aa]; ok {
		return three
	}
	return "Xaa"
}

// formatAASequence converts a string of single-letter amino acid codes to
// concatenated three-letter codes (e.g., "AL" → "AlaLeu").
func formatAASequence(aas string) string {
	if len(aas) == 0 {
		return ""
	}
	var b strings.Builder
	for i := 0; i < len(aas); i++ {
		b.WriteString(aaThree(aas[i]))
	}
	return b.String()
}

// FormatHGVSp formats the HGVS protein notation for a consequence result.
// Returns an empty string for non-coding consequences.
func FormatHGVSp(result *ConsequenceResult) string {
	if result.ProteinPosition < 1 {
		return ""
	}

	// Get the primary consequence (first term before comma)
	conseq := result.Consequence
	if idx := strings.IndexByte(conseq, ','); idx >= 0 {
		conseq = conseq[:idx]
	}

	pos := result.ProteinPosition

	posStr := strconv.FormatInt(pos, 10)

	switch conseq {
	case ConsequenceMissenseVariant:
		return "p." + aaThree(result.RefAA) + posStr + aaThree(result.AltAA)

	case ConsequenceSynonymousVariant:
		return "p." + aaThree(result.RefAA) + posStr + "="

	case ConsequenceStopGained:
		return "p." + aaThree(result.RefAA) + posStr + "Ter"

	case ConsequenceStopLost:
		if result.StopLostExtDist > 0 {
			return "p.Ter" + posStr + aaThree(result.AltAA) + "ext*" + strconv.Itoa(result.StopLostExtDist)
		}
		return "p.Ter" + posStr + aaThree(result.AltAA) + "ext*?"

	case ConsequenceStartLost:
		return "p.Met1?"

	case ConsequenceStopRetained:
		return "p.Ter" + posStr + "="

	case ConsequenceFrameshiftVariant:
		if result.RefAA != 0 {
			altStr := ""
			if result.AltAA != 0 {
				altStr = aaThree(result.AltAA)
			}
			if result.FrameshiftStopDist > 0 {
				return "p." + aaThree(result.RefAA) + posStr + altStr + "fsTer" + strconv.Itoa(result.FrameshiftStopDist)
			}
			return "p." + aaThree(result.RefAA) + posStr + altStr + "fs"
		}
		return "p." + posStr + "fs"

	case ConsequenceInframeDeletion:
		suffix := "del"
		if len(result.InsertedAAs) > 0 {
			suffix = "delins" + formatAASequence(result.InsertedAAs)
		}
		if result.ProteinEndPosition > result.ProteinPosition && result.EndAA != 0 {
			endPosStr := strconv.FormatInt(result.ProteinEndPosition, 10)
			if result.RefAA != 0 {
				return "p." + aaThree(result.RefAA) + posStr + "_" + aaThree(result.EndAA) + endPosStr + suffix
			}
			return "p." + posStr + "_" + endPosStr + suffix
		}
		if result.RefAA != 0 {
			return "p." + aaThree(result.RefAA) + posStr + suffix
		}
		return "p." + posStr + suffix

	case ConsequenceInframeInsertion:
		if result.IsDup {
			if result.ProteinEndPosition > result.ProteinPosition && result.EndAA != 0 {
				endPosStr := strconv.FormatInt(result.ProteinEndPosition, 10)
				return "p." + aaThree(result.RefAA) + posStr + "_" + aaThree(result.EndAA) + endPosStr + "dup"
			}
			return "p." + aaThree(result.RefAA) + posStr + "dup"
		}
		pos1Str := strconv.FormatInt(pos+1, 10)
		insAAs := formatAASequence(result.InsertedAAs)
		if result.RefAA != 0 {
			return "p." + aaThree(result.RefAA) + posStr + "_" + pos1Str + "ins" + insAAs
		}
		return "p." + posStr + "_" + pos1Str + "ins" + insAAs

	case ConsequenceSpliceDonor, ConsequenceSpliceAcceptor:
		return "p.X" + posStr + "_splice"

	default:
		return ""
	}
}
