package annotate

import (
	"fmt"
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

	switch conseq {
	case ConsequenceMissenseVariant:
		return fmt.Sprintf("p.%s%d%s", aaThree(result.RefAA), pos, aaThree(result.AltAA))

	case ConsequenceSynonymousVariant:
		return fmt.Sprintf("p.%s%d=", aaThree(result.RefAA), pos)

	case ConsequenceStopGained:
		return fmt.Sprintf("p.%s%dTer", aaThree(result.RefAA), pos)

	case ConsequenceStopLost:
		if result.StopLostExtDist > 0 {
			return fmt.Sprintf("p.Ter%d%sext*%d", pos, aaThree(result.AltAA), result.StopLostExtDist)
		}
		return fmt.Sprintf("p.Ter%d%sext*?", pos, aaThree(result.AltAA))

	case ConsequenceStartLost:
		return "p.Met1?"

	case ConsequenceStopRetained:
		return fmt.Sprintf("p.Ter%d=", pos)

	case ConsequenceFrameshiftVariant:
		if result.RefAA != 0 {
			altStr := ""
			if result.AltAA != 0 {
				altStr = aaThree(result.AltAA)
			}
			if result.FrameshiftStopDist > 0 {
				return fmt.Sprintf("p.%s%d%sfsTer%d", aaThree(result.RefAA), pos, altStr, result.FrameshiftStopDist)
			}
			return fmt.Sprintf("p.%s%d%sfs", aaThree(result.RefAA), pos, altStr)
		}
		return fmt.Sprintf("p.%dfs", pos)

	case ConsequenceInframeDeletion:
		if result.RefAA != 0 {
			return fmt.Sprintf("p.%s%ddel", aaThree(result.RefAA), pos)
		}
		return fmt.Sprintf("p.%ddel", pos)

	case ConsequenceInframeInsertion:
		if result.RefAA != 0 {
			return fmt.Sprintf("p.%s%d_%dins", aaThree(result.RefAA), pos, pos+1)
		}
		return fmt.Sprintf("p.%d_%dins", pos, pos+1)

	default:
		return ""
	}
}
