package annotate

import "strings"

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
	b.Grow(len(aas) * 3)
	for i := 0; i < len(aas); i++ {
		b.WriteString(aaThree(aas[i]))
	}
	return b.String()
}

// FormatHGVSp formats the HGVS protein notation for a consequence result.
// Returns an empty string for non-coding consequences.
//
// Uses a stack-allocated buffer for common short cases to minimize allocations.
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

	// For most cases we can build in a stack-allocated buffer (≤64 bytes).
	// Only complex cases (delins, insertion with long AA sequences) overflow.
	var buf [64]byte
	n := 0

	switch conseq {
	case ConsequenceMissenseVariant:
		// p.Xxx###Xxx (max ~15 chars)
		n += copy(buf[n:], "p.")
		n += copy(buf[n:], aaThree(result.RefAA))
		n += putInt64(buf[n:], pos)
		n += copy(buf[n:], aaThree(result.AltAA))
		return string(buf[:n])

	case ConsequenceSynonymousVariant:
		// p.Xxx###=
		n += copy(buf[n:], "p.")
		n += copy(buf[n:], aaThree(result.RefAA))
		n += putInt64(buf[n:], pos)
		buf[n] = '='
		n++
		return string(buf[:n])

	case ConsequenceStopGained:
		// p.Xxx###Ter
		n += copy(buf[n:], "p.")
		n += copy(buf[n:], aaThree(result.RefAA))
		n += putInt64(buf[n:], pos)
		n += copy(buf[n:], "Ter")
		return string(buf[:n])

	case ConsequenceStopLost:
		// p.Ter###Xxxext*N or p.Ter###Xxxext*?
		n += copy(buf[n:], "p.Ter")
		n += putInt64(buf[n:], pos)
		n += copy(buf[n:], aaThree(result.AltAA))
		if result.StopLostExtDist > 0 {
			n += copy(buf[n:], "ext*")
			n += putInt64(buf[n:], int64(result.StopLostExtDist))
		} else {
			n += copy(buf[n:], "ext*?")
		}
		return string(buf[:n])

	case ConsequenceStartLost:
		return "p.Met1?"

	case ConsequenceStopRetained:
		// p.Ter###=
		n += copy(buf[n:], "p.Ter")
		n += putInt64(buf[n:], pos)
		buf[n] = '='
		n++
		return string(buf[:n])

	case ConsequenceFrameshiftVariant:
		n += copy(buf[n:], "p.")
		if result.RefAA != 0 {
			n += copy(buf[n:], aaThree(result.RefAA))
		}
		n += putInt64(buf[n:], pos)
		if result.AltAA != 0 {
			n += copy(buf[n:], aaThree(result.AltAA))
		}
		if result.FrameshiftStopDist > 0 {
			n += copy(buf[n:], "fsTer")
			n += putInt64(buf[n:], int64(result.FrameshiftStopDist))
		} else {
			n += copy(buf[n:], "fs")
		}
		return string(buf[:n])

	case ConsequenceInframeDeletion:
		n += copy(buf[n:], "p.")
		if result.RefAA != 0 {
			n += copy(buf[n:], aaThree(result.RefAA))
		}
		n += putInt64(buf[n:], pos)
		if result.ProteinEndPosition > result.ProteinPosition && result.EndAA != 0 {
			buf[n] = '_'
			n++
			n += copy(buf[n:], aaThree(result.EndAA))
			n += putInt64(buf[n:], result.ProteinEndPosition)
		}
		if len(result.InsertedAAs) > 0 {
			n += copy(buf[n:], "delins")
			return string(buf[:n]) + formatAASequence(result.InsertedAAs)
		}
		n += copy(buf[n:], "del")
		return string(buf[:n])

	case ConsequenceInframeInsertion:
		n += copy(buf[n:], "p.")
		if result.IsDup {
			if result.RefAA != 0 {
				n += copy(buf[n:], aaThree(result.RefAA))
			}
			n += putInt64(buf[n:], pos)
			if result.ProteinEndPosition > result.ProteinPosition && result.EndAA != 0 {
				buf[n] = '_'
				n++
				n += copy(buf[n:], aaThree(result.EndAA))
				n += putInt64(buf[n:], result.ProteinEndPosition)
			}
			n += copy(buf[n:], "dup")
			return string(buf[:n])
		}
		if result.RefAA != 0 {
			n += copy(buf[n:], aaThree(result.RefAA))
		}
		n += putInt64(buf[n:], pos)
		buf[n] = '_'
		n++
		n += putInt64(buf[n:], pos+1)
		n += copy(buf[n:], "ins")
		return string(buf[:n]) + formatAASequence(result.InsertedAAs)

	case ConsequenceSpliceDonor, ConsequenceSpliceAcceptor:
		n += copy(buf[n:], "p.X")
		n += putInt64(buf[n:], pos)
		n += copy(buf[n:], "_splice")
		return string(buf[:n])

	default:
		return ""
	}
}
