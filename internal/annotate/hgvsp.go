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
	// Get the primary consequence (first term before comma)
	conseq := result.Consequence
	if idx := strings.IndexByte(conseq, ','); idx >= 0 {
		conseq = conseq[:idx]
	}

	// start_lost is always p.Met1? regardless of position; handle before the
	// ProteinPosition guard to be robust against early-return code paths that
	// may not set ProteinPosition.
	if conseq == ConsequenceStartLost {
		// VEP writes p.Met1? for the great majority (516 of 613 start_lost rows
		// in a VEP111 MSK-IMPACT MAF) but names the affected residue when the
		// change lands past the start codon: c.7del is p.V3?, c.4del is p.A2?.
		// Both are the same p.<residue><position>? form; position 1 just
		// happens to be Met. Only diverge when we have both, since some paths
		// reach here without setting ProteinPosition.
		if result.ProteinPosition > 1 && result.RefAA != 0 {
			var b [32]byte
			n := copy(b[:], "p.")
			n += copy(b[n:], aaThree(result.RefAA))
			n += putInt64(b[n:], result.ProteinPosition)
			b[n] = '?'
			n++
			return string(b[:n])
		}
		return "p.Met1?"
	}

	if result.ProteinPosition < 1 {
		return ""
	}

	pos := result.ProteinPosition

	// For most cases we can build in a stack-allocated buffer (≤64 bytes).
	// Only complex cases (delins, insertion with long AA sequences) overflow.
	var buf [64]byte
	n := 0

	switch conseq {
	case ConsequenceMissenseVariant:
		if result.IsDelIns {
			// Multi-codon MNV: p.AA1pos[_AA2endPos]delins[newAAs]
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
			n += copy(buf[n:], "delins")
			return string(buf[:n]) + formatAASequence(result.InsertedAAs)
		}
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
		// A multi-codon MNV that introduces a stop alters more than the single
		// residue, and VEP describes the whole span rather than collapsing it
		// to a nonsense call: c.1125_1126delinsTT is p.K375_Q376delinsN* and
		// c.1284_1285inv is p.Y428_Y429delins*, not p.K375* / p.Y428*.
		// predictMNVConsequence already recorded the span on the delins fields.
		// A pure insertion that introduces a stop still brackets the insertion
		// point: VEP writes p.E288_T289ins* for c.863_864insGTA, not p.E288*.
		// The delins form below covers the case where the anchor codon also
		// changed; this one is the plain insertion.
		if !result.IsDelIns && len(result.InsertedAAs) > 0 && result.RefAA != 0 && result.EndAA != 0 {
			n += copy(buf[n:], "p.")
			n += copy(buf[n:], aaThree(result.RefAA))
			n += putInt64(buf[n:], pos)
			buf[n] = '_'
			n++
			n += copy(buf[n:], aaThree(result.EndAA))
			n += putInt64(buf[n:], pos+1)
			n += copy(buf[n:], "ins")
			ins := result.InsertedAAs
			for i := 0; i < len(ins); i++ {
				if ins[i] == '*' {
					ins = ins[:i+1] // nothing is translated past a stop
					break
				}
			}
			return string(buf[:n]) + formatAASequence(ins)
		}
		if result.IsDelIns && len(result.InsertedAAs) > 0 {
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
			n += copy(buf[n:], "delins")
			// Nothing is translated past a stop, so c.1284_1285inv is
			// p.Y428_Y429delins*, not p.Y428_Y429delins*N.
			ins := result.InsertedAAs
			for i := 0; i < len(ins); i++ {
				if ins[i] == '*' {
					ins = ins[:i+1]
					break
				}
			}
			return string(buf[:n]) + formatAASequence(ins)
		}
		// p.Xxx###Ter
		n += copy(buf[n:], "p.")
		n += copy(buf[n:], aaThree(result.RefAA))
		n += putInt64(buf[n:], pos)
		n += copy(buf[n:], "Ter")
		return string(buf[:n])

	case ConsequenceStopLost:
		// p.Ter###XxxextTer### or p.Ter###XxxextTer?
		//
		// The extension's terminator is spelled out like every other residue
		// in the three-letter form. genome-nexus shortens this to the p.*###
		// Xext*### that the MAF's HGVSp_Short column carries, which is why the
		// short column agreed while this one did not: every one of the 480
		// extensions in a VEP111 MSK-IMPACT MAF writes extTer here.
		n += copy(buf[n:], "p.Ter")
		n += putInt64(buf[n:], pos)
		n += copy(buf[n:], aaThree(result.AltAA))
		n += copy(buf[n:], "extTer")
		if result.StopLostExtDist > 0 {
			n += putInt64(buf[n:], int64(result.StopLostExtDist))
		} else {
			buf[n] = '?'
			n++
		}
		return string(buf[:n])

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
		// A frameshift whose first new codon is already a stop is written as a
		// plain nonsense change: VEP reports p.Ser70Ter, not p.Ser70TerfsTer1.
		if result.AltAA == '*' && result.FrameshiftStopDist == 1 {
			return string(buf[:n])
		}
		if result.FrameshiftStopDist > 0 {
			n += copy(buf[n:], "fsTer")
			n += putInt64(buf[n:], int64(result.FrameshiftStopDist))
		} else if result.AltAA != 0 {
			// No downstream stop found, but the new residue is known: VEP writes
			// fsTer? (p.Pro315GlnfsTer?), shortening to p.P315Qfs*?.
			n += copy(buf[n:], "fsTer?")
		} else {
			// Neither the new residue nor a stop is known — VEP writes a bare
			// fs, as on splice-site deletions (p.K267fs).
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
		if result.IsDelIns {
			// Anchor codon(s) changed: p.AA1pos[_AA2endPos]delins[newAAs]
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
			n += copy(buf[n:], "delins")
			return string(buf[:n]) + formatAASequence(result.InsertedAAs)
		}
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
		if result.EndAA != 0 {
			n += copy(buf[n:], aaThree(result.EndAA))
		}
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
