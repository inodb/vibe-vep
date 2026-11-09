package annotate

import (
	"fmt"
	"regexp"
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/cache"
)

// Resolving an HGVS c. notation back to a genomic position is the inverse of
// what hgvsc.go does on the way out, and this file is written as that inverse:
// genomicToHGVScPos picks between an exonic, UTR and intronic form, and
// cdsCoordToGenomic undoes each of them.
//
// Only the coordinate is handled here. What the variant *does* at that
// coordinate -- substitute, delete, insert -- is the same question whether it
// was written in c. or g. notation, and is answered by the existing genomic
// edit builders once both ends are known.

// cdsCoord is one position in HGVS c. numbering.
//
// A coordinate names an exonic base and may then step off it into an intron:
// c.1899-1 is one base before the exonic base numbered 1899, and c.-29+2 is
// two bases after the base numbered -29, which is itself 29 bases before the
// start of the coding sequence.
type cdsCoord struct {
	// Pos is the base number as written, negative in the 5'UTR.
	Pos int64
	// UTR3 marks the "*" prefix, counting from the end of the coding
	// sequence rather than its start.
	UTR3 bool
	// Offset steps into the neighbouring intron, in transcript direction:
	// positive past the anchor base, negative before it.
	Offset int64
}

// cdsCoordRe matches a single HGVS c. coordinate: an optional "*" or "-"
// marker, the base number, and an optional signed intron offset.
var cdsCoordRe = regexp.MustCompile(`^(\*)?(-)?(\d+)(?:([+-])(\d+))?$`)

// parseCDSCoord reads one coordinate, e.g. "1799", "-29", "*8", "1899-1",
// "-29+2".
func parseCDSCoord(s string) (cdsCoord, bool) {
	m := cdsCoordRe.FindStringSubmatch(s)
	if m == nil {
		return cdsCoord{}, false
	}

	n, err := strconv.ParseInt(m[3], 10, 64)
	if err != nil {
		return cdsCoord{}, false
	}

	c := cdsCoord{Pos: n, UTR3: m[1] == "*"}
	if m[2] == "-" {
		c.Pos = -n
	}
	if m[4] != "" {
		off, err := strconv.ParseInt(m[5], 10, 64)
		if err != nil {
			return cdsCoord{}, false
		}
		if m[4] == "-" {
			off = -off
		}
		c.Offset = off
	}
	return c, true
}

// transcriptPosToGenomic is the inverse of GenomicToTranscriptPos: it turns a
// 1-based position along the spliced transcript back into a genomic
// coordinate. Returns 0 when the position falls outside the transcript.
func transcriptPosToGenomic(txPos int64, t *cache.Transcript) int64 {
	if txPos < 1 || len(t.Exons) == 0 || len(t.ExonCumBases) != len(t.Exons) {
		return 0
	}

	for i := range t.Exons {
		exon := &t.Exons[i]
		length := exon.End - exon.Start + 1
		before := t.ExonCumBases[i]
		if txPos <= before || txPos > before+length {
			continue
		}
		offset := txPos - before - 1
		if t.Strand == 1 {
			return exon.Start + offset
		}
		return exon.End - offset
	}
	return 0
}

// cdsStartTranscriptPos returns where the coding sequence begins along the
// spliced transcript, and cdsEndTranscriptPos where it ends. On the reverse
// strand the coding sequence starts at the genomically higher coordinate.
func cdsStartTranscriptPos(t *cache.Transcript) int64 {
	if t.Strand == 1 {
		return GenomicToTranscriptPos(t.CDSStart, t)
	}
	return GenomicToTranscriptPos(t.CDSEnd, t)
}

func cdsEndTranscriptPos(t *cache.Transcript) int64 {
	if t.Strand == 1 {
		return GenomicToTranscriptPos(t.CDSEnd, t)
	}
	return GenomicToTranscriptPos(t.CDSStart, t)
}

// cdsCoordToGenomic resolves a c. coordinate to a genomic position.
//
// The anchor comes first -- the exonic base the coordinate names -- and the
// intron offset is then applied in transcript direction, which on the reverse
// strand runs towards lower genomic coordinates.
func cdsCoordToGenomic(c cdsCoord, t *cache.Transcript) (int64, error) {
	anchor, err := cdsAnchorToGenomic(c, t)
	if err != nil {
		return 0, err
	}

	if c.Offset == 0 {
		return anchor, nil
	}
	if t.Strand == 1 {
		return anchor + c.Offset, nil
	}
	return anchor - c.Offset, nil
}

// cdsAnchorToGenomic resolves the exonic base a coordinate names, ignoring any
// intron offset.
func cdsAnchorToGenomic(c cdsCoord, t *cache.Transcript) (int64, error) {
	switch {
	case c.UTR3:
		// "*N" counts N bases past the last base of the coding sequence.
		end := cdsEndTranscriptPos(t)
		if end == 0 {
			return 0, fmt.Errorf("transcript %s has no coding sequence to count from", t.ID)
		}
		pos := transcriptPosToGenomic(end+c.Pos, t)
		if pos == 0 {
			return 0, fmt.Errorf("c.*%d is past the end of transcript %s", c.Pos, t.ID)
		}
		return pos, nil

	case c.Pos < 0:
		// "-N" counts N bases before the first base of the coding sequence.
		start := cdsStartTranscriptPos(t)
		if start == 0 {
			return 0, fmt.Errorf("transcript %s has no coding sequence to count from", t.ID)
		}
		pos := transcriptPosToGenomic(start+c.Pos, t)
		if pos == 0 {
			return 0, fmt.Errorf("c.%d is before the start of transcript %s", c.Pos, t.ID)
		}
		return pos, nil

	default:
		pos := CDSToGenomic(c.Pos, t)
		if pos == 0 {
			return 0, fmt.Errorf("c.%d is outside the coding sequence of transcript %s", c.Pos, t.ID)
		}
		return pos, nil
	}
}

// cdsChangeRe splits an HGVS c. change into its coordinates and the edit
// applied there: one or two positions, then "A>G", "del", "ins…", "delins…"
// or "dup".
var cdsChangeRe = regexp.MustCompile(`^([*\-]?\d+(?:[+-]\d+)?)(?:_([*\-]?\d+(?:[+-]\d+)?))?(.*)$`)

// cdsChangeToGenomic rewrites a transcript-relative change as the genomic
// change it describes, which the genomic resolver then reads as usual. A
// variant does the same thing whichever way it is written; only the
// coordinates and, on the reverse strand, the strand of the sequence differ.
func cdsChangeToGenomic(cdsChange string, t *cache.Transcript) (string, error) {
	m := cdsChangeRe.FindStringSubmatch(cdsChange)
	if m == nil {
		return "", fmt.Errorf("cannot read HGVS c. change %q", cdsChange)
	}

	startCoord, ok := parseCDSCoord(m[1])
	if !ok {
		return "", fmt.Errorf("cannot read coordinate %q in %q", m[1], cdsChange)
	}
	start, err := cdsCoordToGenomic(startCoord, t)
	if err != nil {
		return "", err
	}

	end := start
	if m[2] != "" {
		endCoord, ok := parseCDSCoord(m[2])
		if !ok {
			return "", fmt.Errorf("cannot read coordinate %q in %q", m[2], cdsChange)
		}
		if end, err = cdsCoordToGenomic(endCoord, t); err != nil {
			return "", err
		}
	}

	// Transcript order runs backwards along the genome on the reverse strand,
	// so the pair has to be put back into genomic order.
	if start > end {
		start, end = end, start
	}

	edit, err := genomicEdit(m[3], t.IsReverseStrand())
	if err != nil {
		return "", fmt.Errorf("%w in %q", err, cdsChange)
	}

	if start == end && !strings.HasPrefix(edit, "ins") {
		return strconv.FormatInt(start, 10) + edit, nil
	}
	return strconv.FormatInt(start, 10) + "_" + strconv.FormatInt(end, 10) + edit, nil
}

// genomicEdit rewrites the edit half of a c. change for the genome, turning
// any sequence it carries round when the transcript is on the reverse strand.
func genomicEdit(edit string, reverse bool) (string, error) {
	switch {
	case reSubstitutionEdit.MatchString(edit):
		m := reSubstitutionEdit.FindStringSubmatch(edit)
		ref, alt := m[1], m[2]
		if reverse {
			ref, alt = reverseComplement(ref), reverseComplement(alt)
		}
		return ref + ">" + alt, nil

	case strings.HasPrefix(edit, "delins"):
		return "delins" + orientSequence(edit[len("delins"):], reverse), nil

	case strings.HasPrefix(edit, "del"):
		return "del" + orientSequence(edit[len("del"):], reverse), nil

	case strings.HasPrefix(edit, "ins"):
		return "ins" + orientSequence(edit[len("ins"):], reverse), nil

	case strings.HasPrefix(edit, "dup"):
		return "dup" + orientSequence(edit[len("dup"):], reverse), nil

	default:
		return "", fmt.Errorf("unsupported edit %q", edit)
	}
}

// reSubstitutionEdit matches the "A>G" half of a substitution.
var reSubstitutionEdit = regexp.MustCompile(`^([ACGTN]+)>([ACGTN]+)$`)

// orientSequence turns a sequence round for the opposite strand, leaving an
// absent one absent.
func orientSequence(seq string, reverse bool) string {
	if seq == "" || !reverse {
		return seq
	}
	return reverseComplement(seq)
}

// reverseComplement returns the other strand's reading of a sequence.
func reverseComplement(seq string) string {
	out := make([]byte, 0, len(seq))
	for i := len(seq) - 1; i >= 0; i-- {
		out = append(out, Complement(seq[i]))
	}
	return string(out)
}
