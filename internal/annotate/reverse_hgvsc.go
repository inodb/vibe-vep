package annotate

import (
	"fmt"
	"regexp"
	"strconv"

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
