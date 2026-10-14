package cache

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strings"
)

// IntronFlankBases is how much intron sequence is kept on each side of every
// coding exon.
//
// HGVS places a deletion as far 3' as possible, and testing each step reads the
// base just past the deleted span. When a deletion sits at an exon edge that
// base is intronic, which vibe-vep otherwise has no way to see: transcript
// sequence is spliced, and the MAF's reference allele stops at the deleted
// bases. NPM1 c.510_524del shifts 14 bases to the exon boundary and stalls
// there, where VEP reports c.511_524+1del with offset 15.
//
// The flank must span the deletion's intronic tail AND the shift, not just
// the shift: a deletion already 40 bases into an intron needs base 41 to test
// its first step. Measured against VEP111 on the MSK-IMPACT set the required
// depth is median 69, so 100 reaches about half of them; going deeper is a
// straight size trade (see the sizing table in the commit history).
const IntronFlankBases = 100

// GenomeFlankLoader fills in the intron flanks of coding exons from a genomic
// FASTA. The genome itself is not retained: one chromosome is held at a time,
// the windows are copied out, and the buffer is reused.
type GenomeFlankLoader struct {
	path string
}

// NewGenomeFlankLoader creates a loader reading the given (optionally gzipped)
// genomic FASTA.
func NewGenomeFlankLoader(path string) *GenomeFlankLoader {
	return &GenomeFlankLoader{path: path}
}

// Attach fills IntronBefore/IntronAfter on every coding exon of every
// transcript. Flanks are stored in CODING orientation so callers need not
// re-derive strand.
func (l *GenomeFlankLoader) Attach(byChrom map[string][]*Transcript) error {
	f, err := os.Open(l.path)
	if err != nil {
		return fmt.Errorf("open genome FASTA: %w", err)
	}
	defer f.Close()

	var r io.Reader = f
	if strings.HasSuffix(l.path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return fmt.Errorf("gunzip genome FASTA: %w", err)
		}
		defer gz.Close()
		r = gz
	}

	// GENCODE ships the genome unwrapped: chromosome 1 is a single line of
	// 249,250,621 characters, which no line scanner will read. Stream fixed
	// chunks instead and split on newlines by hand, bulk-copying the runs
	// between them.
	const chunk = 1 << 20
	buf := make([]byte, chunk)

	var cur string
	var header []byte
	inHeader := false
	lineStart := true
	seq := make([]byte, 0, 1<<28) // reused across chromosomes

	flush := func() {
		if cur != "" {
			fillFlanks(seq, byChrom[cur])
		}
		seq = seq[:0]
	}

	for {
		n, readErr := r.Read(buf)
		data := buf[:n]
		for len(data) > 0 {
			if inHeader {
				if j := bytes.IndexByte(data, '\n'); j >= 0 {
					header = append(header, data[:j]...)
					data = data[j+1:]
					inHeader = false
					lineStart = true
					flush()
					name := string(header)
					if i := bytes.IndexAny([]byte(name), " \t"); i >= 0 {
						name = name[:i]
					}
					cur = normalizeChrom(name)
					if _, wanted := byChrom[cur]; !wanted {
						cur = "" // skip scaffolds and patches
					}
					header = header[:0]
				} else {
					header = append(header, data...)
					data = nil
				}
				continue
			}
			if lineStart && data[0] == '>' {
				inHeader = true
				lineStart = false
				data = data[1:]
				continue
			}
			j := bytes.IndexByte(data, '\n')
			if j < 0 {
				if cur != "" {
					seq = append(seq, data...)
				}
				lineStart = false
				data = nil
				continue
			}
			if cur != "" {
				seq = append(seq, data[:j]...)
			}
			lineStart = true
			data = data[j+1:]
		}
		if readErr == io.EOF {
			break
		}
		if readErr != nil {
			return fmt.Errorf("read genome FASTA: %w", readErr)
		}
	}
	flush()
	return nil
}

// fillFlanks copies the intron windows around each coding exon out of one
// chromosome's sequence. seq is 0-based; exon coordinates are 1-based.
func fillFlanks(seq []byte, transcripts []*Transcript) {
	n := int64(len(seq))
	if n == 0 {
		return
	}
	for _, t := range transcripts {
		for i := range t.Exons {
			e := &t.Exons[i]
			if !e.IsCoding() {
				continue
			}
			// Genomic window immediately 5' and 3' of the exon.
			loStart, loEnd := e.Start-IntronFlankBases, e.Start-1
			hiStart, hiEnd := e.End+1, e.End+IntronFlankBases
			before := clampSlice(seq, loStart, loEnd, n)
			after := clampSlice(seq, hiStart, hiEnd, n)
			if t.Strand != 1 {
				// Coding orientation: the genomic 3' side precedes the exon.
				before, after = ReverseComplementDNA(after), ReverseComplementDNA(before)
			}
			var nb, na int
			e.IntronBeforePacked, nb = PackDNA2Bit(before)
			e.IntronAfterPacked, na = PackDNA2Bit(after)
			e.IntronBeforeLen, e.IntronAfterLen = uint8(nb), uint8(na)
		}
	}
}

// clampSlice returns seq[start..end] (1-based, inclusive) uppercased, clipped to
// the chromosome, or "" when the window falls entirely outside.
func clampSlice(seq []byte, start, end, n int64) string {
	if start < 1 {
		start = 1
	}
	if end > n {
		end = n
	}
	if start > end {
		return ""
	}
	return strings.ToUpper(string(seq[start-1 : end]))
}

// ReverseComplementDNA reverse-complements a DNA string. Unknown bases are
// passed through unchanged so masked regions stay visible rather than silently
// becoming A.
func ReverseComplementDNA(s string) string {
	if s == "" {
		return ""
	}
	out := make([]byte, len(s))
	for i := 0; i < len(s); i++ {
		switch c := s[len(s)-1-i]; c {
		case 'A':
			out[i] = 'T'
		case 'T':
			out[i] = 'A'
		case 'C':
			out[i] = 'G'
		case 'G':
			out[i] = 'C'
		case 'a':
			out[i] = 't'
		case 't':
			out[i] = 'a'
		case 'c':
			out[i] = 'g'
		case 'g':
			out[i] = 'c'
		default:
			out[i] = c
		}
	}
	return string(out)
}
