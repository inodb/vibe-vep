// Package genomeref keeps the reference genome in a form the server can read
// one window at a time.
//
// HGVS names a change as far 3' along the reference as the sequence allows,
// and testing each step reads the bases just past it. A transcript-relative
// name can ask the cached transcript sequence for those bases; a genomic one
// -- hgvsg -- cannot, because the shift runs along the chromosome and the
// variant need not lie in any transcript at all.
//
// Two bits to a base puts GRCh37's primary chromosomes in about 740 MB, which
// is read a window at a time rather than mapped, so the server's own memory
// does not grow. Runs of N are kept as intervals instead of bases: there are
// a few hundred of them and nothing is annotated inside one, but a shift must
// still stop at the edge rather than walk into a run that would decode as A.
package genomeref

import (
	"encoding/binary"
	"encoding/gob"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"
)

// FileName is the index's name inside an assembly directory.
const FileName = "genome_reference.bin"

// magic ends the file, after the index, so neither the index nor the N runs
// need a size bound decided up front.
const magic = "VVGREF01"

// bases maps a two-bit code to its base and back.
const bases = "ACGT"

// NRun is a half-open run of unknown bases, 1-based, [Start, End].
type NRun struct {
	Start, End int64
}

// chromIndex locates one chromosome's packed bases within the file.
type chromIndex struct {
	Name    string
	Length  int64
	DataOff int64
	NRuns   []NRun
}

// Store reads reference bases out of a built index.
type Store struct {
	f      *os.File
	chroms map[string]*chromIndex
}

// Build writes a reference index from a (optionally gzipped) genomic FASTA,
// keeping only the chromosomes want accepts. Names are normalized the way the
// transcript cache normalizes them, so "chr13" and "13" are the same request.
func Build(fastaPath, outPath string, want func(chrom string) bool,
	stream func(path string, want func(string) bool, emit func(string, []byte)) error) error {

	f, err := os.Create(outPath)
	if err != nil {
		return fmt.Errorf("create %s: %w", outPath, err)
	}
	defer f.Close()

	var index []chromIndex
	var off int64
	var writeErr error

	err = stream(fastaPath, want, func(chrom string, seq []byte) {
		if writeErr != nil {
			return
		}
		packed, runs := pack(seq)
		if _, err := f.Write(packed); err != nil {
			writeErr = fmt.Errorf("write %s: %w", chrom, err)
			return
		}
		index = append(index, chromIndex{
			Name:    chrom,
			Length:  int64(len(seq)),
			DataOff: off,
			NRuns:   runs,
		})
		off += int64(len(packed))
	})
	if err != nil {
		return err
	}
	if writeErr != nil {
		return writeErr
	}
	if len(index) == 0 {
		return fmt.Errorf("%s held none of the wanted chromosomes", fastaPath)
	}

	indexOff := off
	if err := gob.NewEncoder(f).Encode(index); err != nil {
		return fmt.Errorf("write index: %w", err)
	}
	var footer [16]byte
	binary.LittleEndian.PutUint64(footer[:8], uint64(indexOff))
	copy(footer[8:], magic)
	if _, err := f.Write(footer[:]); err != nil {
		return fmt.Errorf("write footer: %w", err)
	}
	return f.Sync()
}

// pack encodes one chromosome two bits to a base, collecting the runs of
// anything that is not A, C, G or T as intervals. Those positions are packed
// as A so the stride stays fixed; Range puts the N back.
func pack(seq []byte) ([]byte, []NRun) {
	buf := make([]byte, (len(seq)+3)/4)
	var runs []NRun
	runStart := int64(-1)
	for i := 0; i < len(seq); i++ {
		c := seq[i] &^ 0x20 // uppercase, which is what soft-masking differs by
		code := strings.IndexByte(bases, c)
		if code < 0 {
			if runStart < 0 {
				runStart = int64(i) + 1
			}
			continue
		}
		if runStart >= 0 {
			runs = append(runs, NRun{runStart, int64(i)})
			runStart = -1
		}
		buf[i/4] |= byte(code) << (uint(3-i%4) * 2)
	}
	if runStart >= 0 {
		runs = append(runs, NRun{runStart, int64(len(seq))})
	}
	return buf, runs
}

// Open reads an index's header. The bases themselves stay on disk.
func Open(path string) (*Store, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	size, err := f.Seek(0, io.SeekEnd)
	if err != nil {
		f.Close()
		return nil, err
	}
	if size < 16 {
		f.Close()
		return nil, fmt.Errorf("%s is too short to be a reference index", path)
	}
	var footer [16]byte
	if _, err := f.ReadAt(footer[:], size-16); err != nil {
		f.Close()
		return nil, err
	}
	if string(footer[8:]) != magic {
		f.Close()
		return nil, fmt.Errorf("%s is not a reference index", path)
	}
	indexOff := int64(binary.LittleEndian.Uint64(footer[:8]))
	if indexOff < 0 || indexOff > size-16 {
		f.Close()
		return nil, fmt.Errorf("%s names an index past its own end", path)
	}

	var index []chromIndex
	section := io.NewSectionReader(f, indexOff, size-16-indexOff)
	if err := gob.NewDecoder(section).Decode(&index); err != nil {
		f.Close()
		return nil, fmt.Errorf("read index of %s: %w", path, err)
	}

	s := &Store{f: f, chroms: make(map[string]*chromIndex, len(index))}
	for i := range index {
		s.chroms[index[i].Name] = &index[i]
	}
	return s, nil
}

// Close releases the file.
func (s *Store) Close() error {
	if s == nil || s.f == nil {
		return nil
	}
	return s.f.Close()
}

// Chromosomes reports how many chromosomes the index holds.
func (s *Store) Chromosomes() int {
	if s == nil {
		return 0
	}
	return len(s.chroms)
}

// ErrNoChromosome says the index does not hold the chromosome asked for,
// which is the normal answer for a scaffold or a patch.
var ErrNoChromosome = errors.New("chromosome not in reference index")

// Range returns the reference bases from start to end, 1-based and inclusive,
// upper case, with unknown bases as N. A window reaching past either end of
// the chromosome is clipped; one entirely outside it comes back empty.
func (s *Store) Range(chrom string, start, end int64) ([]byte, error) {
	if s == nil || s.f == nil {
		return nil, ErrNoChromosome
	}
	ci, ok := s.chroms[chrom]
	if !ok {
		if ci, ok = s.chroms[strings.TrimPrefix(chrom, "chr")]; !ok {
			return nil, fmt.Errorf("%s: %w", chrom, ErrNoChromosome)
		}
	}
	if start < 1 {
		start = 1
	}
	if end > ci.Length {
		end = ci.Length
	}
	if start > end {
		return nil, nil
	}

	first, last := (start-1)/4, (end-1)/4
	raw := make([]byte, last-first+1)
	if _, err := s.f.ReadAt(raw, ci.DataOff+first); err != nil {
		return nil, fmt.Errorf("read %s:%d-%d: %w", chrom, start, end, err)
	}

	out := make([]byte, end-start+1)
	for i := range out {
		p := start - 1 + int64(i) // 0-based position in the chromosome
		b := raw[p/4-first]
		out[i] = bases[(b>>(uint(3-p%4)*2))&3]
	}
	for _, r := range ci.NRuns {
		if r.End < start {
			continue
		}
		if r.Start > end {
			break
		}
		lo, hi := max(r.Start, start), min(r.End, end)
		for p := lo; p <= hi; p++ {
			out[p-start] = 'N'
		}
	}
	return out, nil
}
