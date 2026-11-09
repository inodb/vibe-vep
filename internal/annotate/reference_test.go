package annotate

import (
	"errors"
	"testing"
)

// countingReference records how often the underlying reference is read, so
// the window can be shown to collapse a walk into one read.
type countingReference struct {
	first int64
	seq   string
	reads int
}

func (c *countingReference) Range(chrom string, start, end int64) ([]byte, error) {
	c.reads++
	if chrom != "1" {
		return nil, errors.New("no such chromosome")
	}
	last := c.first + int64(len(c.seq)) - 1
	if start < c.first {
		start = c.first
	}
	if end > last {
		end = last
	}
	if start > end {
		return nil, nil
	}
	return []byte(c.seq[start-c.first : end-c.first+1]), nil
}

func TestBufferedReferenceReadsOnceForAWalk(t *testing.T) {
	src := &countingReference{first: 1, seq: string(make([]byte, 0))}
	seq := make([]byte, 1000)
	for i := range seq {
		seq[i] = "ACGT"[i%4]
	}
	src.seq = string(seq)

	b := buffered(src)
	for pos := int64(400); pos <= 500; pos++ {
		got, err := b.Range("1", pos, pos)
		if err != nil || len(got) != 1 || got[0] != seq[pos-1] {
			t.Fatalf("Range(1,%d,%d) = %q, %v", pos, pos, got, err)
		}
	}
	if src.reads != 1 {
		t.Errorf("underlying reads = %d, want 1 for a walk inside one window", src.reads)
	}

	// Walking out of the window fetches another, and still answers right.
	if got, _ := b.Range("1", 900, 900); len(got) != 1 || got[0] != seq[899] {
		t.Errorf("Range past the window = %q, want %q", got, seq[899:900])
	}
	if src.reads != 2 {
		t.Errorf("underlying reads = %d, want 2", src.reads)
	}
}

func TestBufferedReferenceAnswersAtTheChromosomeEnd(t *testing.T) {
	src := &countingReference{first: 1, seq: "ACGTACGTAC"}
	b := buffered(src)
	if got, _ := b.Range("1", 9, 10); string(got) != "AC" {
		t.Errorf("Range at the end = %q, want %q", got, "AC")
	}
	if got, _ := b.Range("1", 1, 4); string(got) != "ACGT" {
		t.Errorf("Range at the start = %q, want %q", got, "ACGT")
	}
}

func TestBufferedLeavesNilAlone(t *testing.T) {
	if buffered(nil) != nil {
		t.Error("buffered(nil) should stay nil so callers can still test for one")
	}
}
