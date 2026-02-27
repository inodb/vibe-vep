package annotate

import (
	"fmt"
	"testing"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// mockLookup returns no transcripts (intergenic) for all positions.
type mockLookup struct{}

func (m *mockLookup) FindTranscripts(string, int64) []*cache.Transcript { return nil }

func makeItems(n int) <-chan WorkItem {
	ch := make(chan WorkItem, n)
	for i := range n {
		ch <- WorkItem{
			Seq: i,
			Variant: &vcf.Variant{
				Chrom: "1",
				Pos:   int64(100 + i),
				Ref:   "A",
				Alt:   "T",
			},
			Extra: i,
		}
	}
	close(ch)
	return ch
}

func TestParallelAnnotate_OrderPreservation(t *testing.T) {
	ann := NewAnnotator(&mockLookup{})

	items := makeItems(200)
	results := ann.ParallelAnnotate(items, 8)

	var collected []int
	err := OrderedCollect(results, func(r WorkResult) error {
		require.NoError(t, r.Err)
		collected = append(collected, r.Seq)
		return nil
	})
	require.NoError(t, err)

	assert.Len(t, collected, 200)
	for i, seq := range collected {
		assert.Equal(t, i, seq, "result %d out of order", i)
	}
}

func TestParallelAnnotate_SingleWorker(t *testing.T) {
	ann := NewAnnotator(&mockLookup{})

	items := makeItems(50)
	results := ann.ParallelAnnotate(items, 1)

	var collected []int
	err := OrderedCollect(results, func(r WorkResult) error {
		collected = append(collected, r.Seq)
		return nil
	})
	require.NoError(t, err)

	assert.Len(t, collected, 50)
	for i, seq := range collected {
		assert.Equal(t, i, seq)
	}
}

func TestParallelAnnotate_ExtraPreserved(t *testing.T) {
	ann := NewAnnotator(&mockLookup{})

	items := makeItems(10)
	results := ann.ParallelAnnotate(items, 4)

	err := OrderedCollect(results, func(r WorkResult) error {
		// Extra was set to the sequence number in makeItems
		assert.Equal(t, r.Seq, r.Extra.(int))
		return nil
	})
	require.NoError(t, err)
}

func TestParallelAnnotate_EmptyInput(t *testing.T) {
	ann := NewAnnotator(&mockLookup{})

	ch := make(chan WorkItem)
	close(ch)
	results := ann.ParallelAnnotate(ch, 4)

	count := 0
	err := OrderedCollect(results, func(r WorkResult) error {
		count++
		return nil
	})
	require.NoError(t, err)
	assert.Equal(t, 0, count)
}

func TestOrderedCollect_EarlyError(t *testing.T) {
	ann := NewAnnotator(&mockLookup{})

	items := makeItems(100)
	results := ann.ParallelAnnotate(items, 4)

	count := 0
	err := OrderedCollect(results, func(r WorkResult) error {
		count++
		if count == 5 {
			return fmt.Errorf("stop at 5")
		}
		return nil
	})
	require.Error(t, err)
	assert.Equal(t, 5, count)
}

func TestParallelAnnotate_ProducesAnnotations(t *testing.T) {
	ann := NewAnnotator(&mockLookup{})

	items := makeItems(5)
	results := ann.ParallelAnnotate(items, 2)

	err := OrderedCollect(results, func(r WorkResult) error {
		require.NoError(t, r.Err)
		// mockLookup returns no transcripts → intergenic annotation
		require.Len(t, r.Anns, 1)
		assert.Equal(t, ConsequenceIntergenicVariant, r.Anns[0].Consequence)
		return nil
	})
	require.NoError(t, err)
}
