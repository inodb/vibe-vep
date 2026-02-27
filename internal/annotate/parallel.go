package annotate

import (
	"runtime"
	"sync"

	"github.com/inodb/vibe-vep/internal/vcf"
)

// WorkItem holds a parsed variant ready for annotation.
type WorkItem struct {
	Seq     int
	Variant *vcf.Variant
	Extra   any // caller-specific data (e.g. *maf.MAFAnnotation)
}

// WorkResult holds the annotation output for a single variant.
type WorkResult struct {
	Seq     int
	Variant *vcf.Variant
	Anns    []*Annotation
	Err     error
	Extra   any
}

// ParallelAnnotate annotates work items using a pool of workers.
// Results are sent to the returned channel in arrival order (not sequence order).
// Use OrderedCollect to consume results in sequence-number order.
// If workers is 0, runtime.NumCPU() is used.
func (a *Annotator) ParallelAnnotate(items <-chan WorkItem, workers int) <-chan WorkResult {
	if workers <= 0 {
		workers = runtime.NumCPU()
	}

	results := make(chan WorkResult, 2*workers)

	var wg sync.WaitGroup
	wg.Add(workers)

	for range workers {
		go func() {
			defer wg.Done()
			for item := range items {
				anns, err := a.Annotate(item.Variant)
				results <- WorkResult{
					Seq:     item.Seq,
					Variant: item.Variant,
					Anns:    anns,
					Err:     err,
					Extra:   item.Extra,
				}
			}
		}()
	}

	go func() {
		wg.Wait()
		close(results)
	}()

	return results
}

// OrderedCollect calls fn for each result in sequence-number order.
// It buffers out-of-order results in a pending map and emits them
// as soon as the next expected sequence number is available.
// Blocks until the results channel is closed.
func OrderedCollect(results <-chan WorkResult, fn func(WorkResult) error) error {
	pending := make(map[int]WorkResult)
	nextSeq := 0

	for r := range results {
		pending[r.Seq] = r

		for {
			rr, ok := pending[nextSeq]
			if !ok {
				break
			}
			delete(pending, nextSeq)
			nextSeq++
			if err := fn(rr); err != nil {
				// Drain remaining results to unblock workers.
				for range results {
				}
				return err
			}
		}
	}

	return nil
}
