package cache

import "sort"

// IntervalTree provides O(log n + k) overlap queries using a sorted-slice approach.
// Transcripts are loaded once and never modified after build.
type IntervalTree struct {
	intervals []interval
	maxEnd    []int64 // maxEnd[i] = max(End) for intervals[i:]
}

type interval struct {
	start      int64
	end        int64
	transcript *Transcript
}

// BuildIntervalTree creates an interval tree from a slice of transcripts.
func BuildIntervalTree(transcripts []*Transcript) *IntervalTree {
	if len(transcripts) == 0 {
		return &IntervalTree{}
	}

	intervals := make([]interval, len(transcripts))
	for i, t := range transcripts {
		intervals[i] = interval{start: t.Start, end: t.End, transcript: t}
	}

	sort.Slice(intervals, func(i, j int) bool {
		return intervals[i].start < intervals[j].start
	})

	// Build suffix-max array: maxEnd[i] = max(end) for intervals[i:]
	maxEnd := make([]int64, len(intervals))
	maxEnd[len(intervals)-1] = intervals[len(intervals)-1].end
	for i := len(intervals) - 2; i >= 0; i-- {
		maxEnd[i] = intervals[i].end
		if maxEnd[i+1] > maxEnd[i] {
			maxEnd[i] = maxEnd[i+1]
		}
	}

	return &IntervalTree{intervals: intervals, maxEnd: maxEnd}
}

// FindOverlaps returns all transcripts whose [Start, End] range contains pos.
func (t *IntervalTree) FindOverlaps(pos int64) []*Transcript {
	if len(t.intervals) == 0 {
		return nil
	}

	var result []*Transcript

	// Binary search: find rightmost interval with start <= pos.
	// All candidates must have start <= pos, so we only need to scan
	// from index 0 to that boundary.
	hi := sort.Search(len(t.intervals), func(i int) bool {
		return t.intervals[i].start > pos
	})
	// hi is the first index with start > pos; candidates are [0, hi).

	for i := hi - 1; i >= 0; i-- {
		// Prune: maxEnd[i] is the max end for intervals[i:].
		// If maxEnd[i] < pos, no interval from 0..i can contain pos.
		if t.maxEnd[i] < pos {
			break
		}
		if t.intervals[i].end >= pos {
			result = append(result, t.intervals[i].transcript)
		}
	}

	return result
}
