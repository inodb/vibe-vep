// Package cache provides VEP cache loading functionality.
package cache

import "sort"

// CDSRegion represents a contiguous coding segment within an exon,
// with its pre-computed cumulative CDS offset for O(log n) lookup.
type CDSRegion struct {
	GenomicStart int64 // exon CDSStart
	GenomicEnd   int64 // exon CDSEnd
	CDSOffset    int64 // cumulative CDS bases before this region (0-based)
}

// Transcript represents a specific gene isoform.
type Transcript struct {
	ID              string // Transcript ID (e.g., ENST00000311936)
	GeneID          string // Parent gene ID
	GeneName        string // Parent gene symbol
	Chrom           string // Chromosome
	Start           int64  // Transcript start (1-based)
	End             int64  // Transcript end (1-based, inclusive)
	Strand          int8   // +1 or -1
	Biotype         string // Transcript biotype
	IsCanonical     bool   // Ensembl canonical flag
	IsMANESelect    bool   // MANE Select transcript
	Exons           []Exon // Exons sorted ascending by genomic Start
	CDSStart        int64  // CDS start (genomic, 1-based), 0 if non-coding
	CDSEnd          int64  // CDS end (genomic, 1-based), 0 if non-coding
	CDSSequence     string // Coding DNA sequence (loaded on demand)
	UTR3Sequence    string // 3'UTR sequence immediately following CDSSequence (for stop scanning)
	ProteinSequence string // Translated protein sequence (loaded on demand)
	CDSRegions      []CDSRegion // Pre-computed CDS regions sorted ascending by GenomicStart
	ExonCumBases    []int64     // Cumulative exonic bases before each exon (transcript order)
}

// Exon represents a single exon within a transcript.
type Exon struct {
	Number   int   // Exon number (1-based, biological transcript order)
	Start    int64 // Genomic start (1-based)
	End      int64 // Genomic end (1-based, inclusive)
	CDSStart int64 // CDS portion start, 0 if entirely non-coding
	CDSEnd   int64 // CDS portion end, 0 if entirely non-coding
	Frame    int   // Reading frame (0, 1, or 2), -1 if non-coding
}

// BuildCDSIndex pre-computes CDS region offsets and exonic base counts
// for O(log n) position lookups. Sorts exons ascending by genomic Start.
// Must be called after exons are loaded.
func (t *Transcript) BuildCDSIndex() {
	n := len(t.Exons)
	if n == 0 {
		return
	}

	// Sort exons by ascending genomic Start (Exon.Number preserves transcript order).
	sort.Slice(t.Exons, func(i, j int) bool {
		return t.Exons[i].Start < t.Exons[j].Start
	})

	// Build CDS regions for protein-coding transcripts.
	if t.IsProteinCoding() {
		var regions []CDSRegion
		if t.Strand == 1 {
			// Forward strand: CDS pos 1 at lowest genomic coordinate.
			var cum int64
			for i := range t.Exons {
				e := &t.Exons[i]
				if !e.IsCoding() {
					continue
				}
				regions = append(regions, CDSRegion{
					GenomicStart: e.CDSStart,
					GenomicEnd:   e.CDSEnd,
					CDSOffset:    cum,
				})
				cum += e.CDSEnd - e.CDSStart + 1
			}
		} else {
			// Reverse strand: CDS pos 1 at highest genomic coordinate.
			// Iterate from highest to lowest to assign ascending CDSOffset.
			var cum int64
			for i := n - 1; i >= 0; i-- {
				e := &t.Exons[i]
				if !e.IsCoding() {
					continue
				}
				regions = append(regions, CDSRegion{
					GenomicStart: e.CDSStart,
					GenomicEnd:   e.CDSEnd,
					CDSOffset:    cum,
				})
				cum += e.CDSEnd - e.CDSStart + 1
			}
			// Re-sort ascending by GenomicStart for binary search.
			sort.Slice(regions, func(i, j int) bool {
				return regions[i].GenomicStart < regions[j].GenomicStart
			})
		}
		t.CDSRegions = regions
	}

	// Build exon cumulative base counts (transcript order).
	t.ExonCumBases = make([]int64, n)
	if t.Strand == 1 {
		// Forward: transcript 5' is lowest genomic coordinate.
		var cum int64
		for i := 0; i < n; i++ {
			t.ExonCumBases[i] = cum
			cum += t.Exons[i].End - t.Exons[i].Start + 1
		}
	} else {
		// Reverse: transcript 5' is highest genomic coordinate.
		var cum int64
		for i := n - 1; i >= 0; i-- {
			t.ExonCumBases[i] = cum
			cum += t.Exons[i].End - t.Exons[i].Start + 1
		}
	}
}

// IsProteinCoding returns true if the transcript has a coding sequence.
// This includes protein_coding, nonsense_mediated_decay, IG/TR gene segments,
// protein_coding_LoF, and any other biotype with CDS features in GENCODE.
func (t *Transcript) IsProteinCoding() bool {
	return t.CDSStart > 0 && t.CDSEnd > 0
}

// IsForwardStrand returns true if the transcript is on the forward strand.
func (t *Transcript) IsForwardStrand() bool {
	return t.Strand == 1
}

// IsReverseStrand returns true if the transcript is on the reverse strand.
func (t *Transcript) IsReverseStrand() bool {
	return t.Strand == -1
}

// Contains returns true if the given position is within the transcript boundaries.
func (t *Transcript) Contains(pos int64) bool {
	return pos >= t.Start && pos <= t.End
}

// ContainsCDS returns true if the given position is within the CDS boundaries.
func (t *Transcript) ContainsCDS(pos int64) bool {
	if !t.IsProteinCoding() {
		return false
	}
	return pos >= t.CDSStart && pos <= t.CDSEnd
}

// FindExonIdx returns the index of the exon containing the given genomic position,
// or -1 if not in an exon. Uses binary search. Handles both ascending (post-BuildCDSIndex)
// and descending (legacy) exon ordering.
func (t *Transcript) FindExonIdx(pos int64) int {
	n := len(t.Exons)
	if n == 0 {
		return -1
	}
	ascending := n < 2 || t.Exons[0].Start <= t.Exons[n-1].Start
	lo, hi := 0, n-1
	for lo <= hi {
		mid := lo + (hi-lo)/2
		e := &t.Exons[mid]
		if pos >= e.Start && pos <= e.End {
			return mid
		}
		if ascending {
			if pos < e.Start {
				hi = mid - 1
			} else {
				lo = mid + 1
			}
		} else {
			if pos > e.End {
				hi = mid - 1
			} else {
				lo = mid + 1
			}
		}
	}
	return -1
}

// FindExon returns the exon containing the given genomic position, or nil if not in an exon.
func (t *Transcript) FindExon(pos int64) *Exon {
	idx := t.FindExonIdx(pos)
	if idx < 0 {
		return nil
	}
	return &t.Exons[idx]
}

// FindNearestExonIdx returns the index of the exon nearest to pos using binary search.
// Returns the index of the exon containing pos, or the nearest exon boundary.
// Handles both ascending (post-BuildCDSIndex) and descending (legacy) exon ordering.
func (t *Transcript) FindNearestExonIdx(pos int64) int {
	n := len(t.Exons)
	if n == 0 {
		return 0
	}
	ascending := n < 2 || t.Exons[0].Start <= t.Exons[n-1].Start
	lo, hi := 0, n-1
	for lo <= hi {
		mid := lo + (hi-lo)/2
		e := &t.Exons[mid]
		if pos >= e.Start && pos <= e.End {
			return mid // inside exon
		}
		if ascending {
			if pos < e.Start {
				hi = mid - 1
			} else {
				lo = mid + 1
			}
		} else {
			if pos > e.End {
				hi = mid - 1
			} else {
				lo = mid + 1
			}
		}
	}
	// pos is between exons. Return the closer one.
	if lo >= n {
		return n - 1
	}
	if hi < 0 {
		return 0
	}
	distHi := pos - t.Exons[hi].End
	if distHi < 0 {
		distHi = -distHi
	}
	distLo := t.Exons[lo].Start - pos
	if distLo < 0 {
		distLo = -distLo
	}
	if distHi <= distLo {
		return hi
	}
	return lo
}

// IsCoding returns true if the exon contains coding sequence.
func (e *Exon) IsCoding() bool {
	return e.CDSStart > 0 && e.CDSEnd > 0
}
