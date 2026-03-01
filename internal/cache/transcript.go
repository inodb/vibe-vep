// Package cache provides VEP cache loading functionality.
package cache

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
	Exons           []Exon // Ordered exons
	CDSStart        int64  // CDS start (genomic, 1-based), 0 if non-coding
	CDSEnd          int64  // CDS end (genomic, 1-based), 0 if non-coding
	CDSSequence     string // Coding DNA sequence (loaded on demand)
	UTR3Sequence    string // 3'UTR sequence immediately following CDSSequence (for stop scanning)
	ProteinSequence string // Translated protein sequence (loaded on demand)
}

// Exon represents a single exon within a transcript.
type Exon struct {
	Number   int   // Exon number (1-based)
	Start    int64 // Genomic start (1-based)
	End      int64 // Genomic end (1-based, inclusive)
	CDSStart int64 // CDS portion start, 0 if entirely non-coding
	CDSEnd   int64 // CDS portion end, 0 if entirely non-coding
	Frame    int   // Reading frame (0, 1, or 2), -1 if non-coding
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

// FindExon returns the exon containing the given genomic position, or nil if not in an exon.
// Uses binary search. Handles both forward-strand (ascending Start) and
// reverse-strand (descending Start) exon ordering.
func (t *Transcript) FindExon(pos int64) *Exon {
	n := len(t.Exons)
	if n == 0 {
		return nil
	}
	// Detect ordering: forward-strand exons are ascending, reverse-strand descending.
	ascending := n < 2 || t.Exons[0].Start <= t.Exons[n-1].Start
	lo, hi := 0, n-1
	for lo <= hi {
		mid := lo + (hi-lo)/2
		e := &t.Exons[mid]
		if pos >= e.Start && pos <= e.End {
			return e
		}
		if ascending {
			if pos < e.Start {
				hi = mid - 1
			} else {
				lo = mid + 1
			}
		} else {
			// Descending: higher Start values come first
			if pos > e.End {
				hi = mid - 1
			} else {
				lo = mid + 1
			}
		}
	}
	return nil
}

// FindNearestExonIdx returns the index of the exon nearest to pos using binary search.
// Returns the index of the exon containing pos, or the nearest exon boundary.
// Handles both ascending (forward-strand) and descending (reverse-strand) exon ordering.
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
