// Package cache provides VEP cache loading functionality.
package cache

import (
	"sort"
	"strings"
)

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
	ProteinID       string // Ensembl protein ID (e.g., ENSP00000493376), empty if non-coding
	HGNCId          string // HGNC identifier (e.g., HGNC:14825), empty if not available
	EntrezGeneID    string // NCBI Entrez gene ID (e.g., "673"), empty if not available
	RefSeqIDs       []string // RefSeq mRNA accessions (versioned), in GENCODE metadata order
	Chrom           string // Chromosome
	Start           int64  // Transcript start (1-based)
	End             int64  // Transcript end (1-based, inclusive)
	Strand          int8   // +1 or -1
	Biotype         string // Transcript biotype
	IsCanonicalMSK     bool // MSK canonical transcript
	IsCanonicalEnsembl bool // Ensembl canonical transcript (from GTF tag)
	IsMANESelect    bool   // MANE Select transcript
	Exons           []Exon // Exons sorted ascending by genomic Start
	CDSStart        int64  // CDS start (genomic, 1-based), 0 if non-coding
	CDSEnd          int64  // CDS end (genomic, 1-based), 0 if non-coding
	CDSSequence     string // Coding DNA sequence (loaded on demand)
	// UTR5Packed holds up to UTR5TailBases of 5'UTR immediately preceding the
	// CDS, 2-bit packed like the intron flanks. A deletion that starts
	// upstream of the start codon needs these bases to 3'-shift; without
	// them c.-28_12del cannot become c.-27_13del.
	UTR5Packed      []byte
	UTR5Len         uint16
	// CDSStartOffset is how many bases of the first codon are missing from a
	// 5'-incomplete CDS (GENCODE cds_start_NF): 0, 1 or 2. VEP numbers c.1 from
	// the first base of that notional complete codon, so CDSSequence is padded
	// by this many bases and CDS positions are shifted to match.
	CDSStartOffset  int
	UTR3Sequence    string // 3'UTR sequence immediately following CDSSequence (for stop scanning)
	ProteinLength   int    // Protein length in amino acids (persisted, unlike ProteinSequence)
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
	// IntronBeforePacked/IntronAfterPacked hold up to IntronFlankBases of the
	// adjacent intron in CODING orientation, empty for non-coding exons. They
	// exist so a deletion at an exon edge can still be 3'-shifted past the
	// boundary; see the note on IntronFlankBases.
	//
	// Stored two bits per base, four to a byte, because these are carried for
	// every coding exon of every transcript: 1.4M flanks, where plain bytes
	// cost four times as much cache. The lengths are separate since the final
	// byte is usually partial.
	IntronBeforePacked []byte
	IntronAfterPacked  []byte
	IntronBeforeLen    uint8
	IntronAfterLen     uint8
}

// UTR5TailBases is how much 5'UTR is retained. Measured against VEP111, a
// deletion reaching upstream of the start codon starts at most 330 bases out.
const UTR5TailBases = 500

// UTR5BaseAt returns the base at c.-k (k >= 1), the kth base upstream of the
// start codon.
func (t *Transcript) UTR5BaseAt(k int) (byte, bool) {
	// The stored tail ends at c.-1, so c.-k sits that many from the end.
	return unpackBase(t.UTR5Packed, int(t.UTR5Len), int(t.UTR5Len)-k)
}

// packedBases is the 2-bit alphabet; anything else (an N in an assembly gap)
// cannot be encoded and truncates the flank, which is the right behaviour
// anyway: a shift cannot be validated across a base we do not know.
const packedBases = "ACGT"

// PackDNA2Bit encodes a DNA string two bits per base, stopping at the first
// base outside ACGT. It returns the packed bytes and how many bases they hold.
func PackDNA2Bit(s string) ([]byte, int) {
	n := 0
	for n < len(s) {
		if strings.IndexByte(packedBases, s[n]) < 0 {
			break
		}
		n++
	}
	if n == 0 {
		return nil, 0
	}
	buf := make([]byte, (n+3)/4)
	for i := 0; i < n; i++ {
		code := byte(strings.IndexByte(packedBases, s[i]))
		buf[i/4] |= code << (uint(3-i%4) * 2)
	}
	return buf, n
}

// IntronAfterBase returns the i-th base (0-based) of the intron flank 3' of the
// exon in coding orientation.
func (e *Exon) IntronAfterBase(i int) (byte, bool) {
	return unpackBase(e.IntronAfterPacked, int(e.IntronAfterLen), i)
}

// IntronBeforeBase returns the i-th base (0-based) of the intron flank 5' of
// the exon in coding orientation.
func (e *Exon) IntronBeforeBase(i int) (byte, bool) {
	return unpackBase(e.IntronBeforePacked, int(e.IntronBeforeLen), i)
}

func unpackBase(buf []byte, n, i int) (byte, bool) {
	if i < 0 || i >= n {
		return 0, false
	}
	return packedBases[(buf[i/4]>>(uint(3-i%4)*2))&3], true
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
