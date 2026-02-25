// Package cache provides VEP cache loading functionality.
package cache

// Transcript represents a specific gene isoform.
type Transcript struct {
	ID              string  // Transcript ID (e.g., ENST00000311936)
	GeneID          string  // Parent gene ID
	GeneName        string  // Parent gene symbol
	Chrom           string  // Chromosome
	Start           int64   // Transcript start (1-based)
	End             int64   // Transcript end (1-based, inclusive)
	Strand          int8    // +1 or -1
	Biotype         string  // Transcript biotype
	IsCanonical     bool    // Ensembl canonical flag
	IsMANESelect    bool    // MANE Select transcript
	Exons           []Exon  // Ordered exons
	CDSStart        int64   // CDS start (genomic, 1-based), 0 if non-coding
	CDSEnd          int64   // CDS end (genomic, 1-based), 0 if non-coding
	CDSSequence     string  // Coding DNA sequence (loaded on demand)
	UTR3Sequence    string  // 3'UTR sequence immediately following CDSSequence (for stop scanning)
	ProteinSequence string  // Translated protein sequence (loaded on demand)
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
func (t *Transcript) FindExon(pos int64) *Exon {
	for i := range t.Exons {
		if pos >= t.Exons[i].Start && pos <= t.Exons[i].End {
			return &t.Exons[i]
		}
	}
	return nil
}

// IsCoding returns true if the exon contains coding sequence.
func (e *Exon) IsCoding() bool {
	return e.CDSStart > 0 && e.CDSEnd > 0
}
