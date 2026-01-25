// Package cache provides VEP cache loading functionality.
package cache

// Gene represents a genomic region with associated transcripts.
type Gene struct {
	ID          string        // Gene identifier (e.g., ENSG00000133703)
	Name        string        // Gene symbol (e.g., KRAS)
	Chrom       string        // Chromosome
	Start       int64         // Gene start position (1-based)
	End         int64         // Gene end position (1-based, inclusive)
	Strand      int8          // +1 (forward) or -1 (reverse)
	Biotype     string        // Gene biotype (e.g., protein_coding)
	Transcripts []*Transcript // Associated transcripts
}

// IsForwardStrand returns true if the gene is on the forward strand.
func (g *Gene) IsForwardStrand() bool {
	return g.Strand == 1
}

// IsReverseStrand returns true if the gene is on the reverse strand.
func (g *Gene) IsReverseStrand() bool {
	return g.Strand == -1
}

// Contains returns true if the given position is within the gene boundaries.
func (g *Gene) Contains(pos int64) bool {
	return pos >= g.Start && pos <= g.End
}
