package output

import (
	"bufio"
	"io"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// MAFWriter writes annotated variants in MAF format, preserving all original columns.
// All vibe-vep output is appended as namespaced vibe.* columns — original columns
// are never overwritten.
type MAFWriter struct {
	w          *bufio.Writer
	headerLine string
	columns    maf.ColumnIndices
	sources    []annotate.AnnotationSource
}

// NewMAFWriter creates a new MAF writer that preserves all original columns.
func NewMAFWriter(w io.Writer, headerLine string, columns maf.ColumnIndices) *MAFWriter {
	return &MAFWriter{
		w:          bufio.NewWriter(w),
		headerLine: headerLine,
		columns:    columns,
	}
}

// SetSources registers annotation sources whose columns will be appended.
func (m *MAFWriter) SetSources(sources []annotate.AnnotationSource) {
	m.sources = sources
}

// WriteHeader writes the original MAF header line plus vibe.* columns.
func (m *MAFWriter) WriteHeader() error {
	header := m.headerLine

	// Core prediction columns
	for _, col := range annotate.CoreColumns {
		header += "\tvibe." + col.Name
	}

	// Annotation source columns
	for _, src := range m.sources {
		for _, col := range src.Columns() {
			header += "\tvibe." + src.Name() + "." + col.Name
		}
	}

	_, err := m.w.WriteString(header + "\n")
	return err
}

// WriteRow writes a MAF row with vibe.* namespaced columns appended.
// Original MAF columns are preserved exactly as-is.
func (m *MAFWriter) WriteRow(rawFields []string, ann *annotate.Annotation, v *vcf.Variant) error {
	row := make([]string, len(rawFields), len(rawFields)+7+len(m.sources)*3)
	copy(row, rawFields)

	// Core prediction columns
	if ann != nil {
		row = append(row,
			ann.GeneName,    // vibe.hugo_symbol
			ann.Consequence, // vibe.consequence
			SOToMAFClassification(ann.Consequence, v), // vibe.variant_classification
			ann.TranscriptID,        // vibe.transcript_id
			ann.HGVSc,               // vibe.hgvsc
			ann.HGVSp,               // vibe.hgvsp
			hgvspToShort(ann.HGVSp), // vibe.hgvsp_short
		)
	} else {
		row = append(row, "", "", "", "", "", "", "")
	}

	// Annotation source columns
	for _, src := range m.sources {
		for _, col := range src.Columns() {
			val := ""
			if ann != nil {
				val = ann.GetExtra(src.Name(), col.Name)
			}
			row = append(row, val)
		}
	}

	_, err := m.w.WriteString(strings.Join(row, "\t") + "\n")
	return err
}

// Flush flushes any buffered data to the underlying writer.
func (m *MAFWriter) Flush() error {
	return m.w.Flush()
}

// SOToMAFClassification converts an SO consequence term to a MAF Variant_Classification.
// The variant is used to distinguish Frame_Shift_Del/Ins and In_Frame_Del/Ins.
func SOToMAFClassification(consequence string, v *vcf.Variant) string {
	// Use the first (highest-impact) term if comma-separated
	primary := consequence
	if idx := strings.Index(consequence, ","); idx >= 0 {
		primary = consequence[:idx]
	}
	primary = strings.TrimSpace(primary)

	isDel := v != nil && len(v.Ref) > len(v.Alt)

	switch primary {
	case "missense_variant":
		return "Missense_Mutation"
	case "stop_gained":
		return "Nonsense_Mutation"
	case "synonymous_variant":
		return "Silent"
	case "frameshift_variant":
		if isDel {
			return "Frame_Shift_Del"
		}
		return "Frame_Shift_Ins"
	case "inframe_deletion":
		return "In_Frame_Del"
	case "inframe_insertion":
		return "In_Frame_Ins"
	case "splice_donor_variant", "splice_acceptor_variant":
		return "Splice_Site"
	case "splice_region_variant":
		return "Splice_Region"
	case "stop_lost":
		return "Nonstop_Mutation"
	case "start_lost":
		return "Translation_Start_Site"
	case "3_prime_UTR_variant", "3_prime_utr_variant":
		return "3'UTR"
	case "5_prime_UTR_variant", "5_prime_utr_variant":
		return "5'UTR"
	case "intron_variant":
		return "Intron"
	case "intergenic_variant":
		return "IGR"
	case "downstream_gene_variant":
		return "3'Flank"
	case "upstream_gene_variant":
		return "5'Flank"
	case "non_coding_transcript_exon_variant":
		return "RNA"
	default:
		return primary
	}
}
