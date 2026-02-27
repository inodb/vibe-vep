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
type MAFWriter struct {
	w            *bufio.Writer
	headerLine   string
	columns      maf.ColumnIndices
	extraColumns []string // additional column names appended after original columns
}

// NewMAFWriter creates a new MAF writer that preserves all original columns.
func NewMAFWriter(w io.Writer, headerLine string, columns maf.ColumnIndices) *MAFWriter {
	return &MAFWriter{
		w:          bufio.NewWriter(w),
		headerLine: headerLine,
		columns:    columns,
	}
}

// AddExtraColumn registers an additional column to append after the original MAF columns.
func (m *MAFWriter) AddExtraColumn(name string) {
	m.extraColumns = append(m.extraColumns, name)
}

// WriteHeader writes the original MAF header line plus any extra columns.
func (m *MAFWriter) WriteHeader() error {
	header := m.headerLine
	for _, col := range m.extraColumns {
		header += "\t" + col
	}
	_, err := m.w.WriteString(header + "\n")
	return err
}

// WriteRow writes a MAF row, updating annotation columns with VEP predictions.
// Original values are preserved when the VEP annotation has no prediction.
func (m *MAFWriter) WriteRow(rawFields []string, ann *annotate.Annotation, v *vcf.Variant) error {
	row := make([]string, len(rawFields))
	copy(row, rawFields)

	if ann != nil {
		m.updateField(row, m.columns.HugoSymbol, ann.GeneName)
		m.updateField(row, m.columns.Consequence, ann.Consequence)
		m.updateField(row, m.columns.TranscriptID, ann.TranscriptID)
		m.updateField(row, m.columns.HGVSc, ann.HGVSc)
		m.updateField(row, m.columns.HGVSp, ann.HGVSp)
		m.updateField(row, m.columns.HGVSpShort, hgvspToShort(ann.HGVSp))
		m.updateField(row, m.columns.VariantClassification, SOToMAFClassification(ann.Consequence, v))
	}

	// Append extra column values
	for _, col := range m.extraColumns {
		val := ""
		if ann != nil {
			switch col {
			case "Gene_Type":
				val = ann.GeneType
			}
		}
		row = append(row, val)
	}

	_, err := m.w.WriteString(strings.Join(row, "\t") + "\n")
	return err
}

// Flush flushes any buffered data to the underlying writer.
func (m *MAFWriter) Flush() error {
	return m.w.Flush()
}

// updateField sets row[idx] to value only if idx is valid and value is non-empty.
func (m *MAFWriter) updateField(row []string, idx int, value string) {
	if idx >= 0 && idx < len(row) && value != "" {
		row[idx] = value
	}
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
