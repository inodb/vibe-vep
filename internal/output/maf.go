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
// When replace=false (default), all vibe-vep output is appended as namespaced vibe.* columns.
// When replace=true, core columns are overwritten in-place.
type MAFWriter struct {
	w          *bufio.Writer
	headerLine string
	columns    maf.ColumnIndices
	sources    []annotate.AnnotationSource
	sourceKeys []string // pre-built Extra map keys for source columns
	replace    bool
	allEffects bool // include vibe.all_effects column (default true)
}

// NewMAFWriter creates a new MAF writer that preserves all original columns.
func NewMAFWriter(w io.Writer, headerLine string, columns maf.ColumnIndices) *MAFWriter {
	return &MAFWriter{
		w:          bufio.NewWriter(w),
		headerLine: headerLine,
		columns:    columns,
		allEffects: true,
	}
}

// SetAllEffects enables or disables the vibe.all_effects column.
func (m *MAFWriter) SetAllEffects(enabled bool) {
	m.allEffects = enabled
}

// SetSources registers annotation sources whose columns will be appended.
func (m *MAFWriter) SetSources(sources []annotate.AnnotationSource) {
	m.sources = sources
	m.sourceKeys = buildSourceKeys(sources)
}

// SetReplace enables in-place overwrite mode. When true, core columns
// (Hugo_Symbol, Consequence, Variant_Classification, Transcript_ID, HGVSc,
// HGVSp, HGVSp_Short) are overwritten at their original positions instead of
// being appended as vibe.* columns.
func (m *MAFWriter) SetReplace(replace bool) {
	m.replace = replace
}

// WriteHeader writes the MAF header line.
// In default mode, appends vibe.* core columns + source columns.
// In replace mode, keeps original header unchanged and appends only source columns.
func (m *MAFWriter) WriteHeader() error {
	header := m.headerLine

	if !m.replace {
		// Core prediction columns with vibe. prefix
		for _, col := range annotate.CoreColumns {
			header += "\tvibe." + col.Name
		}
	}

	// all_effects column (before source columns)
	if m.allEffects {
		if m.replace {
			header += "\tall_effects"
		} else {
			header += "\tvibe.all_effects"
		}
	}

	// Annotation source columns
	for _, src := range m.sources {
		name := src.Name()
		for _, col := range src.Columns() {
			if name == "" {
				if m.replace {
					header += "\t" + col.Name
				} else {
					header += "\tvibe." + col.Name
				}
			} else {
				if m.replace {
					header += "\t" + name + "." + col.Name
				} else {
					header += "\tvibe." + name + "." + col.Name
				}
			}
		}
	}

	_, err := m.w.WriteString(header + "\n")
	return err
}

// WriteRow writes a MAF row.
// In default mode, original columns are preserved and vibe.* columns appended.
// In replace mode, core columns are overwritten at their original indices.
func (m *MAFWriter) WriteRow(rawFields []string, ann *annotate.Annotation, allAnns []*annotate.Annotation, v *vcf.Variant) error {
	if m.replace {
		return m.writeRowReplace(rawFields, ann, allAnns, v)
	}
	return m.writeRowAppend(rawFields, ann, allAnns, v)
}

// writeRowAppend writes a row in default (append) mode.
func (m *MAFWriter) writeRowAppend(rawFields []string, ann *annotate.Annotation, allAnns []*annotate.Annotation, v *vcf.Variant) error {
	row := make([]string, len(rawFields), len(rawFields)+10+len(m.sources)*3)
	copy(row, rawFields)

	// Core prediction columns
	if ann != nil {
		canonMSK := ""
		if ann.IsCanonicalMSK {
			canonMSK = "YES"
		}
		canonEns := ""
		if ann.IsCanonicalEnsembl {
			canonEns = "YES"
		}
		row = append(row,
			ann.GeneName,                              // vibe.hugo_symbol
			ann.Consequence,                           // vibe.consequence
			SOToMAFClassification(ann.Consequence, v), // vibe.variant_classification
			ann.TranscriptID,                          // vibe.transcript_id
			ann.HGVSc,                                 // vibe.hgvsc
			ann.HGVSp,                                 // vibe.hgvsp
			HGVSpToShort(ann.HGVSp),                   // vibe.hgvsp_short
			canonMSK,                                   // vibe.canonical_mskcc
			canonEns,                                   // vibe.canonical_ensembl
		)
	} else {
		row = append(row, "", "", "", "", "", "", "", "", "")
	}

	// all_effects column
	if m.allEffects {
		row = append(row, FormatAllEffects(allAnns))
	}

	// Annotation source columns
	for _, key := range m.sourceKeys {
		val := ""
		if ann != nil {
			val = ann.GetExtraKey(key)
		}
		row = append(row, val)
	}

	_, err := m.w.WriteString(strings.Join(row, "\t") + "\n")
	return err
}

// writeRowReplace writes a row in replace mode, overwriting core columns in-place.
func (m *MAFWriter) writeRowReplace(rawFields []string, ann *annotate.Annotation, allAnns []*annotate.Annotation, v *vcf.Variant) error {
	row := make([]string, len(rawFields), len(rawFields)+1+len(m.sources)*3)
	copy(row, rawFields)

	if ann != nil {
		// Overwrite core columns at their original positions
		setIfPresent(row, m.columns.HugoSymbol, ann.GeneName)
		setIfPresent(row, m.columns.Consequence, ann.Consequence)
		setIfPresent(row, m.columns.VariantClassification, SOToMAFClassification(ann.Consequence, v))
		setIfPresent(row, m.columns.TranscriptID, ann.TranscriptID)
		setIfPresent(row, m.columns.HGVSc, ann.HGVSc)
		setIfPresent(row, m.columns.HGVSp, ann.HGVSp)
		setIfPresent(row, m.columns.HGVSpShort, HGVSpToShort(ann.HGVSp))
	}

	// all_effects column
	if m.allEffects {
		row = append(row, FormatAllEffects(allAnns))
	}

	// Annotation source columns appended at end
	for _, key := range m.sourceKeys {
		val := ""
		if ann != nil {
			val = ann.GetExtraKey(key)
		}
		row = append(row, val)
	}

	_, err := m.w.WriteString(strings.Join(row, "\t") + "\n")
	return err
}

// setIfPresent sets row[idx] = val if idx >= 0 and within bounds.
func setIfPresent(row []string, idx int, val string) {
	if idx >= 0 && idx < len(row) {
		row[idx] = val
	}
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
