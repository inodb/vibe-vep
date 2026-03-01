package output

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// CSQ sub-field names in VEP convention order.
var csqFields = []string{
	"Allele",
	"Consequence",
	"IMPACT",
	"SYMBOL",
	"Gene",
	"Feature_type",
	"Feature",
	"BIOTYPE",
	"EXON",
	"INTRON",
	"HGVSc",
	"HGVSp",
	"cDNA_position",
	"CDS_position",
	"Protein_position",
	"Amino_acids",
	"Codons",
	"CANONICAL",
}

// VCFWriter writes annotations in VCF format with a CSQ INFO field.
// Annotations are buffered per variant and flushed when the variant changes.
type VCFWriter struct {
	w           *bufio.Writer
	headerLines []string // original VCF header lines (## and #CHROM)
	sources     []annotate.AnnotationSource
	sourceKeys  []string // pre-built Extra map keys for source columns

	// Buffered state for the current variant.
	currentChrom string                 // chromosome for grouping
	currentPos   int64                  // position for grouping
	hasVariant   bool                   // whether we have a buffered variant
	currentVars  []*vcf.Variant         // variants seen for this key (may differ in alt)
	annotations  []*annotate.Annotation // buffered annotations
	alts         []string               // unique alt alleles seen
}

// NewVCFWriter creates a new VCF output writer.
func NewVCFWriter(w io.Writer, headerLines []string) *VCFWriter {
	return &VCFWriter{
		w:           bufio.NewWriter(w),
		headerLines: headerLines,
	}
}

// SetSources registers annotation sources whose fields will be included in CSQ.
func (vw *VCFWriter) SetSources(sources []annotate.AnnotationSource) {
	vw.sources = sources
	vw.sourceKeys = buildSourceKeys(sources)
}

// WriteHeader writes the original VCF header lines with an inserted CSQ INFO line.
func (vw *VCFWriter) WriteHeader() error {
	allFields := make([]string, len(csqFields))
	copy(allFields, csqFields)
	for _, src := range vw.sources {
		for _, col := range src.Columns() {
			allFields = append(allFields, src.Name()+"_"+col.Name)
		}
	}

	csqLine := fmt.Sprintf(
		"##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from vibe-vep. Format: %s\">",
		strings.Join(allFields, "|"),
	)

	for _, line := range vw.headerLines {
		if strings.HasPrefix(line, "#CHROM") {
			// Insert CSQ INFO line before #CHROM
			if _, err := vw.w.WriteString(csqLine + "\n"); err != nil {
				return err
			}
		}
		if _, err := vw.w.WriteString(line + "\n"); err != nil {
			return err
		}
	}
	return nil
}

// Write buffers an annotation for the given variant. When a new variant is
// encountered (different chrom/pos), the previous variant's VCF line is flushed.
func (vw *VCFWriter) Write(v *vcf.Variant, ann *annotate.Annotation) error {
	if vw.hasVariant && (vw.currentChrom != v.Chrom || vw.currentPos != v.Pos) {
		if err := vw.flushVariant(); err != nil {
			return err
		}
	}

	if !vw.hasVariant || vw.currentChrom != v.Chrom || vw.currentPos != v.Pos {
		vw.currentChrom = v.Chrom
		vw.currentPos = v.Pos
		vw.hasVariant = true
		vw.currentVars = nil
		vw.annotations = nil
		vw.alts = nil
	}

	vw.currentVars = append(vw.currentVars, v)
	vw.annotations = append(vw.annotations, ann)

	// Track unique alts
	found := false
	for _, a := range vw.alts {
		if a == v.Alt {
			found = true
			break
		}
	}
	if !found {
		vw.alts = append(vw.alts, v.Alt)
	}

	return nil
}

// Flush writes any buffered variant and flushes the underlying writer.
func (vw *VCFWriter) Flush() error {
	if vw.hasVariant {
		if err := vw.flushVariant(); err != nil {
			return err
		}
	}
	return vw.w.Flush()
}

// flushVariant writes the buffered variant as a VCF line with CSQ annotations.
func (vw *VCFWriter) flushVariant() error {
	if len(vw.currentVars) == 0 {
		return nil
	}

	// Use the first variant for base fields
	v := vw.currentVars[0]

	// Reconstruct ALT (may be multi-allelic)
	alt := strings.Join(vw.alts, ",")

	// Reconstruct INFO field from raw string
	info := vw.formatInfo(v.RawInfo)

	// Build CSQ value and append to INFO
	var lb strings.Builder
	lb.Grow(256)

	// Write VCF line fields
	lb.WriteString(v.Chrom)
	lb.WriteByte('\t')
	lb.WriteString(strconv.FormatInt(v.Pos, 10))
	lb.WriteByte('\t')
	lb.WriteString(v.ID)
	lb.WriteByte('\t')
	lb.WriteString(v.Ref)
	lb.WriteByte('\t')
	lb.WriteString(alt)
	lb.WriteByte('\t')
	if v.Qual != 0 {
		lb.WriteString(strconv.FormatFloat(v.Qual, 'g', -1, 64))
	} else {
		lb.WriteByte('.')
	}
	lb.WriteByte('\t')
	lb.WriteString(v.Filter)
	lb.WriteByte('\t')
	if info == "." {
		lb.WriteString("CSQ=")
	} else {
		lb.WriteString(info)
		lb.WriteString(";CSQ=")
	}
	for i, ann := range vw.annotations {
		if i > 0 {
			lb.WriteByte(',')
		}
		vw.writeCSQEntry(&lb, ann)
	}

	// Append FORMAT + sample columns if present
	if v.SampleColumns != "" {
		lb.WriteByte('\t')
		lb.WriteString(v.SampleColumns)
	}

	lb.WriteByte('\n')
	if _, err := vw.w.WriteString(lb.String()); err != nil {
		return err
	}

	// Reset buffer
	vw.hasVariant = false
	vw.currentVars = nil
	vw.annotations = nil
	vw.alts = nil

	return nil
}

// formatInfo strips any existing CSQ field from the raw INFO string.
func (vw *VCFWriter) formatInfo(rawInfo string) string {
	if rawInfo == "" || rawInfo == "." {
		return "."
	}

	// Fast path: no CSQ field present
	if !strings.Contains(rawInfo, "CSQ") {
		return rawInfo
	}

	// Strip CSQ=... from the INFO string
	var b strings.Builder
	for rest := rawInfo; rest != ""; {
		semi := strings.IndexByte(rest, ';')
		var field string
		if semi >= 0 {
			field = rest[:semi]
			rest = rest[semi+1:]
		} else {
			field = rest
			rest = ""
		}
		if strings.HasPrefix(field, "CSQ=") || field == "CSQ" {
			continue
		}
		if b.Len() > 0 {
			b.WriteByte(';')
		}
		b.WriteString(field)
	}

	if b.Len() == 0 {
		return "."
	}
	return b.String()
}

// writeCSQEntry writes a single annotation as a pipe-delimited CSQ entry to a builder.
func (vw *VCFWriter) writeCSQEntry(b *strings.Builder, ann *annotate.Annotation) {
	// Write core fields separated by |
	b.WriteString(ann.Allele)
	b.WriteByte('|')
	b.WriteString(ann.Consequence)
	b.WriteByte('|')
	b.WriteString(ann.Impact)
	b.WriteByte('|')
	b.WriteString(ann.GeneName)
	b.WriteByte('|')
	b.WriteString(ann.GeneID)
	b.WriteByte('|')
	if ann.TranscriptID != "" {
		b.WriteString("Transcript")
	}
	b.WriteByte('|')
	b.WriteString(ann.TranscriptID)
	b.WriteByte('|')
	b.WriteString(ann.Biotype)
	b.WriteByte('|')
	b.WriteString(ann.ExonNumber)
	b.WriteByte('|')
	b.WriteString(ann.IntronNumber)
	b.WriteByte('|')
	b.WriteString(ann.HGVSc)
	b.WriteByte('|')
	b.WriteString(ann.HGVSp)
	b.WriteByte('|')
	if ann.CDNAPosition > 0 {
		b.WriteString(strconv.FormatInt(ann.CDNAPosition, 10))
	}
	b.WriteByte('|')
	if ann.CDSPosition > 0 {
		b.WriteString(strconv.FormatInt(ann.CDSPosition, 10))
	}
	b.WriteByte('|')
	if ann.ProteinPosition > 0 {
		b.WriteString(strconv.FormatInt(ann.ProteinPosition, 10))
	}
	b.WriteByte('|')
	b.WriteString(ann.AminoAcidChange)
	b.WriteByte('|')
	b.WriteString(ann.CodonChange)
	b.WriteByte('|')
	if ann.IsCanonical {
		b.WriteString("YES")
	}

	// Append annotation source fields from Extra map using pre-built keys
	for _, key := range vw.sourceKeys {
		b.WriteByte('|')
		b.WriteString(ann.GetExtraKey(key))
	}
}

// buildSourceKeys pre-computes the Extra map keys for all source columns.
func buildSourceKeys(sources []annotate.AnnotationSource) []string {
	var keys []string
	for _, src := range sources {
		name := src.Name()
		for _, col := range src.Columns() {
			keys = append(keys, name+"."+col.Name)
		}
	}
	return keys
}
