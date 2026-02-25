package output

import (
	"bufio"
	"fmt"
	"io"
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

	// Buffered state for the current variant.
	currentKey  string                 // "chrom:pos" key for grouping
	currentVars []*vcf.Variant         // variants seen for this key (may differ in alt)
	annotations []*annotate.Annotation // buffered annotations
	alts        []string               // unique alt alleles seen
}

// NewVCFWriter creates a new VCF output writer.
func NewVCFWriter(w io.Writer, headerLines []string) *VCFWriter {
	return &VCFWriter{
		w:           bufio.NewWriter(w),
		headerLines: headerLines,
	}
}

// WriteHeader writes the original VCF header lines with an inserted CSQ INFO line.
func (vw *VCFWriter) WriteHeader() error {
	csqLine := fmt.Sprintf(
		"##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from vibe-vep. Format: %s\">",
		strings.Join(csqFields, "|"),
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
	key := fmt.Sprintf("%s:%d", v.Chrom, v.Pos)

	if vw.currentKey != "" && vw.currentKey != key {
		if err := vw.flushVariant(); err != nil {
			return err
		}
	}

	if vw.currentKey == "" || vw.currentKey != key {
		vw.currentKey = key
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
	if vw.currentKey != "" {
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

	// Format QUAL
	qual := "."
	if v.Qual != 0 {
		qual = fmt.Sprintf("%g", v.Qual)
	}

	// Reconstruct INFO field
	info := vw.formatInfo(v.Info)

	// Build CSQ value
	var csqEntries []string
	for _, ann := range vw.annotations {
		csqEntries = append(csqEntries, formatCSQEntry(ann))
	}
	csqValue := strings.Join(csqEntries, ",")

	// Append CSQ to INFO
	if info == "." {
		info = "CSQ=" + csqValue
	} else {
		info = info + ";CSQ=" + csqValue
	}

	// Write VCF line
	line := fmt.Sprintf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s",
		v.Chrom, v.Pos, v.ID, v.Ref, alt, qual, v.Filter, info)

	if _, err := vw.w.WriteString(line + "\n"); err != nil {
		return err
	}

	// Reset buffer
	vw.currentKey = ""
	vw.currentVars = nil
	vw.annotations = nil
	vw.alts = nil

	return nil
}

// formatInfo reconstructs the INFO field string from a parsed map,
// excluding any existing CSQ field.
func (vw *VCFWriter) formatInfo(info map[string]interface{}) string {
	if len(info) == 0 {
		return "."
	}

	var parts []string
	for k, v := range info {
		if k == "CSQ" {
			continue // Don't carry forward old CSQ
		}
		if b, ok := v.(bool); ok && b {
			parts = append(parts, k)
		} else {
			parts = append(parts, fmt.Sprintf("%s=%v", k, v))
		}
	}

	if len(parts) == 0 {
		return "."
	}
	return strings.Join(parts, ";")
}

// formatCSQEntry formats a single annotation as a pipe-delimited CSQ entry.
func formatCSQEntry(ann *annotate.Annotation) string {
	canonical := ""
	if ann.IsCanonical {
		canonical = "YES"
	}

	featureType := ""
	if ann.TranscriptID != "" {
		featureType = "Transcript"
	}

	cdsPos := ""
	if ann.CDSPosition > 0 {
		cdsPos = fmt.Sprintf("%d", ann.CDSPosition)
	}

	proteinPos := ""
	if ann.ProteinPosition > 0 {
		proteinPos = fmt.Sprintf("%d", ann.ProteinPosition)
	}

	cdnaPos := ""
	if ann.CDNAPosition > 0 {
		cdnaPos = fmt.Sprintf("%d", ann.CDNAPosition)
	}

	fields := []string{
		ann.Allele,
		ann.Consequence,
		ann.Impact,
		ann.GeneName,
		ann.GeneID,
		featureType,
		ann.TranscriptID,
		ann.Biotype,
		ann.ExonNumber,
		ann.IntronNumber,
		ann.HGVSc,
		ann.HGVSp,
		cdnaPos,
		cdsPos,
		proteinPos,
		ann.AminoAcidChange,
		ann.CodonChange,
		canonical,
	}

	return strings.Join(fields, "|")
}
