package output

import (
	"bufio"
	"io"
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// MAF column headers for VCF→MAF conversion.
var vcf2mafColumns = []string{
	"Hugo_Symbol",
	"Entrez_Gene_Id",
	"Center",
	"NCBI_Build",
	"Chromosome",
	"Start_Position",
	"End_Position",
	"Strand",
	"Variant_Classification",
	"Variant_Type",
	"Reference_Allele",
	"Tumor_Seq_Allele1",
	"Tumor_Seq_Allele2",
	"dbSNP_RS",
	"dbSNP_Val_Status",
	"Tumor_Sample_Barcode",
	"Matched_Norm_Sample_Barcode",
	"HGVSc",
	"HGVSp",
	"HGVSp_Short",
	"Transcript_ID",
	"Exon_Number",
	"Consequence",
	"IMPACT",
	"BIOTYPE",
	"CANONICAL",
	"Protein_position",
	"Amino_acids",
	"Codons",
}

// VCF2MAFWriter converts annotated VCF variants to MAF format.
type VCF2MAFWriter struct {
	w             *bufio.Writer
	assembly      string
	tumorSampleID string
	sources       []annotate.AnnotationSource
	sourceKeys    []string // pre-built Extra map keys for source columns
	headerWritten bool
}

// NewVCF2MAFWriter creates a new VCF→MAF converter.
func NewVCF2MAFWriter(w io.Writer, assembly, tumorSampleID string) *VCF2MAFWriter {
	return &VCF2MAFWriter{
		w:             bufio.NewWriter(w),
		assembly:      assembly,
		tumorSampleID: tumorSampleID,
	}
}

// SetSources registers annotation sources whose columns will be appended.
func (m *VCF2MAFWriter) SetSources(sources []annotate.AnnotationSource) {
	m.sources = sources
	m.sourceKeys = buildSourceKeys(sources)
}

// WriteHeader writes the MAF header line.
func (m *VCF2MAFWriter) WriteHeader() error {
	cols := make([]string, len(vcf2mafColumns))
	copy(cols, vcf2mafColumns)

	// Append source columns
	for _, src := range m.sources {
		for _, col := range src.Columns() {
			cols = append(cols, src.Name()+"_"+col.Name)
		}
	}

	_, err := m.w.WriteString(strings.Join(cols, "\t") + "\n")
	m.headerWritten = true
	return err
}

// WriteRow writes a single MAF row from a VCF variant and its best annotation.
func (m *VCF2MAFWriter) WriteRow(v *vcf.Variant, ann *annotate.Annotation) error {
	ref, alt, start, end := VCFToMAFAlleles(v.Pos, v.Ref, v.Alt)
	variantType := VariantType(ref, alt)

	var b strings.Builder
	b.Grow(512)

	first := true
	writeField := func(s string) {
		if !first {
			b.WriteByte('\t')
		}
		first = false
		b.WriteString(s)
	}

	writeInt := func(n int64) {
		b.WriteByte('\t')
		first = false
		b.WriteString(strconv.FormatInt(n, 10))
	}

	if ann != nil {
		writeField(ann.GeneName) // Hugo_Symbol
	} else {
		writeField("") // Hugo_Symbol
	}
	writeField("")         // Entrez_Gene_Id
	writeField("")         // Center
	writeField(m.assembly) // NCBI_Build
	writeField(v.Chrom)    // Chromosome
	writeInt(start)        // Start_Position
	writeInt(end)          // End_Position
	writeField("+")        // Strand

	if ann != nil {
		writeField(SOToMAFClassification(ann.Consequence, v)) // Variant_Classification
	} else {
		writeField("") // Variant_Classification
	}
	writeField(variantType) // Variant_Type
	writeField(ref)         // Reference_Allele
	writeField(ref)         // Tumor_Seq_Allele1
	writeField(alt)         // Tumor_Seq_Allele2

	// dbSNP RS ID
	if v.ID != "" && v.ID != "." {
		writeField(v.ID)
	} else {
		writeField("")
	}
	writeField("")              // dbSNP_Val_Status
	writeField(m.tumorSampleID) // Tumor_Sample_Barcode
	writeField("")              // Matched_Norm_Sample_Barcode

	if ann != nil {
		writeField(ann.HGVSc)              // HGVSc
		writeField(ann.HGVSp)              // HGVSp
		writeField(hgvspToShort(ann.HGVSp)) // HGVSp_Short
		writeField(ann.TranscriptID)        // Transcript_ID
		writeField(ann.ExonNumber)          // Exon_Number
		writeField(ann.Consequence)         // Consequence
		writeField(ann.Impact)              // IMPACT
		writeField(ann.Biotype)             // BIOTYPE
		if ann.IsCanonical {
			writeField("YES")
		} else {
			writeField("")
		}
		if ann.ProteinPosition > 0 {
			writeInt(ann.ProteinPosition)
		} else {
			writeField("")
		}
		writeField(ann.AminoAcidChange) // Amino_acids
		writeField(ann.CodonChange)     // Codons
	} else {
		for range 12 {
			writeField("")
		}
	}

	// Append source columns using pre-built keys
	for _, key := range m.sourceKeys {
		if ann != nil {
			writeField(ann.GetExtraKey(key))
		} else {
			writeField("")
		}
	}

	b.WriteByte('\n')
	_, err := m.w.WriteString(b.String())
	return err
}

// Flush flushes any buffered data.
func (m *VCF2MAFWriter) Flush() error {
	return m.w.Flush()
}

// VCFToMAFAlleles converts VCF-convention alleles/position to MAF convention.
// VCF uses a shared prefix base for indels; MAF strips the prefix and adjusts positions.
func VCFToMAFAlleles(vcfPos int64, vcfRef, vcfAlt string) (ref, alt string, start, end int64) {
	// SNV or MNV: no prefix stripping needed
	if len(vcfRef) == len(vcfAlt) {
		return vcfRef, vcfAlt, vcfPos, vcfPos + int64(len(vcfRef)) - 1
	}

	// Find shared prefix length
	prefixLen := 0
	minLen := len(vcfRef)
	if len(vcfAlt) < minLen {
		minLen = len(vcfAlt)
	}
	for i := 0; i < minLen; i++ {
		if vcfRef[i] != vcfAlt[i] {
			break
		}
		prefixLen++
	}

	ref = vcfRef[prefixLen:]
	alt = vcfAlt[prefixLen:]

	if len(ref) == 0 {
		// Insertion: ref is empty after stripping prefix
		ref = "-"
		// MAF insertion: start = last base of prefix, end = start + 1
		start = vcfPos + int64(prefixLen) - 1
		end = start + 1
	} else if len(alt) == 0 {
		// Deletion: alt is empty after stripping prefix
		alt = "-"
		start = vcfPos + int64(prefixLen)
		end = start + int64(len(ref)) - 1
	} else {
		// Complex indel (delins)
		start = vcfPos + int64(prefixLen)
		end = start + int64(len(ref)) - 1
	}

	return ref, alt, start, end
}

// VariantType returns the MAF variant type: SNP, DNP, TNP, ONP, INS, or DEL.
func VariantType(ref, alt string) string {
	if ref == "-" {
		return "INS"
	}
	if alt == "-" {
		return "DEL"
	}
	switch len(ref) {
	case 1:
		if len(alt) == 1 {
			return "SNP"
		}
	case 2:
		if len(alt) == 2 {
			return "DNP"
		}
	case 3:
		if len(alt) == 3 {
			return "TNP"
		}
	}
	if len(ref) == len(alt) {
		return "ONP"
	}
	// Complex cases
	if len(ref) > len(alt) {
		return "DEL"
	}
	return "INS"
}
