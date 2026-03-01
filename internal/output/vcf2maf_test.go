package output

import (
	"bytes"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

func TestVCFToMAFAlleles(t *testing.T) {
	tests := []struct {
		name                     string
		vcfPos                   int64
		vcfRef, vcfAlt           string
		wantRef, wantAlt         string
		wantStart, wantEnd       int64
	}{
		{
			name:   "SNV",
			vcfPos: 100, vcfRef: "C", vcfAlt: "A",
			wantRef: "C", wantAlt: "A", wantStart: 100, wantEnd: 100,
		},
		{
			name:   "DNP",
			vcfPos: 100, vcfRef: "CG", vcfAlt: "AT",
			wantRef: "CG", wantAlt: "AT", wantStart: 100, wantEnd: 101,
		},
		{
			name:   "deletion",
			vcfPos: 100, vcfRef: "ACG", vcfAlt: "A",
			wantRef: "CG", wantAlt: "-", wantStart: 101, wantEnd: 102,
		},
		{
			name:   "insertion",
			vcfPos: 100, vcfRef: "A", vcfAlt: "ATG",
			wantRef: "-", wantAlt: "TG", wantStart: 100, wantEnd: 101,
		},
		{
			name:   "complex indel",
			vcfPos: 100, vcfRef: "ACG", vcfAlt: "AT",
			wantRef: "CG", wantAlt: "T", wantStart: 101, wantEnd: 102,
		},
		{
			name:   "single base deletion",
			vcfPos: 100, vcfRef: "AT", vcfAlt: "A",
			wantRef: "T", wantAlt: "-", wantStart: 101, wantEnd: 101,
		},
		{
			name:   "single base insertion",
			vcfPos: 100, vcfRef: "A", vcfAlt: "AT",
			wantRef: "-", wantAlt: "T", wantStart: 100, wantEnd: 101,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			ref, alt, start, end := VCFToMAFAlleles(tt.vcfPos, tt.vcfRef, tt.vcfAlt)
			assert.Equal(t, tt.wantRef, ref, "ref")
			assert.Equal(t, tt.wantAlt, alt, "alt")
			assert.Equal(t, tt.wantStart, start, "start")
			assert.Equal(t, tt.wantEnd, end, "end")
		})
	}
}

func TestVariantType(t *testing.T) {
	tests := []struct {
		ref, alt string
		want     string
	}{
		{"C", "A", "SNP"},
		{"CG", "AT", "DNP"},
		{"CGA", "ATT", "TNP"},
		{"CGAT", "ATTG", "ONP"},
		{"-", "TG", "INS"},
		{"CG", "-", "DEL"},
	}
	for _, tt := range tests {
		t.Run(tt.want, func(t *testing.T) {
			assert.Equal(t, tt.want, VariantType(tt.ref, tt.alt))
		})
	}
}

func TestVCF2MAFWriter_FullRow(t *testing.T) {
	var buf bytes.Buffer
	w := NewVCF2MAFWriter(&buf, "GRCh38", "TUMOR_01")
	require.NoError(t, w.WriteHeader())

	v := &vcf.Variant{
		Chrom:  "12",
		Pos:    25245351,
		ID:     "rs121913529",
		Ref:    "C",
		Alt:    "A",
		Qual:   100,
		Filter: "PASS",
		Info:   map[string]interface{}{},
	}

	ann := &annotate.Annotation{
		GeneName:        "KRAS",
		Consequence:     "missense_variant",
		Impact:          "MODERATE",
		TranscriptID:    "ENST00000311936",
		Biotype:         "protein_coding",
		IsCanonical:     true,
		ExonNumber:      "2/5",
		HGVSc:           "c.34G>T",
		HGVSp:           "p.Gly12Cys",
		ProteinPosition: 12,
		AminoAcidChange: "G/C",
		CodonChange:     "Ggt/Tgt",
		Allele:          "A",
	}

	require.NoError(t, w.WriteRow(v, ann))
	require.NoError(t, w.Flush())

	out := buf.String()
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	require.Len(t, lines, 2) // header + 1 data row

	// Check header
	assert.True(t, strings.HasPrefix(lines[0], "Hugo_Symbol\t"))

	// Parse data row
	fields := strings.Split(lines[1], "\t")
	assert.Equal(t, "KRAS", fields[0])             // Hugo_Symbol
	assert.Equal(t, "GRCh38", fields[3])            // NCBI_Build
	assert.Equal(t, "12", fields[4])                // Chromosome
	assert.Equal(t, "25245351", fields[5])          // Start_Position
	assert.Equal(t, "25245351", fields[6])          // End_Position
	assert.Equal(t, "+", fields[7])                 // Strand
	assert.Equal(t, "Missense_Mutation", fields[8]) // Variant_Classification
	assert.Equal(t, "SNP", fields[9])               // Variant_Type
	assert.Equal(t, "C", fields[10])                // Reference_Allele
	assert.Equal(t, "C", fields[11])                // Tumor_Seq_Allele1
	assert.Equal(t, "A", fields[12])                // Tumor_Seq_Allele2
	assert.Equal(t, "rs121913529", fields[13])      // dbSNP_RS
	assert.Equal(t, "TUMOR_01", fields[15])         // Tumor_Sample_Barcode
	assert.Equal(t, "c.34G>T", fields[17])          // HGVSc
	assert.Equal(t, "p.Gly12Cys", fields[18])       // HGVSp
	assert.Equal(t, "p.G12C", fields[19])           // HGVSp_Short
	assert.Equal(t, "ENST00000311936", fields[20])  // Transcript_ID
	assert.Equal(t, "2/5", fields[21])              // Exon_Number
	assert.Equal(t, "missense_variant", fields[22]) // Consequence
	assert.Equal(t, "MODERATE", fields[23])          // IMPACT
	assert.Equal(t, "protein_coding", fields[24])    // BIOTYPE
	assert.Equal(t, "YES", fields[25])               // CANONICAL
	assert.Equal(t, "12", fields[26])                // Protein_position
}

func TestVCF2MAFWriter_Deletion(t *testing.T) {
	var buf bytes.Buffer
	w := NewVCF2MAFWriter(&buf, "GRCh38", "SAMPLE")
	require.NoError(t, w.WriteHeader())

	v := &vcf.Variant{
		Chrom: "1", Pos: 100, ID: ".", Ref: "ACG", Alt: "A",
		Info: map[string]interface{}{},
	}
	ann := &annotate.Annotation{
		GeneName:    "TEST",
		Consequence: "frameshift_variant",
		Impact:      "HIGH",
	}

	require.NoError(t, w.WriteRow(v, ann))
	require.NoError(t, w.Flush())

	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	fields := strings.Split(lines[1], "\t")
	assert.Equal(t, "101", fields[5])           // Start
	assert.Equal(t, "102", fields[6])           // End
	assert.Equal(t, "DEL", fields[9])           // Variant_Type
	assert.Equal(t, "CG", fields[10])           // Reference_Allele
	assert.Equal(t, "-", fields[12])            // Tumor_Seq_Allele2
	assert.Equal(t, "Frame_Shift_Del", fields[8]) // Variant_Classification
}

func TestVCF2MAFWriter_Insertion(t *testing.T) {
	var buf bytes.Buffer
	w := NewVCF2MAFWriter(&buf, "GRCh38", "SAMPLE")
	require.NoError(t, w.WriteHeader())

	v := &vcf.Variant{
		Chrom: "1", Pos: 100, ID: ".", Ref: "A", Alt: "ATG",
		Info: map[string]interface{}{},
	}
	ann := &annotate.Annotation{
		GeneName:    "TEST",
		Consequence: "frameshift_variant",
		Impact:      "HIGH",
	}

	require.NoError(t, w.WriteRow(v, ann))
	require.NoError(t, w.Flush())

	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	fields := strings.Split(lines[1], "\t")
	assert.Equal(t, "100", fields[5])           // Start
	assert.Equal(t, "101", fields[6])           // End
	assert.Equal(t, "INS", fields[9])           // Variant_Type
	assert.Equal(t, "-", fields[10])            // Reference_Allele
	assert.Equal(t, "TG", fields[12])           // Tumor_Seq_Allele2
	assert.Equal(t, "Frame_Shift_Ins", fields[8]) // Variant_Classification
}

func TestVCF2MAFWriter_NilAnnotation(t *testing.T) {
	var buf bytes.Buffer
	w := NewVCF2MAFWriter(&buf, "GRCh38", "SAMPLE")
	require.NoError(t, w.WriteHeader())

	v := &vcf.Variant{
		Chrom: "1", Pos: 100, ID: ".", Ref: "A", Alt: "T",
		Info: map[string]interface{}{},
	}

	require.NoError(t, w.WriteRow(v, nil))
	require.NoError(t, w.Flush())

	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	fields := strings.Split(lines[1], "\t")
	assert.Equal(t, "", fields[0])  // Hugo_Symbol empty
	assert.Equal(t, "1", fields[4]) // Chromosome
	assert.Equal(t, "SNP", fields[9])
}
