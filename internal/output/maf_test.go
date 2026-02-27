package output

import (
	"bytes"
	"strings"
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/vcf"
)

func TestMAFWriter_Header(t *testing.T) {
	var buf bytes.Buffer
	header := "Hugo_Symbol\tChromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2"
	w := NewMAFWriter(&buf, header, maf.ColumnIndices{})
	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}
	got := buf.String()
	if got != header+"\n" {
		t.Errorf("header = %q, want %q", got, header+"\n")
	}
}

func TestMAFWriter_PreservesAllColumns(t *testing.T) {
	// 20-column input row
	fields := make([]string, 20)
	for i := range fields {
		fields[i] = "orig_" + string(rune('A'+i))
	}

	var buf bytes.Buffer
	cols := maf.ColumnIndices{
		Chromosome:            -1,
		StartPosition:         -1,
		EndPosition:           -1,
		ReferenceAllele:       -1,
		TumorSeqAllele2:       -1,
		HugoSymbol:            -1,
		Consequence:           -1,
		HGVSpShort:            -1,
		TranscriptID:          -1,
		VariantType:           -1,
		NCBIBuild:             -1,
		HGVSc:                 -1,
		VariantClassification: -1,
		HGVSp:                 -1,
	}
	w := NewMAFWriter(&buf, "header", cols)
	// nil annotation = no updates
	if err := w.WriteRow(fields, nil, nil); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	got := strings.TrimRight(buf.String(), "\n")
	parts := strings.Split(got, "\t")
	if len(parts) != 20 {
		t.Fatalf("expected 20 columns, got %d", len(parts))
	}
	for i, p := range parts {
		want := fields[i]
		if p != want {
			t.Errorf("column %d = %q, want %q", i, p, want)
		}
	}
}

func TestMAFWriter_UpdatesAnnotationColumns(t *testing.T) {
	cols := maf.ColumnIndices{
		Chromosome:            4,
		StartPosition:         5,
		EndPosition:           6,
		ReferenceAllele:       7,
		TumorSeqAllele2:       8,
		HugoSymbol:            0,
		Consequence:           1,
		HGVSpShort:            2,
		TranscriptID:          3,
		VariantType:           -1,
		NCBIBuild:             -1,
		HGVSc:                 9,
		VariantClassification: 10,
		HGVSp:                 11,
	}

	fields := []string{
		"OLD_GENE",  // 0: Hugo_Symbol
		"old_conseq", // 1: Consequence
		"p.O1X",     // 2: HGVSp_Short
		"ENST0001",  // 3: Transcript_ID
		"12",        // 4: Chromosome
		"100",       // 5: Start_Position
		"100",       // 6: End_Position
		"C",         // 7: Reference_Allele
		"T",         // 8: Tumor_Seq_Allele2
		"c.1A>T",    // 9: HGVSc
		"Silent",    // 10: Variant_Classification
		"p.Old1New", // 11: HGVSp
	}

	ann := &annotate.Annotation{
		GeneName:     "KRAS",
		Consequence:  "missense_variant",
		TranscriptID: "ENST00000311936",
		HGVSp:        "p.Gly12Cys",
		HGVSc:        "c.34G>T",
	}
	v := &vcf.Variant{Ref: "C", Alt: "T"}

	var buf bytes.Buffer
	w := NewMAFWriter(&buf, "header", cols)
	if err := w.WriteRow(fields, ann, v); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	parts := strings.Split(strings.TrimRight(buf.String(), "\n"), "\t")

	checks := map[int]string{
		0:  "KRAS",
		1:  "missense_variant",
		2:  "p.G12C",
		3:  "ENST00000311936",
		9:  "c.34G>T",
		10: "Missense_Mutation",
		11: "p.Gly12Cys",
	}
	for idx, want := range checks {
		if parts[idx] != want {
			t.Errorf("column %d = %q, want %q", idx, parts[idx], want)
		}
	}
}

func TestMAFWriter_PreservesOriginalWhenEmpty(t *testing.T) {
	cols := maf.ColumnIndices{
		Chromosome:            -1,
		StartPosition:         -1,
		EndPosition:           -1,
		ReferenceAllele:       -1,
		TumorSeqAllele2:       -1,
		HugoSymbol:            0,
		Consequence:           1,
		HGVSpShort:            2,
		TranscriptID:          3,
		VariantType:           -1,
		NCBIBuild:             -1,
		HGVSc:                 4,
		VariantClassification: 5,
		HGVSp:                 6,
	}

	fields := []string{
		"ORIGINAL_GENE",
		"original_consequence",
		"p.O1X",
		"ENST0001",
		"c.1A>T",
		"Missense_Mutation",
		"p.Orig",
	}

	// Intergenic annotation — most fields empty
	ann := &annotate.Annotation{
		Consequence: "intergenic_variant",
	}
	v := &vcf.Variant{Ref: "C", Alt: "T"}

	var buf bytes.Buffer
	w := NewMAFWriter(&buf, "header", cols)
	if err := w.WriteRow(fields, ann, v); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	parts := strings.Split(strings.TrimRight(buf.String(), "\n"), "\t")

	// Gene, HGVSp_Short, Transcript_ID, HGVSc, HGVSp should be preserved (VEP has empty values)
	if parts[0] != "ORIGINAL_GENE" {
		t.Errorf("Hugo_Symbol = %q, want preserved original", parts[0])
	}
	if parts[2] != "p.O1X" {
		t.Errorf("HGVSp_Short = %q, want preserved original", parts[2])
	}
	if parts[3] != "ENST0001" {
		t.Errorf("Transcript_ID = %q, want preserved original", parts[3])
	}
	if parts[4] != "c.1A>T" {
		t.Errorf("HGVSc = %q, want preserved original", parts[4])
	}
	if parts[6] != "p.Orig" {
		t.Errorf("HGVSp = %q, want preserved original", parts[6])
	}
	// Consequence and Variant_Classification should be updated
	if parts[1] != "intergenic_variant" {
		t.Errorf("Consequence = %q, want intergenic_variant", parts[1])
	}
	if parts[5] != "IGR" {
		t.Errorf("Variant_Classification = %q, want IGR", parts[5])
	}
}

func TestMAFWriter_ExtraGeneTypeColumn(t *testing.T) {
	cols := maf.ColumnIndices{
		Chromosome:            -1,
		StartPosition:         -1,
		EndPosition:           -1,
		ReferenceAllele:       -1,
		TumorSeqAllele2:       -1,
		HugoSymbol:            0,
		Consequence:           1,
		HGVSpShort:            -1,
		TranscriptID:          -1,
		VariantType:           -1,
		NCBIBuild:             -1,
		HGVSc:                 -1,
		VariantClassification: -1,
		HGVSp:                 -1,
	}

	var buf bytes.Buffer
	w := NewMAFWriter(&buf, "Hugo_Symbol\tConsequence", cols)
	w.AddExtraColumn("Gene_Type")

	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}

	// Row with GeneType
	ann := &annotate.Annotation{
		GeneName:    "KRAS",
		Consequence: "missense_variant",
		GeneType:    "ONCOGENE",
	}
	if err := w.WriteRow([]string{"KRAS", "old"}, ann, &vcf.Variant{Ref: "C", Alt: "T"}); err != nil {
		t.Fatal(err)
	}

	// Row without GeneType
	ann2 := &annotate.Annotation{
		GeneName:    "UNKNOWN",
		Consequence: "intron_variant",
	}
	if err := w.WriteRow([]string{"UNKNOWN", "old"}, ann2, &vcf.Variant{Ref: "C", Alt: "T"}); err != nil {
		t.Fatal(err)
	}

	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	if len(lines) != 3 {
		t.Fatalf("expected 3 lines (header + 2 rows), got %d", len(lines))
	}

	// Header should have extra column
	if !strings.HasSuffix(lines[0], "\tGene_Type") {
		t.Errorf("header missing Gene_Type column: %s", lines[0])
	}

	// First row should have ONCOGENE
	parts := strings.Split(lines[1], "\t")
	if parts[len(parts)-1] != "ONCOGENE" {
		t.Errorf("Gene_Type = %q, want ONCOGENE", parts[len(parts)-1])
	}

	// Second row should have empty Gene_Type
	parts2 := strings.Split(lines[2], "\t")
	if parts2[len(parts2)-1] != "" {
		t.Errorf("Gene_Type = %q, want empty", parts2[len(parts2)-1])
	}
}

func TestSOToMAFClassification(t *testing.T) {
	tests := []struct {
		consequence string
		ref, alt    string
		want        string
	}{
		{"missense_variant", "C", "T", "Missense_Mutation"},
		{"stop_gained", "C", "T", "Nonsense_Mutation"},
		{"synonymous_variant", "C", "T", "Silent"},
		{"frameshift_variant", "CA", "C", "Frame_Shift_Del"},
		{"frameshift_variant", "C", "CA", "Frame_Shift_Ins"},
		{"inframe_deletion", "CGA", "C", "In_Frame_Del"},
		{"inframe_insertion", "C", "CGAT", "In_Frame_Ins"},
		{"splice_donor_variant", "C", "T", "Splice_Site"},
		{"splice_acceptor_variant", "C", "T", "Splice_Site"},
		{"splice_region_variant", "C", "T", "Splice_Region"},
		{"stop_lost", "C", "T", "Nonstop_Mutation"},
		{"start_lost", "C", "T", "Translation_Start_Site"},
		{"3_prime_UTR_variant", "C", "T", "3'UTR"},
		{"5_prime_UTR_variant", "C", "T", "5'UTR"},
		{"intron_variant", "C", "T", "Intron"},
		{"intergenic_variant", "C", "T", "IGR"},
		{"downstream_gene_variant", "C", "T", "3'Flank"},
		{"upstream_gene_variant", "C", "T", "5'Flank"},
		{"non_coding_transcript_exon_variant", "C", "T", "RNA"},
		// Comma-separated: use first term
		{"missense_variant,splice_region_variant", "C", "T", "Missense_Mutation"},
	}

	for _, tt := range tests {
		v := &vcf.Variant{Ref: tt.ref, Alt: tt.alt}
		got := SOToMAFClassification(tt.consequence, v)
		if got != tt.want {
			t.Errorf("SOToMAFClassification(%q, ref=%s alt=%s) = %q, want %q",
				tt.consequence, tt.ref, tt.alt, got, tt.want)
		}
	}
}
