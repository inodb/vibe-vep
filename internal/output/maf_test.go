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
	// Header should have original columns + 7 vibe.* core columns
	wantPrefix := header + "\tvibe.hugo_symbol\tvibe.consequence\tvibe.variant_classification\tvibe.transcript_id\tvibe.hgvsc\tvibe.hgvsp\tvibe.hgvsp_short\n"
	if got != wantPrefix {
		t.Errorf("header = %q, want %q", got, wantPrefix)
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
	// nil annotation = empty vibe.* columns
	if err := w.WriteRow(fields, nil, nil); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	got := strings.TrimRight(buf.String(), "\n")
	parts := strings.Split(got, "\t")
	// 20 original + 7 vibe.* core columns
	if len(parts) != 27 {
		t.Fatalf("expected 27 columns, got %d", len(parts))
	}
	for i := 0; i < 20; i++ {
		want := fields[i]
		if parts[i] != want {
			t.Errorf("column %d = %q, want %q", i, parts[i], want)
		}
	}
	// vibe.* columns should be empty
	for i := 20; i < 27; i++ {
		if parts[i] != "" {
			t.Errorf("vibe column %d = %q, want empty", i, parts[i])
		}
	}
}

func TestMAFWriter_NamespacedColumns(t *testing.T) {
	fields := []string{
		"OLD_GENE",   // 0: Hugo_Symbol
		"old_conseq", // 1: Consequence
		"p.O1X",      // 2: HGVSp_Short
		"ENST0001",   // 3: Transcript_ID
		"12",         // 4: Chromosome
		"100",        // 5: Start_Position
		"100",        // 6: End_Position
		"C",          // 7: Reference_Allele
		"T",          // 8: Tumor_Seq_Allele2
		"c.1A>T",     // 9: HGVSc
		"Silent",     // 10: Variant_Classification
		"p.Old1New",  // 11: HGVSp
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
	w := NewMAFWriter(&buf, "header", cols)
	if err := w.WriteRow(fields, ann, v); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	parts := strings.Split(strings.TrimRight(buf.String(), "\n"), "\t")

	// Original columns should NOT be modified
	if parts[0] != "OLD_GENE" {
		t.Errorf("original Hugo_Symbol = %q, want OLD_GENE (preserved)", parts[0])
	}
	if parts[1] != "old_conseq" {
		t.Errorf("original Consequence = %q, want old_conseq (preserved)", parts[1])
	}
	if parts[10] != "Silent" {
		t.Errorf("original Variant_Classification = %q, want Silent (preserved)", parts[10])
	}

	// vibe.* columns should have predictions
	vibeStart := len(fields)
	checks := map[int]string{
		vibeStart + 0: "KRAS",              // vibe.hugo_symbol
		vibeStart + 1: "missense_variant",  // vibe.consequence
		vibeStart + 2: "Missense_Mutation", // vibe.variant_classification
		vibeStart + 3: "ENST00000311936",   // vibe.transcript_id
		vibeStart + 4: "c.34G>T",           // vibe.hgvsc
		vibeStart + 5: "p.Gly12Cys",        // vibe.hgvsp
		vibeStart + 6: "p.G12C",            // vibe.hgvsp_short
	}
	for idx, want := range checks {
		if idx >= len(parts) {
			t.Errorf("column %d missing, want %q", idx, want)
			continue
		}
		if parts[idx] != want {
			t.Errorf("column %d = %q, want %q", idx, parts[idx], want)
		}
	}
}

func TestMAFWriter_WithSources(t *testing.T) {
	var buf bytes.Buffer
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

	w := NewMAFWriter(&buf, "Hugo_Symbol\tConsequence", cols)
	w.SetSources([]annotate.AnnotationSource{&testSource{
		name:    "oncokb",
		version: "v1",
		columns: []annotate.ColumnDef{{Name: "gene_type", Description: "Gene type"}},
	}})

	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}

	// Row with extra data
	ann := &annotate.Annotation{
		GeneName:    "KRAS",
		Consequence: "missense_variant",
	}
	ann.SetExtra("oncokb", "gene_type", "ONCOGENE")
	if err := w.WriteRow([]string{"KRAS", "old"}, ann, &vcf.Variant{Ref: "C", Alt: "T"}); err != nil {
		t.Fatal(err)
	}

	// Row without extra data
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

	// Header should have vibe.oncokb.gene_type
	if !strings.HasSuffix(lines[0], "\tvibe.oncokb.gene_type") {
		t.Errorf("header missing vibe.oncokb.gene_type: %s", lines[0])
	}

	// First row should have ONCOGENE in last column
	parts := strings.Split(lines[1], "\t")
	if parts[len(parts)-1] != "ONCOGENE" {
		t.Errorf("vibe.oncokb.gene_type = %q, want ONCOGENE", parts[len(parts)-1])
	}

	// Second row should have empty gene_type
	parts2 := strings.Split(lines[2], "\t")
	if parts2[len(parts2)-1] != "" {
		t.Errorf("vibe.oncokb.gene_type = %q, want empty", parts2[len(parts2)-1])
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

// testSource is a minimal AnnotationSource for testing.
type testSource struct {
	name    string
	version string
	columns []annotate.ColumnDef
}

func (s *testSource) Name() string                                         { return s.name }
func (s *testSource) Version() string                                      { return s.version }
func (s *testSource) Columns() []annotate.ColumnDef                        { return s.columns }
func (s *testSource) Annotate(v *vcf.Variant, anns []*annotate.Annotation) {}
