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
	// Header should have original columns + 9 vibe.* core columns + vibe.all_effects
	wantPrefix := header + "\tvibe.hugo_symbol\tvibe.consequence\tvibe.variant_classification\tvibe.transcript_id\tvibe.hgvsc\tvibe.hgvsp\tvibe.hgvsp_short\tvibe.canonical_mskcc\tvibe.canonical_ensembl\tvibe.all_effects\n"
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
	if err := w.WriteRow(fields, nil, nil, nil); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	got := strings.TrimRight(buf.String(), "\n")
	parts := strings.Split(got, "\t")
	// 20 original + 9 vibe.* core columns + 1 vibe.all_effects
	if len(parts) != 30 {
		t.Fatalf("expected 30 columns, got %d", len(parts))
	}
	for i := 0; i < 20; i++ {
		want := fields[i]
		if parts[i] != want {
			t.Errorf("column %d = %q, want %q", i, parts[i], want)
		}
	}
	// vibe.* columns should be empty
	for i := 20; i < 30; i++ {
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
	if err := w.WriteRow(fields, ann, []*annotate.Annotation{ann}, v); err != nil {
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
	if err := w.WriteRow([]string{"KRAS", "old"}, ann, []*annotate.Annotation{ann}, &vcf.Variant{Ref: "C", Alt: "T"}); err != nil {
		t.Fatal(err)
	}

	// Row without extra data
	ann2 := &annotate.Annotation{
		GeneName:    "UNKNOWN",
		Consequence: "intron_variant",
	}
	if err := w.WriteRow([]string{"UNKNOWN", "old"}, ann2, []*annotate.Annotation{ann2}, &vcf.Variant{Ref: "C", Alt: "T"}); err != nil {
		t.Fatal(err)
	}

	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	if len(lines) != 3 {
		t.Fatalf("expected 3 lines (header + 2 rows), got %d", len(lines))
	}

	// Header should have vibe.all_effects then vibe.oncokb.gene_type
	if !strings.HasSuffix(lines[0], "\tvibe.all_effects\tvibe.oncokb.gene_type") {
		t.Errorf("header missing vibe.all_effects or vibe.oncokb.gene_type: %s", lines[0])
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

func TestMAFWriter_AllEffects_Oncogene(t *testing.T) {
	// KRAS G12C: canonical missense + non-coding intron transcript
	var buf bytes.Buffer
	cols := maf.ColumnIndices{
		HugoSymbol: 0, Consequence: 1,
		Chromosome: -1, StartPosition: -1, EndPosition: -1,
		ReferenceAllele: -1, TumorSeqAllele2: -1,
		HGVSpShort: -1, TranscriptID: -1, VariantType: -1,
		NCBIBuild: -1, HGVSc: -1, VariantClassification: -1, HGVSp: -1,
	}
	w := NewMAFWriter(&buf, "Hugo_Symbol\tConsequence", cols)
	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}

	best := &annotate.Annotation{
		GeneName: "KRAS", Consequence: "missense_variant",
		HGVSp: "p.Gly12Cys", TranscriptID: "ENST00000311936.8",
		HGVSc: "c.34G>T", Impact: "MODERATE", IsCanonicalMSK: true,
	}
	intron := &annotate.Annotation{
		GeneName: "KRAS", Consequence: "intron_variant",
		TranscriptID: "ENST00000556131.1", Impact: "MODIFIER",
	}
	noncoding := &annotate.Annotation{
		GeneName: "KRAS", Consequence: "non_coding_transcript_exon_variant",
		TranscriptID: "ENST00000557334.5", Impact: "MODIFIER",
	}
	allAnns := []*annotate.Annotation{best, intron, noncoding}

	if err := w.WriteRow([]string{"KRAS", "old"}, best, allAnns, &vcf.Variant{Ref: "C", Alt: "A"}); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected 2 lines, got %d", len(lines))
	}

	parts := strings.Split(lines[1], "\t")
	// all_effects is last column (after 2 orig + 9 core)
	allEffectsIdx := len(parts) - 1
	allEffects := parts[allEffectsIdx]

	entries := strings.Split(allEffects, ";")
	if len(entries) != 3 {
		t.Fatalf("expected 3 all_effects entries, got %d: %q", len(entries), allEffects)
	}

	// First entry: canonical KRAS missense
	fields := strings.Split(entries[0], ",")
	if len(fields) != 7 {
		t.Fatalf("entry 0: expected 7 fields, got %d: %q", len(fields), entries[0])
	}
	if fields[0] != "KRAS" {
		t.Errorf("entry 0 gene = %q, want KRAS", fields[0])
	}
	if fields[1] != "missense_variant" {
		t.Errorf("entry 0 consequence = %q, want missense_variant", fields[1])
	}
	if fields[2] != "p.G12C" {
		t.Errorf("entry 0 HGVSp_Short = %q, want p.G12C", fields[2])
	}
	if fields[6] != "YES" {
		t.Errorf("entry 0 canonical_msk = %q, want YES", fields[6])
	}

	// Second entry: intron (no protein change, not canonical)
	fields2 := strings.Split(entries[1], ",")
	if fields2[1] != "intron_variant" {
		t.Errorf("entry 1 consequence = %q, want intron_variant", fields2[1])
	}
	if fields2[2] != "" {
		t.Errorf("entry 1 HGVSp_Short = %q, want empty", fields2[2])
	}
	if fields2[6] != "" {
		t.Errorf("entry 1 canonical_msk = %q, want empty", fields2[6])
	}
}

func TestMAFWriter_AllEffects_TSG(t *testing.T) {
	// TP53 R175H: two coding transcripts (canonical + alternate)
	var buf bytes.Buffer
	cols := maf.ColumnIndices{
		HugoSymbol: 0, Consequence: 1,
		Chromosome: -1, StartPosition: -1, EndPosition: -1,
		ReferenceAllele: -1, TumorSeqAllele2: -1,
		HGVSpShort: -1, TranscriptID: -1, VariantType: -1,
		NCBIBuild: -1, HGVSc: -1, VariantClassification: -1, HGVSp: -1,
	}
	w := NewMAFWriter(&buf, "Hugo_Symbol\tConsequence", cols)
	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}

	canonical := &annotate.Annotation{
		GeneName: "TP53", Consequence: "missense_variant",
		HGVSp: "p.Arg175His", TranscriptID: "ENST00000269305.9",
		HGVSc: "c.524G>A", Impact: "MODERATE", IsCanonicalMSK: true,
	}
	alternate := &annotate.Annotation{
		GeneName: "TP53", Consequence: "missense_variant",
		HGVSp: "p.Arg141His", TranscriptID: "ENST00000610292.4",
		HGVSc: "c.422G>A", Impact: "MODERATE",
	}
	allAnns := []*annotate.Annotation{canonical, alternate}

	if err := w.WriteRow([]string{"TP53", "old"}, canonical, allAnns, &vcf.Variant{Ref: "G", Alt: "A"}); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	parts := strings.Split(lines[1], "\t")
	allEffects := parts[len(parts)-1]

	entries := strings.Split(allEffects, ";")
	if len(entries) != 2 {
		t.Fatalf("expected 2 all_effects entries, got %d: %q", len(entries), allEffects)
	}

	// Canonical transcript should show YES
	f1 := strings.Split(entries[0], ",")
	if f1[0] != "TP53" || f1[2] != "p.R175H" || f1[6] != "YES" {
		t.Errorf("entry 0 = %q, want TP53 with p.R175H and canonical YES", entries[0])
	}

	// Alternate transcript should show empty canonical
	f2 := strings.Split(entries[1], ",")
	if f2[0] != "TP53" || f2[2] != "p.R141H" || f2[6] != "" {
		t.Errorf("entry 1 = %q, want TP53 with p.R141H and no canonical", entries[1])
	}
}

func TestMAFWriter_AllEffects_Replace(t *testing.T) {
	// Verify all_effects works in replace mode
	var buf bytes.Buffer
	cols := maf.ColumnIndices{
		HugoSymbol: 0, Consequence: 1,
		Chromosome: -1, StartPosition: -1, EndPosition: -1,
		ReferenceAllele: -1, TumorSeqAllele2: -1,
		HGVSpShort: -1, TranscriptID: -1, VariantType: -1,
		NCBIBuild: -1, HGVSc: -1, VariantClassification: -1, HGVSp: -1,
	}
	w := NewMAFWriter(&buf, "Hugo_Symbol\tConsequence", cols)
	w.SetReplace(true)
	// allEffects defaults to true

	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}

	ann := &annotate.Annotation{
		GeneName: "BRAF", Consequence: "missense_variant",
		HGVSp: "p.Val600Glu", TranscriptID: "ENST00000288602.11",
		HGVSc: "c.1799T>A", Impact: "MODERATE", IsCanonicalMSK: true,
	}

	if err := w.WriteRow([]string{"BRAF", "old"}, ann, []*annotate.Annotation{ann}, &vcf.Variant{Ref: "T", Alt: "A"}); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")

	// Header should have all_effects appended
	if lines[0] != "Hugo_Symbol\tConsequence\tall_effects" {
		t.Errorf("header = %q, want Hugo_Symbol\\tConsequence\\tall_effects", lines[0])
	}

	// Row: 2 replaced columns + all_effects
	parts := strings.Split(lines[1], "\t")
	if len(parts) != 3 {
		t.Fatalf("expected 3 columns, got %d: %v", len(parts), parts)
	}
	if parts[0] != "BRAF" {
		t.Errorf("Hugo_Symbol = %q, want BRAF", parts[0])
	}
	if !strings.Contains(parts[2], "BRAF,missense_variant,p.V600E") {
		t.Errorf("all_effects = %q, want to contain BRAF,missense_variant,p.V600E", parts[2])
	}
}

func TestMAFWriter_AllEffects_NilAnnotation(t *testing.T) {
	// When annotation fails, all_effects should be empty
	var buf bytes.Buffer
	cols := maf.ColumnIndices{
		HugoSymbol: 0, Consequence: -1,
		Chromosome: -1, StartPosition: -1, EndPosition: -1,
		ReferenceAllele: -1, TumorSeqAllele2: -1,
		HGVSpShort: -1, TranscriptID: -1, VariantType: -1,
		NCBIBuild: -1, HGVSc: -1, VariantClassification: -1, HGVSp: -1,
	}
	w := NewMAFWriter(&buf, "Hugo_Symbol", cols)
	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}

	if err := w.WriteRow([]string{"UNKNOWN"}, nil, nil, nil); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	parts := strings.Split(lines[1], "\t")
	// all_effects should be empty (last column)
	allEffects := parts[len(parts)-1]
	if allEffects != "" {
		t.Errorf("all_effects = %q, want empty for nil annotation", allEffects)
	}
}

func TestMAFWriter_AllEffects_Disabled(t *testing.T) {
	var buf bytes.Buffer
	cols := maf.ColumnIndices{
		HugoSymbol: 0, Consequence: -1,
		Chromosome: -1, StartPosition: -1, EndPosition: -1,
		ReferenceAllele: -1, TumorSeqAllele2: -1,
		HGVSpShort: -1, TranscriptID: -1, VariantType: -1,
		NCBIBuild: -1, HGVSc: -1, VariantClassification: -1, HGVSp: -1,
	}
	w := NewMAFWriter(&buf, "Hugo_Symbol", cols)
	w.SetAllEffects(false)

	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}

	ann := &annotate.Annotation{GeneName: "KRAS", Consequence: "missense_variant"}
	if err := w.WriteRow([]string{"KRAS"}, ann, []*annotate.Annotation{ann}, &vcf.Variant{Ref: "C", Alt: "T"}); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	got := buf.String()
	if strings.Contains(got, "all_effects") {
		t.Error("all_effects column should not appear when disabled")
	}
	// 1 orig + 9 core = 10 columns per row (no all_effects)
	lines := strings.Split(strings.TrimRight(got, "\n"), "\n")
	parts := strings.Split(lines[1], "\t")
	if len(parts) != 10 {
		t.Errorf("expected 10 columns (no all_effects), got %d", len(parts))
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

func TestMAFWriter_Replace(t *testing.T) {
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
	header := "Hugo_Symbol\tConsequence\tHGVSp_Short\tTranscript_ID\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\tHGVSc\tVariant_Classification\tHGVSp"
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
	w := NewMAFWriter(&buf, header, cols)
	w.SetReplace(true)
	w.SetAllEffects(false)

	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteRow(fields, ann, []*annotate.Annotation{ann}, v); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected 2 lines (header + 1 row), got %d", len(lines))
	}

	// Header should be unchanged (no vibe.* core columns)
	if lines[0] != header {
		t.Errorf("header changed in replace mode:\n  got:  %s\n  want: %s", lines[0], header)
	}

	parts := strings.Split(lines[1], "\t")
	if len(parts) != 12 {
		t.Fatalf("expected 12 columns (same as input), got %d", len(parts))
	}

	// Core columns should be OVERWRITTEN
	checks := map[int]string{
		0:  "KRAS",              // Hugo_Symbol
		1:  "missense_variant",  // Consequence
		2:  "p.G12C",            // HGVSp_Short
		3:  "ENST00000311936",   // Transcript_ID
		9:  "c.34G>T",           // HGVSc
		10: "Missense_Mutation", // Variant_Classification
		11: "p.Gly12Cys",        // HGVSp
	}
	for idx, want := range checks {
		if parts[idx] != want {
			t.Errorf("column %d = %q, want %q (overwritten)", idx, parts[idx], want)
		}
	}

	// Non-core columns should be preserved
	if parts[4] != "12" {
		t.Errorf("Chromosome = %q, want 12 (preserved)", parts[4])
	}
}

func TestMAFWriter_ReplaceWithSources(t *testing.T) {
	fields := []string{"GENE", "conseq"}

	ann := &annotate.Annotation{
		GeneName:    "KRAS",
		Consequence: "missense_variant",
	}
	ann.SetExtra("oncokb", "gene_type", "ONCOGENE")

	var buf bytes.Buffer
	header := "Hugo_Symbol\tConsequence"
	cols := maf.ColumnIndices{
		HugoSymbol:            0,
		Consequence:           1,
		HGVSpShort:            -1,
		TranscriptID:          -1,
		HGVSc:                 -1,
		VariantClassification: -1,
		HGVSp:                 -1,
		Chromosome:            -1,
		StartPosition:         -1,
		EndPosition:           -1,
		ReferenceAllele:       -1,
		TumorSeqAllele2:       -1,
		VariantType:           -1,
		NCBIBuild:             -1,
	}
	w := NewMAFWriter(&buf, header, cols)
	w.SetReplace(true)
	w.SetAllEffects(false)
	w.SetSources([]annotate.AnnotationSource{&testSource{
		name:    "oncokb",
		version: "v1",
		columns: []annotate.ColumnDef{{Name: "gene_type", Description: "Gene type"}},
	}})

	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteRow(fields, ann, []*annotate.Annotation{ann}, &vcf.Variant{Ref: "C", Alt: "T"}); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")

	// Header: original + source columns (no vibe. prefix in replace mode)
	if lines[0] != "Hugo_Symbol\tConsequence\toncokb.gene_type" {
		t.Errorf("header = %q, want Hugo_Symbol\\tConsequence\\toncokb.gene_type", lines[0])
	}

	// Row: overwritten core + source columns
	parts := strings.Split(lines[1], "\t")
	if parts[0] != "KRAS" {
		t.Errorf("Hugo_Symbol = %q, want KRAS (overwritten)", parts[0])
	}
	if parts[1] != "missense_variant" {
		t.Errorf("Consequence = %q, want missense_variant (overwritten)", parts[1])
	}
	if parts[2] != "ONCOGENE" {
		t.Errorf("oncokb.gene_type = %q, want ONCOGENE", parts[2])
	}
}

func TestMAFWriter_ReplaceNilAnnotation(t *testing.T) {
	fields := []string{"GENE", "conseq"}

	var buf bytes.Buffer
	cols := maf.ColumnIndices{
		HugoSymbol:            0,
		Consequence:           1,
		HGVSpShort:            -1,
		TranscriptID:          -1,
		HGVSc:                 -1,
		VariantClassification: -1,
		HGVSp:                 -1,
		Chromosome:            -1,
		StartPosition:         -1,
		EndPosition:           -1,
		ReferenceAllele:       -1,
		TumorSeqAllele2:       -1,
		VariantType:           -1,
		NCBIBuild:             -1,
	}
	w := NewMAFWriter(&buf, "Hugo_Symbol\tConsequence", cols)
	w.SetReplace(true)
	w.SetAllEffects(false)

	if err := w.WriteRow(fields, nil, nil, nil); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	parts := strings.Split(strings.TrimRight(buf.String(), "\n"), "\t")
	// Original values should be preserved when annotation is nil
	if parts[0] != "GENE" {
		t.Errorf("Hugo_Symbol = %q, want GENE (unchanged with nil annotation)", parts[0])
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
func (s *testSource) MatchLevel() annotate.MatchLevel                      { return annotate.MatchGene }
func (s *testSource) Columns() []annotate.ColumnDef                        { return s.columns }
func (s *testSource) Annotate(v *vcf.Variant, anns []*annotate.Annotation) {}
