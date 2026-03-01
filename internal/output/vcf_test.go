package output

import (
	"bytes"
	"strings"
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

func TestVCFWriter_Header(t *testing.T) {
	headers := []string{
		"##fileformat=VCFv4.2",
		"##reference=GRCh38",
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
	}

	var buf bytes.Buffer
	w := NewVCFWriter(&buf, headers)
	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	out := buf.String()
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")

	// Original ## lines preserved
	if lines[0] != "##fileformat=VCFv4.2" {
		t.Errorf("first line = %q, want ##fileformat=VCFv4.2", lines[0])
	}
	if lines[1] != "##reference=GRCh38" {
		t.Errorf("second line = %q, want ##reference=GRCh38", lines[1])
	}

	// CSQ INFO line inserted before #CHROM
	csqIdx := -1
	chromIdx := -1
	for i, line := range lines {
		if strings.HasPrefix(line, "##INFO=<ID=CSQ") {
			csqIdx = i
		}
		if strings.HasPrefix(line, "#CHROM") {
			chromIdx = i
		}
	}
	if csqIdx < 0 {
		t.Fatal("CSQ INFO line not found in header")
	}
	if chromIdx < 0 {
		t.Fatal("#CHROM line not found in header")
	}
	if csqIdx >= chromIdx {
		t.Errorf("CSQ INFO line (idx %d) should appear before #CHROM line (idx %d)", csqIdx, chromIdx)
	}

	// CSQ line contains expected format
	if !strings.Contains(lines[csqIdx], "Allele|Consequence|IMPACT|SYMBOL") {
		t.Errorf("CSQ line missing expected fields: %s", lines[csqIdx])
	}
}

func TestVCFWriter_SingleAnnotation(t *testing.T) {
	headers := []string{
		"##fileformat=VCFv4.2",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
	}

	var buf bytes.Buffer
	w := NewVCFWriter(&buf, headers)
	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}

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
		Allele:          "A",
		Consequence:     "missense_variant",
		Impact:          "MODERATE",
		GeneName:        "KRAS",
		GeneID:          "ENSG00000133703",
		TranscriptID:    "ENST00000311936",
		Biotype:         "protein_coding",
		ExonNumber:      "2/5",
		HGVSc:           "c.34G>T",
		HGVSp:           "p.Gly12Cys",
		CDSPosition:     34,
		ProteinPosition: 12,
		AminoAcidChange: "G12C",
		CodonChange:     "ggt/Tgt",
		IsCanonical:     true,
	}

	if err := w.Write(v, ann); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	out := buf.String()
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")

	// Last line should be the data line
	dataLine := lines[len(lines)-1]

	// Should start with variant fields
	if !strings.HasPrefix(dataLine, "12\t25245351\trs121913529\tC\tA\t") {
		t.Errorf("data line has wrong prefix: %s", dataLine)
	}

	// Should contain CSQ
	if !strings.Contains(dataLine, "CSQ=") {
		t.Fatalf("data line missing CSQ: %s", dataLine)
	}

	// Parse CSQ value
	fields := strings.Split(dataLine, "\t")
	info := fields[7]
	if !strings.HasPrefix(info, "CSQ=") {
		t.Fatalf("INFO should be CSQ=..., got: %s", info)
	}

	csq := strings.TrimPrefix(info, "CSQ=")
	parts := strings.Split(csq, "|")
	if len(parts) != len(csqFields) {
		t.Fatalf("CSQ has %d fields, want %d", len(parts), len(csqFields))
	}

	// Check specific CSQ fields
	checks := map[int]string{
		0:  "A",                // Allele
		1:  "missense_variant", // Consequence
		2:  "MODERATE",         // IMPACT
		3:  "KRAS",             // SYMBOL
		6:  "ENST00000311936",  // Feature
		7:  "protein_coding",   // BIOTYPE
		10: "c.34G>T",          // HGVSc
		11: "p.Gly12Cys",       // HGVSp
		17: "YES",              // CANONICAL
	}
	for idx, want := range checks {
		if parts[idx] != want {
			t.Errorf("CSQ field %d (%s) = %q, want %q", idx, csqFields[idx], parts[idx], want)
		}
	}
}

func TestVCFWriter_MultipleAnnotations(t *testing.T) {
	headers := []string{
		"##fileformat=VCFv4.2",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
	}

	var buf bytes.Buffer
	w := NewVCFWriter(&buf, headers)
	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}

	v := &vcf.Variant{
		Chrom:  "12",
		Pos:    25245351,
		ID:     ".",
		Ref:    "C",
		Alt:    "A",
		Qual:   0,
		Filter: "PASS",
		Info:   map[string]interface{}{},
	}

	ann1 := &annotate.Annotation{
		Allele:       "A",
		Consequence:  "missense_variant",
		Impact:       "MODERATE",
		GeneName:     "KRAS",
		TranscriptID: "ENST00000311936",
		Biotype:      "protein_coding",
		IsCanonical:  true,
	}
	ann2 := &annotate.Annotation{
		Allele:       "A",
		Consequence:  "missense_variant",
		Impact:       "MODERATE",
		GeneName:     "KRAS",
		TranscriptID: "ENST00000256078",
		Biotype:      "protein_coding",
	}

	if err := w.Write(v, ann1); err != nil {
		t.Fatal(err)
	}
	if err := w.Write(v, ann2); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	out := buf.String()
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	dataLine := lines[len(lines)-1]

	// Should have exactly one data line
	dataLines := 0
	for _, l := range lines {
		if !strings.HasPrefix(l, "#") {
			dataLines++
		}
	}
	if dataLines != 1 {
		t.Errorf("expected 1 data line, got %d", dataLines)
	}

	// CSQ should have two comma-separated entries
	fields := strings.Split(dataLine, "\t")
	info := fields[7]
	csqValue := strings.TrimPrefix(info, "CSQ=")
	entries := strings.Split(csqValue, ",")
	if len(entries) != 2 {
		t.Fatalf("expected 2 CSQ entries, got %d: %s", len(entries), csqValue)
	}

	// First entry should have ENST00000311936
	if !strings.Contains(entries[0], "ENST00000311936") {
		t.Errorf("first CSQ entry missing ENST00000311936: %s", entries[0])
	}
	// Second entry should have ENST00000256078
	if !strings.Contains(entries[1], "ENST00000256078") {
		t.Errorf("second CSQ entry missing ENST00000256078: %s", entries[1])
	}

	// QUAL should be "." when 0
	if fields[5] != "." {
		t.Errorf("QUAL = %q, want \".\"", fields[5])
	}
}

func TestVCFWriter_PreservesOriginalInfo(t *testing.T) {
	headers := []string{
		"##fileformat=VCFv4.2",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
	}

	var buf bytes.Buffer
	w := NewVCFWriter(&buf, headers)
	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}

	v := &vcf.Variant{
		Chrom:  "12",
		Pos:    25245351,
		ID:     ".",
		Ref:    "C",
		Alt:    "A",
		Qual:   100,
		Filter: "PASS",
		Info:   map[string]interface{}{"DP": "50"},
	}

	ann := &annotate.Annotation{
		Allele:       "A",
		Consequence:  "missense_variant",
		Impact:       "MODERATE",
		GeneName:     "KRAS",
		TranscriptID: "ENST00000311936",
		Biotype:      "protein_coding",
	}

	if err := w.Write(v, ann); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	out := buf.String()
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	dataLine := lines[len(lines)-1]
	fields := strings.Split(dataLine, "\t")
	info := fields[7]

	// Should contain both DP and CSQ
	if !strings.Contains(info, "DP=50") {
		t.Errorf("INFO missing DP=50: %s", info)
	}
	if !strings.Contains(info, "CSQ=") {
		t.Errorf("INFO missing CSQ: %s", info)
	}

	// DP should come before CSQ
	dpIdx := strings.Index(info, "DP=50")
	csqIdx := strings.Index(info, "CSQ=")
	if dpIdx > csqIdx {
		t.Errorf("DP should appear before CSQ in INFO: %s", info)
	}
}

func TestVCFWriter_MultipleVariants(t *testing.T) {
	headers := []string{
		"##fileformat=VCFv4.2",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
	}

	var buf bytes.Buffer
	w := NewVCFWriter(&buf, headers)
	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}

	v1 := &vcf.Variant{
		Chrom: "12", Pos: 25245351, ID: ".", Ref: "C", Alt: "A",
		Filter: "PASS", Info: map[string]interface{}{},
	}
	v2 := &vcf.Variant{
		Chrom: "17", Pos: 7577121, ID: ".", Ref: "G", Alt: "A",
		Filter: "PASS", Info: map[string]interface{}{},
	}

	ann1 := &annotate.Annotation{
		Allele: "A", Consequence: "missense_variant", Impact: "MODERATE",
		GeneName: "KRAS", TranscriptID: "ENST00000311936", Biotype: "protein_coding",
	}
	ann2 := &annotate.Annotation{
		Allele: "A", Consequence: "missense_variant", Impact: "MODERATE",
		GeneName: "TP53", TranscriptID: "ENST00000269305", Biotype: "protein_coding",
	}

	// Write variant 1
	if err := w.Write(v1, ann1); err != nil {
		t.Fatal(err)
	}
	// Write variant 2 (should flush variant 1)
	if err := w.Write(v2, ann2); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	out := buf.String()
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")

	// Count data lines
	var dataLines []string
	for _, l := range lines {
		if !strings.HasPrefix(l, "#") {
			dataLines = append(dataLines, l)
		}
	}

	if len(dataLines) != 2 {
		t.Fatalf("expected 2 data lines, got %d", len(dataLines))
	}

	// First line should be KRAS variant
	if !strings.HasPrefix(dataLines[0], "12\t25245351\t") {
		t.Errorf("first data line wrong: %s", dataLines[0])
	}
	if !strings.Contains(dataLines[0], "KRAS") {
		t.Errorf("first data line should contain KRAS: %s", dataLines[0])
	}

	// Second line should be TP53 variant
	if !strings.HasPrefix(dataLines[1], "17\t7577121\t") {
		t.Errorf("second data line wrong: %s", dataLines[1])
	}
	if !strings.Contains(dataLines[1], "TP53") {
		t.Errorf("second data line should contain TP53: %s", dataLines[1])
	}
}
