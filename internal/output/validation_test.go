package output

import (
	"bytes"
	"strings"
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/vcf"
)

func TestValidationWriter_Match(t *testing.T) {
	var buf bytes.Buffer
	w := NewValidationWriter(&buf, true)

	if err := w.WriteHeader(); err != nil {
		t.Fatalf("WriteHeader failed: %v", err)
	}

	variant := &vcf.Variant{
		Chrom: "12",
		Pos:   25398284,
		Ref:   "C",
		Alt:   "A",
	}

	mafAnn := &maf.MAFAnnotation{
		HugoSymbol:   "KRAS",
		Consequence:  "missense_variant",
		HGVSpShort:   "p.G12C",
		TranscriptID: "ENST00000256078",
	}

	vepAnns := []*annotate.Annotation{
		{
			TranscriptID:    "ENST00000256078",
			GeneName:        "KRAS",
			Consequence:     "missense_variant",
			AminoAcidChange: "G/C",
			IsCanonical:     true,
		},
	}

	if err := w.WriteComparison(variant, mafAnn, vepAnns); err != nil {
		t.Fatalf("WriteComparison failed: %v", err)
	}

	if err := w.Flush(); err != nil {
		t.Fatalf("Flush failed: %v", err)
	}

	output := buf.String()
	if !strings.Contains(output, "KRAS") {
		t.Errorf("Expected output to contain KRAS, got: %s", output)
	}
	if !strings.Contains(output, "Y") {
		t.Errorf("Expected match (Y) in output, got: %s", output)
	}

	total, matches, mismatches := w.Summary()
	if total != 1 || matches != 1 || mismatches != 0 {
		t.Errorf("Expected 1 total, 1 match, 0 mismatches; got %d, %d, %d", total, matches, mismatches)
	}
}

func TestValidationWriter_Mismatch(t *testing.T) {
	var buf bytes.Buffer
	w := NewValidationWriter(&buf, true)

	w.WriteHeader()

	variant := &vcf.Variant{
		Chrom: "12",
		Pos:   25398284,
		Ref:   "C",
		Alt:   "A",
	}

	mafAnn := &maf.MAFAnnotation{
		HugoSymbol:  "KRAS",
		Consequence: "missense_variant",
		HGVSpShort:  "p.G12C",
	}

	vepAnns := []*annotate.Annotation{
		{
			Consequence: "synonymous_variant", // Wrong prediction
		},
	}

	w.WriteComparison(variant, mafAnn, vepAnns)
	w.Flush()

	output := buf.String()
	if !strings.Contains(output, "N") {
		t.Errorf("Expected mismatch (N) in output, got: %s", output)
	}

	_, matches, mismatches := w.Summary()
	if matches != 0 || mismatches != 1 {
		t.Errorf("Expected 0 matches, 1 mismatch; got %d, %d", matches, mismatches)
	}
}

func TestNormalizeConsequence(t *testing.T) {
	tests := []struct {
		input    string
		expected string
	}{
		{"missense_variant", "missense_variant"},
		{"Missense_Mutation", "missense_variant"},
		{"MISSENSE_MUTATION", "missense_variant"},
		{"Nonsense_Mutation", "stop_gained"},
		{"stop_gained", "stop_gained"},
		{"Silent", "synonymous_variant"},
		{"Frame_Shift_Del", "frameshift_variant"},
		{"Frame_Shift_Ins", "frameshift_variant"},
		{"In_Frame_Del", "inframe_deletion"},
		{"unknown_type", "unknown_type"},
	}

	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			result := normalizeConsequence(tt.input)
			if result != tt.expected {
				t.Errorf("normalizeConsequence(%q) = %q, want %q", tt.input, result, tt.expected)
			}
		})
	}
}

func TestValidationWriter_MismatchesOnly(t *testing.T) {
	var buf bytes.Buffer
	w := NewValidationWriter(&buf, false) // showAll = false

	w.WriteHeader()

	variant := &vcf.Variant{Chrom: "12", Pos: 100, Ref: "A", Alt: "T"}

	// Write a match - should not appear in output
	mafAnn := &maf.MAFAnnotation{Consequence: "missense_variant"}
	vepAnns := []*annotate.Annotation{{Consequence: "missense_variant"}}
	w.WriteComparison(variant, mafAnn, vepAnns)

	// Write a mismatch - should appear
	mafAnn2 := &maf.MAFAnnotation{Consequence: "stop_gained", HugoSymbol: "TP53"}
	vepAnns2 := []*annotate.Annotation{{Consequence: "missense_variant"}}
	w.WriteComparison(variant, mafAnn2, vepAnns2)

	w.Flush()

	output := buf.String()
	lines := strings.Split(strings.TrimSpace(output), "\n")

	// Should have header + 1 mismatch line only
	if len(lines) != 2 {
		t.Errorf("Expected 2 lines (header + mismatch), got %d: %v", len(lines), lines)
	}

	if !strings.Contains(output, "TP53") {
		t.Errorf("Expected mismatch line with TP53, got: %s", output)
	}
}
