package output

import (
	"bytes"
	"encoding/json"
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

func TestJSONLWriterVEP(t *testing.T) {
	var buf bytes.Buffer
	w := NewJSONLWriter(&buf, "ensembl-vep-jsonl", "GRCh38")
	w.SetInput("7,140453136,140453136,A,T")

	v := &vcf.Variant{Chrom: "7", Pos: 140453136, Ref: "A", Alt: "T"}
	ann := &annotate.Annotation{
		TranscriptID:    "ENST00000288602.11",
		GeneName:        "BRAF",
		GeneID:          "ENSG00000157764",
		Consequence:     "missense_variant",
		Impact:          "MODERATE",
		Allele:          "T",
		Biotype:         "protein_coding",
		ProteinPosition: 600,
		CDSPosition:     1799,
		CDNAPosition:    1906,
		AminoAcidChange: "V600E",
		CodonChange:     "gTg/gAg",
		ExonNumber:      "15/18",
		HGVSc:           "c.1799T>A",
		HGVSp:           "p.Val600Glu",
	}

	if err := w.Write(v, ann); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	line := buf.String()
	if line == "" {
		t.Fatal("empty output")
	}

	// Parse back as VEP JSON.
	var result VEPVariantAnnotation
	if err := json.Unmarshal([]byte(line), &result); err != nil {
		t.Fatalf("invalid JSON: %v\nline: %s", err, line)
	}

	// Top level.
	if result.SeqRegionName != "7" {
		t.Errorf("seq_region_name=%q, want %q", result.SeqRegionName, "7")
	}
	if result.Start != 140453136 {
		t.Errorf("start=%d, want %d", result.Start, 140453136)
	}
	if result.AssemblyName != "GRCh38" {
		t.Errorf("assembly_name=%q, want %q", result.AssemblyName, "GRCh38")
	}
	if result.MostSevereConsequence != "missense_variant" {
		t.Errorf("most_severe_consequence=%q, want %q", result.MostSevereConsequence, "missense_variant")
	}
	if result.AlleleString != "A/T" {
		t.Errorf("allele_string=%q, want %q", result.AlleleString, "A/T")
	}
	if result.Input != "7,140453136,140453136,A,T" {
		t.Errorf("input=%q, want %q", result.Input, "7,140453136,140453136,A,T")
	}

	// Transcript consequence.
	if len(result.TranscriptConsequences) != 1 {
		t.Fatalf("got %d transcript_consequences, want 1", len(result.TranscriptConsequences))
	}
	tc := result.TranscriptConsequences[0]
	if tc.TranscriptID != "ENST00000288602" {
		t.Errorf("transcript_id=%q, want %q (version stripped)", tc.TranscriptID, "ENST00000288602")
	}
	if tc.GeneSymbol != "BRAF" {
		t.Errorf("gene_symbol=%q, want %q", tc.GeneSymbol, "BRAF")
	}
	if len(tc.ConsequenceTerms) != 1 || tc.ConsequenceTerms[0] != "missense_variant" {
		t.Errorf("consequence_terms=%v, want [missense_variant]", tc.ConsequenceTerms)
	}
	if tc.AminoAcids != "V/E" {
		t.Errorf("amino_acids=%q, want %q", tc.AminoAcids, "V/E")
	}
	if tc.ProteinStart != 600 {
		t.Errorf("protein_start=%d, want %d", tc.ProteinStart, 600)
	}
	if tc.CDSStart != 1799 {
		t.Errorf("cds_start=%d, want %d", tc.CDSStart, 1799)
	}
	if tc.HGVSc != "c.1799T>A" {
		t.Errorf("hgvsc=%q, want %q", tc.HGVSc, "c.1799T>A")
	}
}

func TestJSONLWriterVibeVep(t *testing.T) {
	var buf bytes.Buffer
	w := NewJSONLWriter(&buf, "vibe-vep-jsonl", "GRCh38")
	w.SetInput("12,25245350,25245350,C,A")

	v := &vcf.Variant{Chrom: "12", Pos: 25245350, Ref: "C", Alt: "A"}
	ann := &annotate.Annotation{
		TranscriptID:       "ENST00000311936.8",
		GeneName:           "KRAS",
		GeneID:             "ENSG00000133703",
		Consequence:        "missense_variant",
		Impact:             "MODERATE",
		Allele:             "A",
		Biotype:            "protein_coding",
		ProteinPosition:    12,
		AminoAcidChange:    "G12C",
		CodonChange:        "Ggt/Tgt",
		HGVSp:              "p.Gly12Cys",
		IsCanonicalMSK:     true,
		IsCanonicalEnsembl: true,
	}

	if err := w.Write(v, ann); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	var result VibeVepVariantAnnotation
	if err := json.Unmarshal(buf.Bytes(), &result); err != nil {
		t.Fatalf("invalid JSON: %v", err)
	}

	if result.Chromosome != "12" {
		t.Errorf("chromosome=%q, want %q", result.Chromosome, "12")
	}
	if len(result.TranscriptConsequences) != 1 {
		t.Fatalf("got %d consequences, want 1", len(result.TranscriptConsequences))
	}
	tc := result.TranscriptConsequences[0]
	if tc.TranscriptID != "ENST00000311936.8" {
		t.Errorf("transcript_id=%q, want with version", tc.TranscriptID)
	}
	if tc.HugoSymbol != "KRAS" {
		t.Errorf("hugo_symbol=%q, want %q", tc.HugoSymbol, "KRAS")
	}
	if tc.AminoAcidChange != "G12C" {
		t.Errorf("amino_acid_change=%q, want %q", tc.AminoAcidChange, "G12C")
	}
	if !tc.CanonicalMSKCC {
		t.Error("canonical_mskcc should be true")
	}
}

func TestJSONLWriterError(t *testing.T) {
	var buf bytes.Buffer
	w := NewJSONLWriter(&buf, "ensembl-vep-jsonl", "GRCh38")
	w.WriteError("bad input", "parse error: missing chromosome")

	var result map[string]string
	if err := json.Unmarshal(buf.Bytes(), &result); err != nil {
		t.Fatalf("invalid JSON: %v", err)
	}
	if result["error"] != "parse error: missing chromosome" {
		t.Errorf("error=%q", result["error"])
	}
	if result["input"] != "bad input" {
		t.Errorf("input=%q", result["input"])
	}
}

func TestJSONLWriterWarnings(t *testing.T) {
	var buf bytes.Buffer
	w := NewJSONLWriter(&buf, "ensembl-vep-jsonl", "GRCh38")
	w.SetInput("test")
	w.AddWarning("requested transcript version ENST00000311936.99 not found, using ENST00000311936.8 instead")

	v := &vcf.Variant{Chrom: "12", Pos: 25245351, Ref: "C", Alt: "A"}
	ann := &annotate.Annotation{
		TranscriptID: "ENST00000311936.8",
		GeneName:     "KRAS",
		Consequence:  "missense_variant",
		Impact:       "MODERATE",
		Allele:       "A",
	}
	w.Write(v, ann)
	w.Flush()

	var result VEPVariantAnnotation
	json.Unmarshal(buf.Bytes(), &result)

	if len(result.Warnings) != 1 {
		t.Fatalf("got %d warnings, want 1", len(result.Warnings))
	}
	if result.Warnings[0] != "requested transcript version ENST00000311936.99 not found, using ENST00000311936.8 instead" {
		t.Errorf("warning=%q", result.Warnings[0])
	}
}

func TestJSONLWriterNoWarnings(t *testing.T) {
	var buf bytes.Buffer
	w := NewJSONLWriter(&buf, "ensembl-vep-jsonl", "GRCh38")

	v := &vcf.Variant{Chrom: "12", Pos: 25245351, Ref: "C", Alt: "A"}
	ann := &annotate.Annotation{
		TranscriptID: "ENST00000311936.8",
		Consequence:  "missense_variant",
		Impact:       "MODERATE",
		Allele:       "A",
	}
	w.Write(v, ann)
	w.Flush()

	var result VEPVariantAnnotation
	json.Unmarshal(buf.Bytes(), &result)

	if len(result.Warnings) != 0 {
		t.Errorf("got %d warnings, want 0", len(result.Warnings))
	}
	// Verify warnings field is omitted from JSON (not null).
	if bytes.Contains(buf.Bytes(), []byte(`"warnings"`)) {
		t.Error("warnings field should be omitted when empty")
	}
}

func TestJSONLWriterSIFTPolyPhen(t *testing.T) {
	var buf bytes.Buffer
	w := NewJSONLWriter(&buf, "ensembl-vep-jsonl", "GRCh38")

	v := &vcf.Variant{Chrom: "7", Pos: 140453136, Ref: "A", Alt: "T"}
	ann := &annotate.Annotation{
		TranscriptID:    "ENST00000288602.11",
		GeneName:        "BRAF",
		Consequence:     "missense_variant",
		Impact:          "MODERATE",
		Allele:          "T",
		ProteinPosition: 600,
	}
	ann.SetExtraKey("sift.score", "0.000")
	ann.SetExtraKey("sift.prediction", "deleterious")
	ann.SetExtraKey("polyphen.score", "0.935")
	ann.SetExtraKey("polyphen.prediction", "probably_damaging")

	w.Write(v, ann)
	w.Flush()

	var result VEPVariantAnnotation
	json.Unmarshal(buf.Bytes(), &result)

	tc := result.TranscriptConsequences[0]
	if tc.SIFTScore == nil || *tc.SIFTScore != 0.0 {
		t.Errorf("sift_score=%v, want 0.0", tc.SIFTScore)
	}
	if tc.SIFTPrediction != "deleterious" {
		t.Errorf("sift_prediction=%q, want %q", tc.SIFTPrediction, "deleterious")
	}
	if tc.PolyPhenScore == nil || *tc.PolyPhenScore != 0.935 {
		t.Errorf("polyphen_score=%v, want 0.935", tc.PolyPhenScore)
	}
	if tc.PolyPhenPrediction != "probably_damaging" {
		t.Errorf("polyphen_prediction=%q, want %q", tc.PolyPhenPrediction, "probably_damaging")
	}
}
