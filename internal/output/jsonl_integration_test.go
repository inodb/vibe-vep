package output_test

import (
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/input"
	"github.com/inodb/vibe-vep/internal/output"
)

// TestStreamAnnotation verifies that the stream annotation pipeline produces
// valid VEP-compatible JSONL output for well-known cancer variants.
func TestStreamAnnotation(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping stream integration test in short mode")
	}

	gtfPath, fastaPath, canonicalPath := findGENCODEFiles(t, "GRCh38")
	c := cache.New()
	loader := cache.NewGENCODELoader(gtfPath, fastaPath)
	if canonicalPath != "" {
		mskOverrides, ensOverrides, _ := cache.LoadBiomartCanonicals(canonicalPath)
		loader.SetCanonicalOverrides(mskOverrides, ensOverrides)
	}
	if err := loader.Load(c); err != nil {
		t.Fatalf("load GENCODE: %v", err)
	}
	c.BuildIndex()

	ann := annotate.NewAnnotator(c)

	tests := []struct {
		name       string
		input      string
		wantGene   string
		wantConseq string
		wantAA     string // amino_acids in VEP format (e.g. "V/E")
	}{
		{
			name:       "BRAF V600E",
			input:      `{"chromosome":"7","start":140753336,"end":140753336,"referenceAllele":"A","variantAllele":"T"}`,
			wantGene:   "BRAF",
			wantConseq: "missense_variant",
			wantAA:     "V/E",
		},
		{
			name:       "KRAS G12C",
			input:      `{"chromosome":"12","start":25245351,"end":25245351,"referenceAllele":"C","variantAllele":"A"}`,
			wantGene:   "KRAS",
			wantConseq: "missense_variant",
			wantAA:     "G/C",
		},
		{
			name:       "TP53 R175H",
			input:      `{"chromosome":"17","start":7675088,"end":7675088,"referenceAllele":"C","variantAllele":"T"}`,
			wantGene:   "TP53",
			wantConseq: "missense_variant",
		},
		{
			name:       "EGFR L858R",
			input:      `{"chromosome":"7","start":55191822,"end":55191822,"referenceAllele":"T","variantAllele":"G"}`,
			wantGene:   "EGFR",
			wantConseq: "missense_variant",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gl, err := input.ParseGenomicLocation([]byte(tt.input))
			if err != nil {
				t.Fatalf("parse: %v", err)
			}

			v := gl.ToVariant()
			anns, err := ann.Annotate(v)
			if err != nil {
				t.Fatalf("annotate: %v", err)
			}

			// Write as VEP JSONL.
			var buf strings.Builder
			w := output.NewJSONLWriter(&buf, "ensembl-vep-jsonl", "GRCh38")
			w.SetInput(gl.FormatInput())
			for _, a := range anns {
				w.Write(v, a)
			}
			w.Flush()

			// Parse the output.
			var result output.VEPVariantAnnotation
			if err := json.Unmarshal([]byte(buf.String()), &result); err != nil {
				t.Fatalf("invalid JSON output: %v\nraw: %s", err, buf.String())
			}

			// Verify top-level fields.
			if result.AssemblyName != "GRCh38" {
				t.Errorf("assembly_name=%q, want GRCh38", result.AssemblyName)
			}
			if result.MostSevereConsequence == "" {
				t.Error("most_severe_consequence is empty")
			}

			// Find the best transcript consequence for the expected gene.
			// Prefer transcripts with the expected consequence over others.
			var bestTC *output.VEPTranscriptConsequence
			var fallbackTC *output.VEPTranscriptConsequence
			for i, tc := range result.TranscriptConsequences {
				if tc.GeneSymbol != tt.wantGene {
					continue
				}
				if fallbackTC == nil {
					fallbackTC = &result.TranscriptConsequences[i]
				}
				for _, ct := range tc.ConsequenceTerms {
					if ct == tt.wantConseq {
						bestTC = &result.TranscriptConsequences[i]
						break
					}
				}
				if bestTC != nil {
					break
				}
			}
			if bestTC == nil {
				bestTC = fallbackTC
			}
			if bestTC == nil {
				genes := make([]string, 0)
				for _, tc := range result.TranscriptConsequences {
					genes = append(genes, tc.GeneSymbol)
				}
				t.Fatalf("gene %q not found in transcript_consequences (found: %v)", tt.wantGene, genes)
			}

			hasConseq := false
			for _, ct := range bestTC.ConsequenceTerms {
				if ct == tt.wantConseq {
					hasConseq = true
					break
				}
			}
			if !hasConseq {
				t.Errorf("gene %s: consequence_terms=%v, want to contain %q",
					tt.wantGene, bestTC.ConsequenceTerms, tt.wantConseq)
			}
			if tt.wantAA != "" && bestTC.AminoAcids != tt.wantAA {
				t.Errorf("gene %s: amino_acids=%q, want %q", tt.wantGene, bestTC.AminoAcids, tt.wantAA)
			}
			t.Logf("gene=%s transcript=%s consequence=%v amino_acids=%s protein_start=%d hgvsp=%s",
				bestTC.GeneSymbol, bestTC.TranscriptID, bestTC.ConsequenceTerms,
				bestTC.AminoAcids, bestTC.ProteinStart, bestTC.HGVSp)
		})
	}
}

// findGENCODEFiles and findStudyDir are shared test helpers defined in
// validation_benchmark_test.go. This file uses them for integration testing.
func init() {
	// Ensure testdata paths work from the output package test directory.
	_ = filepath.Join("testdata")
	_ = os.Getenv("HOME")
}
