package server

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"

	"github.com/inodb/vibe-vep/internal/input"
	"github.com/inodb/vibe-vep/internal/output"
)

// --- GN Genomic GET ---

func TestGNGenomicGet_KRASG12C(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	// KRAS G12C: chr12:25245351 C>A
	req := httptest.NewRequest(http.MethodGet, "/genome-nexus/grch38/annotation/genomic/12,25245351,25245351,C,A", nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var resp output.GNAnnotation
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v\nbody: %s", err, w.Body.String())
	}
	if resp.AssemblyName != "GRCh38" {
		t.Errorf("assembly: got %q, want GRCh38", resp.AssemblyName)
	}
	if !resp.SuccessfullyAnnotated {
		t.Error("expected successfully_annotated=true")
	}
	if resp.MostSevereConsequence != "missense_variant" {
		t.Errorf("most_severe_consequence: got %q, want missense_variant", resp.MostSevereConsequence)
	}
	if resp.SeqRegionName != "12" {
		t.Errorf("seq_region_name: got %q, want 12", resp.SeqRegionName)
	}
	if resp.Start != 25245351 {
		t.Errorf("start: got %d, want 25245351", resp.Start)
	}
	if resp.AlleleString != "C/A" {
		t.Errorf("allele_string: got %q, want C/A", resp.AlleleString)
	}

	// Check transcript consequences.
	if len(resp.TranscriptConsequences) == 0 {
		t.Fatal("expected transcript consequences, got none")
	}

	// Find canonical transcript.
	var canonical *output.GNTranscriptConsequence
	for i := range resp.TranscriptConsequences {
		if resp.TranscriptConsequences[i].TranscriptID == "ENST00000311936" {
			canonical = &resp.TranscriptConsequences[i]
			break
		}
	}
	if canonical == nil {
		t.Fatal("canonical transcript ENST00000311936 not found")
	}
	if canonical.GeneSymbol != "KRAS" {
		t.Errorf("gene_symbol: got %q, want KRAS", canonical.GeneSymbol)
	}
	if canonical.Impact != "MODERATE" {
		t.Errorf("impact: got %q, want MODERATE", canonical.Impact)
	}
	if !contains(canonical.ConsequenceTerms, "missense_variant") {
		t.Errorf("consequence_terms: got %v, want to contain missense_variant", canonical.ConsequenceTerms)
	}
	if canonical.AminoAcids != "G/C" {
		t.Errorf("amino_acids: got %q, want G/C", canonical.AminoAcids)
	}
	if canonical.ProteinStart != 12 {
		t.Errorf("protein_start: got %d, want 12", canonical.ProteinStart)
	}
	if canonical.HGVSp != "p.Gly12Cys" {
		t.Errorf("hgvsp: got %q, want p.Gly12Cys", canonical.HGVSp)
	}
	if canonical.Canonical != "1" {
		t.Errorf("canonical: got %q, want 1", canonical.Canonical)
	}
	if canonical.Biotype != "protein_coding" {
		t.Errorf("biotype: got %q, want protein_coding", canonical.Biotype)
	}
}

func TestGNGenomicGet_Intergenic(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	req := httptest.NewRequest(http.MethodGet, "/genome-nexus/grch38/annotation/genomic/12,1000000,1000000,A,T", nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var resp output.GNAnnotation
	json.Unmarshal(w.Body.Bytes(), &resp)
	if resp.MostSevereConsequence != "intergenic_variant" {
		t.Errorf("expected intergenic_variant, got %q", resp.MostSevereConsequence)
	}
	if !resp.SuccessfullyAnnotated {
		t.Error("expected successfully_annotated=true even for intergenic")
	}
}

// --- GN Genomic POST ---

func TestGNGenomicPost_Batch(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	body := `[
		{"chromosome":"12","start":25245351,"end":25245351,"referenceAllele":"C","variantAllele":"A"},
		{"chromosome":"12","start":25245350,"end":25245350,"referenceAllele":"C","variantAllele":"T"},
		{"chromosome":"12","start":1000000,"end":1000000,"referenceAllele":"A","variantAllele":"G"}
	]`
	req := httptest.NewRequest(http.MethodPost, "/genome-nexus/grch38/annotation/genomic", strings.NewReader(body))
	req.Header.Set("Content-Type", "application/json")
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var resp []output.GNAnnotation
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v", err)
	}
	if len(resp) != 3 {
		t.Fatalf("expected 3 annotations, got %d", len(resp))
	}

	// First two are in KRAS, third is intergenic.
	if resp[0].MostSevereConsequence == "intergenic_variant" {
		t.Error("variant 0 should not be intergenic (KRAS position)")
	}
	if resp[1].MostSevereConsequence == "intergenic_variant" {
		t.Error("variant 1 should not be intergenic (KRAS position)")
	}
	if resp[2].MostSevereConsequence != "intergenic_variant" {
		t.Errorf("variant 2 should be intergenic, got %q", resp[2].MostSevereConsequence)
	}

	// All should be successfully annotated.
	for i, a := range resp {
		if !a.SuccessfullyAnnotated {
			t.Errorf("variant %d: expected successfully_annotated=true", i)
		}
	}
}

// --- GN HGVS GET ---

func TestGNHGVSGet_KRASG12C(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	// Use genomic coordinate format (parsed as SpecGenomic by ParseVariantSpec).
	req := httptest.NewRequest(http.MethodGet, "/genome-nexus/grch38/annotation/12:25245351:C:A", nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var resp output.GNAnnotation
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v\nbody: %s", err, w.Body.String())
	}
	if resp.MostSevereConsequence != "missense_variant" {
		t.Errorf("most_severe_consequence: got %q, want missense_variant", resp.MostSevereConsequence)
	}
	if !resp.SuccessfullyAnnotated {
		t.Error("expected successfully_annotated=true")
	}

	// Check canonical transcript.
	var canonical *output.GNTranscriptConsequence
	for i := range resp.TranscriptConsequences {
		if resp.TranscriptConsequences[i].TranscriptID == "ENST00000311936" {
			canonical = &resp.TranscriptConsequences[i]
			break
		}
	}
	if canonical == nil {
		t.Fatal("ENST00000311936 not found")
	}
	if canonical.GeneSymbol != "KRAS" {
		t.Errorf("gene_symbol: got %q, want KRAS", canonical.GeneSymbol)
	}
	if canonical.AminoAcids != "G/C" {
		t.Errorf("amino_acids: got %q, want G/C", canonical.AminoAcids)
	}
	if canonical.ProteinStart != 12 {
		t.Errorf("protein_start: got %d, want 12", canonical.ProteinStart)
	}
}

func TestGNHGVSPost_Batch(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	body := `["12:25245351:C:A", "12:25245350:C:T"]`
	req := httptest.NewRequest(http.MethodPost, "/genome-nexus/grch38/annotation", strings.NewReader(body))
	req.Header.Set("Content-Type", "application/json")
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var resp []output.GNAnnotation
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v", err)
	}
	if len(resp) != 2 {
		t.Fatalf("expected 2 annotations, got %d", len(resp))
	}
	for i, a := range resp {
		if !a.SuccessfullyAnnotated {
			t.Errorf("variant %d: expected successfully_annotated=true", i)
		}
		if a.MostSevereConsequence == "intergenic_variant" {
			t.Errorf("variant %d: unexpected intergenic", i)
		}
	}
}

// --- Error cases ---

func TestGNGenomicBadLocation(t *testing.T) {
	srv := newTestServer(t)
	w := httptest.NewRecorder()
	srv.Handler().ServeHTTP(w, httptest.NewRequest(http.MethodGet, "/genome-nexus/grch38/annotation/genomic/bad", nil))

	if w.Code != http.StatusBadRequest {
		t.Fatalf("expected 400, got %d: %s", w.Code, w.Body.String())
	}
}

func TestGNUnknownAssembly(t *testing.T) {
	srv := newTestServer(t)
	w := httptest.NewRecorder()
	srv.Handler().ServeHTTP(w, httptest.NewRequest(http.MethodGet, "/genome-nexus/hg99/annotation/genomic/7,100,100,A,T", nil))

	if w.Code != http.StatusNotFound {
		t.Fatalf("expected 404, got %d", w.Code)
	}
}

func TestGNGenomicPostEmpty(t *testing.T) {
	srv := newTestServer(t)
	w := httptest.NewRecorder()
	req := httptest.NewRequest(http.MethodPost, "/genome-nexus/grch38/annotation/genomic", strings.NewReader("[]"))
	srv.Handler().ServeHTTP(w, req)

	if w.Code != http.StatusBadRequest {
		t.Fatalf("expected 400, got %d", w.Code)
	}
}

func TestGNGenomicPostBadJSON(t *testing.T) {
	srv := newTestServer(t)
	w := httptest.NewRecorder()
	req := httptest.NewRequest(http.MethodPost, "/genome-nexus/grch38/annotation/genomic", strings.NewReader("not json"))
	srv.Handler().ServeHTTP(w, req)

	if w.Code != http.StatusBadRequest {
		t.Fatalf("expected 400, got %d", w.Code)
	}
}

// --- parseCommaGenomicLocation unit tests ---

func TestParseCommaGenomicLocation(t *testing.T) {
	tests := []struct {
		input string
		want  input.GenomicLocation
	}{
		{
			"7,140753336,140753336,A,T",
			input.GenomicLocation{Chromosome: "7", Start: 140753336, End: 140753336, ReferenceAllele: "A", VariantAllele: "T"},
		},
		{
			"12,25245350,25245350,C,A",
			input.GenomicLocation{Chromosome: "12", Start: 25245350, End: 25245350, ReferenceAllele: "C", VariantAllele: "A"},
		},
	}
	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			got, err := parseCommaGenomicLocation(tt.input)
			if err != nil {
				t.Fatal(err)
			}
			if got != tt.want {
				t.Errorf("got %+v, want %+v", got, tt.want)
			}
		})
	}
}

func TestParseCommaGenomicLocationErrors(t *testing.T) {
	tests := []string{
		"bad",
		"7,abc,140753336,A,T",
		"7,140753336,abc,A,T",
		"7,140753336",
	}
	for _, input := range tests {
		t.Run(input, func(t *testing.T) {
			_, err := parseCommaGenomicLocation(input)
			if err == nil {
				t.Fatal("expected error")
			}
		})
	}
}

// TestGNRoutesAtRoot checks that a server with a single assembly answers the
// genome-nexus API at the root as well as under the assembly prefix, so a
// genome-nexus client can be pointed at it without being reconfigured.
func TestGNRoutesAtRoot(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	for _, path := range []string{
		"/annotation/genomic/12,25245351,25245351,C,A",
		"/genome-nexus/grch38/annotation/genomic/12,25245351,25245351,C,A",
	} {
		req := httptest.NewRequest(http.MethodGet, path, nil)
		w := httptest.NewRecorder()
		handler.ServeHTTP(w, req)

		if w.Code != http.StatusOK {
			t.Fatalf("%s: expected 200, got %d: %s", path, w.Code, w.Body.String())
		}
		var resp output.GNAnnotation
		if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
			t.Fatalf("%s: decode: %v", path, err)
		}
		if resp.MostSevereConsequence != "missense_variant" {
			t.Errorf("%s: most_severe_consequence: got %q, want missense_variant",
				path, resp.MostSevereConsequence)
		}
	}
}

// TestServerEndpointsUnderBothMounts checks that the endpoints describing the
// server, not a variant, answer under the assembly prefix as well as at the
// root. A genome-nexus client configured with a base URL ending in
// /genome-nexus/grch37 asks for .../version there.
func TestServerEndpointsUnderBothMounts(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	for _, path := range []string{"/health", "/info", "/version"} {
		for _, prefix := range []string{"", "/genome-nexus/grch38"} {
			req := httptest.NewRequest(http.MethodGet, prefix+path, nil)
			w := httptest.NewRecorder()
			handler.ServeHTTP(w, req)
			if w.Code != http.StatusOK {
				t.Errorf("%s%s: expected 200, got %d", prefix, path, w.Code)
			}
			var body map[string]any
			if err := json.Unmarshal(w.Body.Bytes(), &body); err != nil {
				t.Errorf("%s%s: decode: %v", prefix, path, err)
			}
		}
	}
}

// TestGNPostDegradesPerItem checks that one unresolvable notation costs the
// caller that entry and nothing else. genome-nexus reports a failure per item;
// aborting the batch would lose every good annotation alongside the bad one.
func TestGNPostDegradesPerItem(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	body := `["12:g.25245351C>A","not-a-variant","12:g.25245351C>A"]`
	req := httptest.NewRequest(http.MethodPost, "/annotation", strings.NewReader(body))
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}
	var resp []output.GNAnnotation
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v\nbody: %s", err, w.Body.String())
	}
	if len(resp) != 3 {
		t.Fatalf("expected one entry per input, got %d", len(resp))
	}
	if !resp[0].SuccessfullyAnnotated || !resp[2].SuccessfullyAnnotated {
		t.Error("the good entries should still be annotated")
	}
	if resp[1].SuccessfullyAnnotated {
		t.Error("the bad entry should be reported as not annotated")
	}
	if resp[1].OriginalVariantQuery != "not-a-variant" {
		t.Errorf("the bad entry should carry its own query, got %q", resp[1].OriginalVariantQuery)
	}
}
