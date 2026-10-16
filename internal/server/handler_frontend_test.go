package server

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"

	"github.com/inodb/vibe-vep/internal/datasource/pfam"
	"github.com/inodb/vibe-vep/internal/output"
)

// TestFieldsParsingRepeated tests that ?fields=a&fields=b works (frontend pattern).
func TestFieldsParsingRepeated(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	// Frontend sends each field as a separate query param
	req := httptest.NewRequest(http.MethodGet,
		"/genome-nexus/grch38/annotation/genomic/12,25245351,25245351,C,A?fields=annotation_summary&fields=clinvar&fields=hotspots",
		nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var resp output.GNAnnotation
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v", err)
	}
	if resp.AnnotationSummary == nil {
		t.Error("expected annotation_summary to be present with repeated ?fields= params")
	}
}

// TestFieldsParsingCommaSeparated tests that ?fields=a,b works.
func TestFieldsParsingCommaSeparated(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	req := httptest.NewRequest(http.MethodGet,
		"/genome-nexus/grch38/annotation/genomic/12,25245351,25245351,C,A?fields=annotation_summary,clinvar,hotspots",
		nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var resp output.GNAnnotation
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v", err)
	}
	if resp.AnnotationSummary == nil {
		t.Error("expected annotation_summary to be present with comma-separated fields")
	}
}

// TestAnnotationSummaryFields tests the annotation_summary response has required fields.
func TestAnnotationSummaryFields(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	req := httptest.NewRequest(http.MethodGet,
		"/genome-nexus/grch38/annotation/genomic/12,25245351,25245351,C,A?fields=annotation_summary",
		nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	var resp output.GNAnnotation
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v", err)
	}

	s := resp.AnnotationSummary
	if s == nil {
		t.Fatal("annotation_summary is nil")
	}
	if s.VariantType != "SNP" {
		t.Errorf("variantType: got %q, want SNP", s.VariantType)
	}
	if s.AssemblyName != "GRCh38" {
		t.Errorf("assemblyName: got %q, want GRCh38", s.AssemblyName)
	}
	if s.TranscriptConsequenceSummary == nil {
		t.Fatal("transcriptConsequenceSummary is nil")
	}
	tcs := s.TranscriptConsequenceSummary
	if tcs.HugoGeneSymbol != "KRAS" {
		t.Errorf("hugoGeneSymbol: got %q, want KRAS", tcs.HugoGeneSymbol)
	}
	if tcs.VariantClassification != "Missense_Mutation" {
		t.Errorf("variantClassification: got %q, want Missense_Mutation", tcs.VariantClassification)
	}
	if tcs.ConsequenceTerms != "missense_variant" {
		t.Errorf("consequenceTerms: got %q, want missense_variant", tcs.ConsequenceTerms)
	}
}

// TestCanonicalTranscriptFirst tests that transcript_consequences[0] is canonical.
func TestCanonicalTranscriptFirst(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	req := httptest.NewRequest(http.MethodGet,
		"/genome-nexus/grch38/annotation/genomic/12,25245351,25245351,C,A",
		nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	var resp output.GNAnnotation
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v", err)
	}

	if len(resp.TranscriptConsequences) == 0 {
		t.Fatal("no transcript consequences")
	}
	if resp.TranscriptConsequences[0].Canonical != "1" {
		t.Error("expected transcript_consequences[0] to be canonical")
	}
}

// TestGNErrorIncludesAssembly tests that error responses include assembly_name.
func TestGNErrorIncludesAssembly(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	// Use a variant notation that can't be resolved (insertion needing ref base lookup)
	req := httptest.NewRequest(http.MethodGet,
		"/genome-nexus/grch38/annotation/99:g.1_2insACGT",
		nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	var resp map[string]interface{}
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v", err)
	}
	if resp["assembly_name"] != "GRCh38" {
		t.Errorf("expected assembly_name=GRCh38 in error response, got %v", resp["assembly_name"])
	}
	if sa, ok := resp["successfully_annotated"].(bool); !ok || sa != false {
		t.Errorf("expected successfully_annotated=false in error response, got %v", resp["successfully_annotated"])
	}
}

// TestPfamDomainPost tests POST /pfam/domain with a JSON array of accessions.
func TestPfamDomainPost(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	// Add a PFAM store with test data
	store := pfam.NewEmptyStore()
	store.AddDomain("PF00069", "Pkinase", "Protein kinase domain")
	store.AddDomain("PF07714", "Pkinase_Tyr", "Protein tyrosine kinase")
	srv.SetPfamStore("GRCh38", store)

	handler := srv.Handler()

	body := `["PF00069","PF07714"]`
	req := httptest.NewRequest(http.MethodPost,
		"/genome-nexus/grch38/pfam/domain",
		strings.NewReader(body))
	req.Header.Set("Content-Type", "application/json")
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var domains []map[string]string
	if err := json.Unmarshal(w.Body.Bytes(), &domains); err != nil {
		t.Fatalf("decode: %v", err)
	}
	if len(domains) != 2 {
		t.Fatalf("expected 2 domains, got %d", len(domains))
	}
	if domains[0]["name"] != "Pkinase" {
		t.Errorf("domain[0].name: got %q, want Pkinase", domains[0]["name"])
	}
}

// TestPfamDomainGet tests GET /pfam/domain/{accession}.
func TestPfamDomainGet(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	store := pfam.NewEmptyStore()
	store.AddDomain("PF07714", "Pkinase_Tyr", "Protein tyrosine kinase")
	srv.SetPfamStore("GRCh38", store)

	handler := srv.Handler()

	req := httptest.NewRequest(http.MethodGet,
		"/genome-nexus/grch38/pfam/domain/PF07714",
		nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var domain map[string]string
	if err := json.Unmarshal(w.Body.Bytes(), &domain); err != nil {
		t.Fatalf("decode: %v", err)
	}
	if domain["pfamAccession"] != "PF07714" {
		t.Errorf("pfamAccession: got %q", domain["pfamAccession"])
	}
	if domain["name"] != "Pkinase_Tyr" {
		t.Errorf("name: got %q", domain["name"])
	}
}

// TestSignalFieldIncluded tests that ?fields=signal returns signalAnnotation.
func TestSignalFieldIncluded(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	// With signal field
	req := httptest.NewRequest(http.MethodGet,
		"/genome-nexus/grch38/annotation/genomic/12,25245351,25245351,C,A?fields=signal",
		nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var resp map[string]json.RawMessage
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v", err)
	}

	// signalAnnotation should be present (even if empty annotation array)
	sigRaw, ok := resp["signalAnnotation"]
	if !ok {
		t.Fatal("expected signalAnnotation in response when ?fields=signal")
	}

	var sig output.GNSignalAnnotation
	if err := json.Unmarshal(sigRaw, &sig); err != nil {
		t.Fatalf("decode signalAnnotation: %v", err)
	}
	if sig.License == "" {
		t.Error("expected non-empty license in signalAnnotation")
	}
	// Annotation should be an array (possibly empty)
	if sig.Annotation == nil {
		t.Error("expected annotation array, got nil")
	}
}

// TestSignalFieldNotIncludedByDefault tests that signalAnnotation is absent without ?fields=signal.
func TestSignalFieldNotIncludedByDefault(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	req := httptest.NewRequest(http.MethodGet,
		"/genome-nexus/grch38/annotation/genomic/12,25245351,25245351,C,A",
		nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	var resp map[string]json.RawMessage
	json.Unmarshal(w.Body.Bytes(), &resp)

	if _, ok := resp["signalAnnotation"]; ok {
		t.Error("signalAnnotation should not be present without ?fields=signal")
	}
}

// TestSignalFieldWithRepeatedParams tests signal with repeated ?fields= params.
func TestSignalFieldWithRepeatedParams(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	req := httptest.NewRequest(http.MethodGet,
		"/genome-nexus/grch38/annotation/genomic/12,25245351,25245351,C,A?fields=annotation_summary&fields=signal&fields=clinvar",
		nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	var resp map[string]json.RawMessage
	json.Unmarshal(w.Body.Bytes(), &resp)

	if _, ok := resp["signalAnnotation"]; !ok {
		t.Error("expected signalAnnotation with repeated ?fields= including signal")
	}
	if _, ok := resp["annotation_summary"]; !ok {
		t.Error("expected annotation_summary with repeated ?fields=")
	}
}

// TestEnsemblTranscriptFilterIsAlternativeNotUnion pins genome-nexus's actual
// semantics: the four lists in the POST body are alternatives, taken in the
// order transcriptIds, geneIds, proteinIds, hugoSymbols, and the first
// non-empty one wins. Sending two of them must not union the results.
func TestEnsemblTranscriptFilterIsAlternativeNotUnion(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	post := func(body string) []map[string]any {
		t.Helper()
		req := httptest.NewRequest(http.MethodPost, "/ensembl/transcript", strings.NewReader(body))
		w := httptest.NewRecorder()
		handler.ServeHTTP(w, req)
		if w.Code != http.StatusOK {
			t.Fatalf("%s: expected 200, got %d", body, w.Code)
		}
		var out []map[string]any
		if err := json.Unmarshal(w.Body.Bytes(), &out); err != nil {
			t.Fatalf("%s: decode: %v (body %s)", body, err, w.Body.String())
		}
		return out
	}

	byGene := post(`{"hugoSymbols":["KRAS"]}`)
	if len(byGene) == 0 {
		t.Fatal("expected transcripts for KRAS")
	}
	byTx := post(`{"transcriptIds":["ENST00000311936"]}`)
	if len(byTx) != 1 {
		t.Fatalf("expected exactly one transcript, got %d", len(byTx))
	}

	// transcriptIds outranks hugoSymbols, so the combined request must equal
	// the transcriptIds-only request, not the two concatenated.
	both := post(`{"transcriptIds":["ENST00000311936"],"hugoSymbols":["KRAS"]}`)
	if len(both) != len(byTx) {
		t.Errorf("combined filter returned %d, want %d (transcriptIds should win)", len(both), len(byTx))
	}

	// An empty filter is an empty list, and empty must serialise as [].
	req := httptest.NewRequest(http.MethodPost, "/ensembl/transcript", strings.NewReader(`{}`))
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)
	if got := strings.TrimSpace(w.Body.String()); got != "[]" {
		t.Errorf("empty filter body = %q, want []", got)
	}
}

// TestEnsemblTranscriptGetIsConjunctive pins the GET form, where the
// parameters narrow each other rather than being alternatives, and no
// parameters at all means an empty list.
func TestEnsemblTranscriptGetIsConjunctive(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	get := func(query string) []map[string]any {
		t.Helper()
		req := httptest.NewRequest(http.MethodGet, "/ensembl/transcript"+query, nil)
		w := httptest.NewRecorder()
		handler.ServeHTTP(w, req)
		if w.Code != http.StatusOK {
			t.Fatalf("%s: expected 200, got %d", query, w.Code)
		}
		var out []map[string]any
		if err := json.Unmarshal(w.Body.Bytes(), &out); err != nil {
			t.Fatalf("%s: decode: %v", query, err)
		}
		return out
	}

	if n := len(get("")); n != 0 {
		t.Errorf("no parameters returned %d transcripts, want 0", n)
	}
	byGene := get("?hugoSymbol=KRAS")
	if len(byGene) == 0 {
		t.Fatal("expected transcripts for KRAS")
	}
	// A gene that does not carry that symbol narrows to nothing.
	if n := len(get("?hugoSymbol=KRAS&geneId=ENSG00000000000")); n != 0 {
		t.Errorf("contradictory filters returned %d, want 0", n)
	}
}

// TestCanonicalBatchEndpoints checks the three batch lookups answer in the
// same shape as their single-value forms, skip identifiers they cannot
// resolve rather than failing, and serialise empty as [].
func TestCanonicalBatchEndpoints(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	post := func(path, body string) ([]map[string]any, string) {
		t.Helper()
		req := httptest.NewRequest(http.MethodPost, path, strings.NewReader(body))
		w := httptest.NewRecorder()
		handler.ServeHTTP(w, req)
		if w.Code != http.StatusOK {
			t.Fatalf("%s %s: expected 200, got %d", path, body, w.Code)
		}
		raw := strings.TrimSpace(w.Body.String())
		var out []map[string]any
		if err := json.Unmarshal([]byte(raw), &out); err != nil {
			t.Fatalf("%s: decode: %v (body %s)", path, err, raw)
		}
		return out, raw
	}

	// A known gene resolves; an unknown one is skipped, not fatal.
	got, _ := post("/ensembl/canonical-transcript/hgnc", `["KRAS","NOT_A_GENE"]`)
	if len(got) != 1 {
		t.Fatalf("canonical-transcript batch returned %d, want 1", len(got))
	}
	if got[0]["transcriptId"] == "" {
		t.Error("expected a transcript id in the response")
	}

	got, _ = post("/ensembl/canonical-gene/hgnc", `["KRAS"]`)
	if len(got) != 1 || got[0]["hugoSymbol"] != "KRAS" {
		t.Errorf("canonical-gene batch = %v, want one KRAS entry", got)
	}

	// Nothing resolvable still answers [], never null.
	_, raw := post("/ensembl/canonical-gene/entrez", `["0"]`)
	if raw != "[]" {
		t.Errorf("entrez batch with no matches = %q, want []", raw)
	}
}
