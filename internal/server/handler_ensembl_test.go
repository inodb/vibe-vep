package server

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"go.uber.org/zap"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/output"
)

// findTestCacheDir locates the testdata/cache directory.
func findTestCacheDir(t *testing.T) string {
	t.Helper()
	paths := []string{
		filepath.Join("testdata", "cache"),
		filepath.Join("..", "..", "testdata", "cache"),
	}
	for _, p := range paths {
		if _, err := os.Stat(p); err == nil {
			return p
		}
	}
	t.Fatal("Test cache directory not found")
	return ""
}

// newTestServer creates a server with an empty cache (all variants intergenic).
func newTestServer(t *testing.T) *Server {
	t.Helper()
	logger := zap.NewNop()
	srv := New(logger, "test")

	c := cache.New()
	ann := annotate.NewAnnotator(c)
	srv.AddAssembly("GRCh38", c, ann, nil)

	return srv
}

// newTestServerWithKRAS creates a server with KRAS transcripts loaded on chr12.
func newTestServerWithKRAS(t *testing.T) *Server {
	t.Helper()
	logger := zap.NewNop()
	srv := New(logger, "test")

	cacheDir := findTestCacheDir(t)
	c := cache.New()
	loader := cache.NewLoader(cacheDir, "homo_sapiens", "GRCh38")
	if err := loader.Load(c, "12"); err != nil {
		t.Fatalf("loading chr12 transcripts: %v", err)
	}
	c.BuildIndex()

	ann := annotate.NewAnnotator(c)
	srv.AddAssembly("GRCh38", c, ann, nil)

	return srv
}

// --- Health & Info ---

func TestHealthEndpoint(t *testing.T) {
	srv := newTestServer(t)
	w := httptest.NewRecorder()
	srv.Handler().ServeHTTP(w, httptest.NewRequest(http.MethodGet, "/health", nil))

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d", w.Code)
	}
	var resp map[string]string
	json.Unmarshal(w.Body.Bytes(), &resp)
	if resp["status"] != "ok" {
		t.Fatalf("expected status ok, got %q", resp["status"])
	}
}

func TestInfoEndpoint(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	w := httptest.NewRecorder()
	srv.Handler().ServeHTTP(w, httptest.NewRequest(http.MethodGet, "/info", nil))

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d", w.Code)
	}
	var resp map[string]any
	json.Unmarshal(w.Body.Bytes(), &resp)
	if resp["version"] != "test" {
		t.Fatalf("expected version test, got %v", resp["version"])
	}
	assemblies := resp["assemblies"].([]any)
	if len(assemblies) != 1 {
		t.Fatalf("expected 1 assembly, got %d", len(assemblies))
	}
	asm := assemblies[0].(map[string]any)
	if tc := asm["transcript_count"].(float64); tc == 0 {
		t.Fatal("expected non-zero transcript count")
	}
}

// --- Ensembl Region GET ---

func TestEnsemblRegionGet_InKRAS(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	// Region GET only provides the alt allele (ref is inferred from reference genome
	// in the real Ensembl VEP). Without a reference genome, ref is empty, so the
	// consequence may differ from missense. We validate that the endpoint returns
	// a valid annotation hitting the KRAS gene.
	req := httptest.NewRequest(http.MethodGet, "/ensembl/grch38/vep/human/region/12:25245351-25245351:1/A", nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var resp []output.VEPVariantAnnotation
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v\nbody: %s", err, w.Body.String())
	}
	if len(resp) != 1 {
		t.Fatalf("expected 1 variant annotation, got %d", len(resp))
	}

	va := resp[0]
	if va.AssemblyName != "GRCh38" {
		t.Errorf("assembly: got %q, want GRCh38", va.AssemblyName)
	}
	if va.SeqRegionName != "12" {
		t.Errorf("seq_region_name: got %q, want 12", va.SeqRegionName)
	}
	// Should NOT be intergenic — we're inside the KRAS gene.
	if va.MostSevereConsequence == "intergenic_variant" {
		t.Error("expected non-intergenic consequence for KRAS position")
	}
	if len(va.TranscriptConsequences) == 0 {
		t.Fatal("expected transcript consequences, got none")
	}

	// Find the canonical transcript.
	var canonical *output.VEPTranscriptConsequence
	for i := range va.TranscriptConsequences {
		if va.TranscriptConsequences[i].TranscriptID == "ENST00000311936" {
			canonical = &va.TranscriptConsequences[i]
			break
		}
	}
	if canonical == nil {
		t.Fatal("canonical transcript ENST00000311936 not found in response")
	}
	if canonical.GeneSymbol != "KRAS" {
		t.Errorf("gene_symbol: got %q, want KRAS", canonical.GeneSymbol)
	}
	if canonical.Impact == "" {
		t.Error("expected non-empty impact")
	}
}

func TestEnsemblRegionGet_Intergenic(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	// Position far from KRAS — intergenic.
	req := httptest.NewRequest(http.MethodGet, "/ensembl/grch38/vep/human/region/12:1000000-1000000:1/T", nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var resp []output.VEPVariantAnnotation
	json.Unmarshal(w.Body.Bytes(), &resp)
	if resp[0].MostSevereConsequence != "intergenic_variant" {
		t.Errorf("expected intergenic_variant, got %q", resp[0].MostSevereConsequence)
	}
}

// --- Ensembl HGVS GET ---

func TestEnsemblHGVSGet_KRASG12C(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	// HGVSc notation resolves to genomic coordinates including ref allele.
	req := httptest.NewRequest(http.MethodGet, "/ensembl/grch38/vep/human/hgvs/ENST00000311936:c.34G>T", nil)
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var resp []output.VEPVariantAnnotation
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v\nbody: %s", err, w.Body.String())
	}
	if len(resp) != 1 {
		t.Fatalf("expected 1 annotation, got %d", len(resp))
	}

	va := resp[0]
	if va.MostSevereConsequence != "missense_variant" {
		t.Errorf("most_severe_consequence: got %q, want missense_variant", va.MostSevereConsequence)
	}
	if va.AssemblyName != "GRCh38" {
		t.Errorf("assembly: got %q, want GRCh38", va.AssemblyName)
	}

	// Find canonical transcript.
	var canonical *output.VEPTranscriptConsequence
	for i := range va.TranscriptConsequences {
		if va.TranscriptConsequences[i].TranscriptID == "ENST00000311936" {
			canonical = &va.TranscriptConsequences[i]
			break
		}
	}
	if canonical == nil {
		t.Fatal("ENST00000311936 not found")
	}
	if canonical.GeneSymbol != "KRAS" {
		t.Errorf("gene_symbol: got %q, want KRAS", canonical.GeneSymbol)
	}
	if canonical.Impact != "MODERATE" {
		t.Errorf("impact: got %q, want MODERATE", canonical.Impact)
	}
	if !contains(canonical.ConsequenceTerms, "missense_variant") {
		t.Errorf("consequence_terms: %v missing missense_variant", canonical.ConsequenceTerms)
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
}

func TestEnsemblHGVSPost_Batch(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	body := `{"hgvs_notations": ["ENST00000311936:c.34G>T", "ENST00000311936:c.35G>A"]}`
	req := httptest.NewRequest(http.MethodPost, "/ensembl/grch38/vep/human/hgvs", strings.NewReader(body))
	req.Header.Set("Content-Type", "application/json")
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var resp []output.VEPVariantAnnotation
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v", err)
	}
	if len(resp) != 2 {
		t.Fatalf("expected 2 annotations, got %d", len(resp))
	}
	for i, va := range resp {
		if va.MostSevereConsequence == "" {
			t.Errorf("variant %d: empty most_severe_consequence", i)
		}
		if va.MostSevereConsequence == "intergenic_variant" {
			t.Errorf("variant %d: unexpected intergenic", i)
		}
	}
}

// --- Ensembl Region POST ---

func TestEnsemblRegionPost_Batch(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	body := `{"variants": ["12:25245351-25245351:1/A", "12:25245350-25245350:1/T"]}`
	req := httptest.NewRequest(http.MethodPost, "/ensembl/grch38/vep/human/region", strings.NewReader(body))
	req.Header.Set("Content-Type", "application/json")
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)

	if w.Code != http.StatusOK {
		t.Fatalf("expected 200, got %d: %s", w.Code, w.Body.String())
	}

	var resp []output.VEPVariantAnnotation
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatalf("decode: %v", err)
	}
	if len(resp) != 2 {
		t.Fatalf("expected 2 annotations, got %d", len(resp))
	}
	// Both should be in KRAS.
	for i, va := range resp {
		if va.MostSevereConsequence == "" {
			t.Errorf("variant %d: empty most_severe_consequence", i)
		}
		if len(va.TranscriptConsequences) == 0 {
			t.Errorf("variant %d: no transcript consequences", i)
		}
	}
}

// --- Error cases ---

func TestEnsemblUnknownAssembly(t *testing.T) {
	srv := newTestServer(t)
	w := httptest.NewRecorder()
	srv.Handler().ServeHTTP(w, httptest.NewRequest(http.MethodGet, "/ensembl/hg99/vep/human/region/7:100-100:1/T", nil))

	if w.Code != http.StatusNotFound {
		t.Fatalf("expected 404, got %d", w.Code)
	}
}

func TestEnsemblRegionBadFormat(t *testing.T) {
	srv := newTestServer(t)
	w := httptest.NewRecorder()
	srv.Handler().ServeHTTP(w, httptest.NewRequest(http.MethodGet, "/ensembl/grch38/vep/human/region/bad/T", nil))

	if w.Code != http.StatusBadRequest {
		t.Fatalf("expected 400, got %d: %s", w.Code, w.Body.String())
	}
}

func TestEnsemblRegionPostBadJSON(t *testing.T) {
	srv := newTestServer(t)
	w := httptest.NewRecorder()
	req := httptest.NewRequest(http.MethodPost, "/ensembl/grch38/vep/human/region", strings.NewReader("not json"))
	srv.Handler().ServeHTTP(w, req)

	if w.Code != http.StatusBadRequest {
		t.Fatalf("expected 400, got %d", w.Code)
	}
}

func TestEnsemblRegionPostEmpty(t *testing.T) {
	srv := newTestServer(t)
	w := httptest.NewRecorder()
	req := httptest.NewRequest(http.MethodPost, "/ensembl/grch38/vep/human/region", strings.NewReader(`{"variants":[]}`))
	srv.Handler().ServeHTTP(w, req)

	if w.Code != http.StatusBadRequest {
		t.Fatalf("expected 400, got %d", w.Code)
	}
}

// --- parseEnsemblRegion unit tests ---

func TestParseEnsemblRegion(t *testing.T) {
	tests := []struct {
		region string
		allele string
		chrom  string
		pos    int64
		alt    string
	}{
		{"7:140753336-140753336:1", "T", "7", 140753336, "T"},
		{"chr12:25245350-25245350:1", "A", "12", 25245350, "A"},
		{"1:100-100", "G", "1", 100, "G"},
	}
	for _, tt := range tests {
		t.Run(tt.region, func(t *testing.T) {
			v, err := parseEnsemblRegion(tt.region, tt.allele)
			if err != nil {
				t.Fatal(err)
			}
			if v.Chrom != tt.chrom {
				t.Errorf("chrom: got %q, want %q", v.Chrom, tt.chrom)
			}
			if v.Pos != tt.pos {
				t.Errorf("pos: got %d, want %d", v.Pos, tt.pos)
			}
			if v.Alt != tt.alt {
				t.Errorf("alt: got %q, want %q", v.Alt, tt.alt)
			}
		})
	}
}

func TestParseEnsemblRegionErrors(t *testing.T) {
	tests := []struct {
		name   string
		region string
	}{
		{"too many colons", "1:2:3:4:5"},
		{"no range", "chr1"},
		{"bad start", "1:abc-100:1"},
		{"bad end", "1:100-abc:1"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			_, err := parseEnsemblRegion(tt.region, "A")
			if err == nil {
				t.Fatal("expected error")
			}
		})
	}
}

func contains(ss []string, s string) bool {
	for _, v := range ss {
		if v == s {
			return true
		}
	}
	return false
}
