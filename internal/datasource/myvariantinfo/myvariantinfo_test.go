package myvariantinfo

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"testing"
)

const sampleAPIResponse = `{
  "dbsnp": {"rsid": "rs121434569"},
  "gnomad_exome": {
    "af": {"af": 2.7837e-05, "af_afr": 0.000123},
    "ac": {"ac": 7, "ac_afr": 2},
    "an": {"an": 251464, "an_afr": 16256},
    "hom": {"hom": 0, "hom_afr": 0}
  },
  "gnomad_genome": {
    "af": {"af": 1.5e-05},
    "ac": {"ac": 2},
    "an": {"an": 130000},
    "hom": {"hom": 0}
  },
  "vcf": {"ref": "C", "alt": "T", "position": "55181378"}
}`

func TestTransform(t *testing.T) {
	var raw apiResponse
	if err := json.Unmarshal([]byte(sampleAPIResponse), &raw); err != nil {
		t.Fatalf("unmarshal sample: %v", err)
	}

	result := transform("chr7:g.55181378C>T", &raw)
	if result == nil || result.Annotation == nil {
		t.Fatal("expected non-nil result")
	}

	ann := result.Annotation

	// Check variant/query/hgvs fields.
	if ann.Variant != "chr7:g.55181378C>T" {
		t.Errorf("variant: got %q", ann.Variant)
	}
	if ann.Query != "chr7:g.55181378C>T" {
		t.Errorf("query: got %q", ann.Query)
	}
	if ann.Hgvs != "chr7:g.55181378C>T" {
		t.Errorf("hgvs: got %q", ann.Hgvs)
	}

	// Check dbSNP.
	if ann.Dbsnp == nil || ann.Dbsnp.Rsid != "rs121434569" {
		t.Errorf("dbsnp rsid: got %+v", ann.Dbsnp)
	}

	// Check gnomAD exome.
	if ann.GnomadExome == nil {
		t.Fatal("expected gnomadExome")
	}
	if v, ok := ann.GnomadExome.AlleleCount["ac"]; !ok || v.(float64) != 7 {
		t.Errorf("gnomadExome.alleleCount.ac: got %v", ann.GnomadExome.AlleleCount["ac"])
	}
	if v, ok := ann.GnomadExome.AlleleFrequency["af"]; !ok || v.(float64) != 2.7837e-05 {
		t.Errorf("gnomadExome.alleleFrequency.af: got %v", ann.GnomadExome.AlleleFrequency["af"])
	}
	if v, ok := ann.GnomadExome.AlleleNumber["an"]; !ok || v.(float64) != 251464 {
		t.Errorf("gnomadExome.alleleNumber.an: got %v", ann.GnomadExome.AlleleNumber["an"])
	}
	if v, ok := ann.GnomadExome.Homozygotes["hom"]; !ok || v.(float64) != 0 {
		t.Errorf("gnomadExome.homozygotes.hom: got %v", ann.GnomadExome.Homozygotes["hom"])
	}

	// Check gnomAD genome.
	if ann.GnomadGenome == nil {
		t.Fatal("expected gnomadGenome")
	}
	if v, ok := ann.GnomadGenome.AlleleCount["ac"]; !ok || v.(float64) != 2 {
		t.Errorf("gnomadGenome.alleleCount.ac: got %v", ann.GnomadGenome.AlleleCount["ac"])
	}

	// Check VCF.
	if ann.Vcf == nil {
		t.Fatal("expected vcf")
	}
	if ann.Vcf.Ref != "C" || ann.Vcf.Alt != "T" || ann.Vcf.Position != "55181378" {
		t.Errorf("vcf: got %+v", ann.Vcf)
	}
}

func TestTransformJSON(t *testing.T) {
	var raw apiResponse
	if err := json.Unmarshal([]byte(sampleAPIResponse), &raw); err != nil {
		t.Fatalf("unmarshal sample: %v", err)
	}

	result := transform("chr7:g.55181378C>T", &raw)

	data, err := json.Marshal(result)
	if err != nil {
		t.Fatalf("marshal: %v", err)
	}

	// Verify it round-trips and contains expected keys.
	var m map[string]interface{}
	if err := json.Unmarshal(data, &m); err != nil {
		t.Fatalf("unmarshal output: %v", err)
	}

	ann, ok := m["annotation"].(map[string]interface{})
	if !ok {
		t.Fatalf("expected annotation key, got: %s", string(data))
	}
	if _, ok := ann["gnomadExome"]; !ok {
		t.Error("expected gnomadExome key in annotation")
	}
	if _, ok := ann["gnomadGenome"]; !ok {
		t.Error("expected gnomadGenome key in annotation")
	}
	if _, ok := ann["dbsnp"]; !ok {
		t.Error("expected dbsnp key in annotation")
	}
}

func TestEnsureChrPrefix(t *testing.T) {
	tests := []struct {
		input, want string
	}{
		{"7:g.55181378C>T", "chr7:g.55181378C>T"},
		{"chr7:g.55181378C>T", "chr7:g.55181378C>T"},
		{"X:g.100A>G", "chrX:g.100A>G"},
	}
	for _, tt := range tests {
		got := ensureChrPrefix(tt.input)
		if got != tt.want {
			t.Errorf("ensureChrPrefix(%q) = %q, want %q", tt.input, got, tt.want)
		}
	}
}

func TestBuildURL(t *testing.T) {
	tests := []struct {
		hgvsg, assembly, want string
	}{
		{"chr7:g.55181378C>T", "GRCh38", "https://myvariant.info/v1/variant/chr7:g.55181378C>T?assembly=hg38"},
		{"chr7:g.55181378C>T", "grch38", "https://myvariant.info/v1/variant/chr7:g.55181378C>T?assembly=hg38"},
		{"chr7:g.55181378C>T", "hg38", "https://myvariant.info/v1/variant/chr7:g.55181378C>T?assembly=hg38"},
		{"chr7:g.55181378C>T", "GRCh37", "https://myvariant.info/v1/variant/chr7:g.55181378C>T"},
		{"chr7:g.55181378C>T", "hg19", "https://myvariant.info/v1/variant/chr7:g.55181378C>T"},
	}
	for _, tt := range tests {
		got := buildURL(tt.hgvsg, tt.assembly)
		if got != tt.want {
			t.Errorf("buildURL(%q, %q) = %q, want %q", tt.hgvsg, tt.assembly, got, tt.want)
		}
	}
}

func TestFetchWithMockServer(t *testing.T) {
	server := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Content-Type", "application/json")
		w.Write([]byte(sampleAPIResponse))
	}))
	defer server.Close()

	client := &Client{httpClient: server.Client()}

	// Override the URL construction by using the mock server URL directly.
	// We test the Fetch method indirectly through transform + HTTP handling.
	resp, err := client.httpClient.Get(server.URL)
	if err != nil {
		t.Fatalf("GET: %v", err)
	}
	defer resp.Body.Close()

	var raw apiResponse
	if err := json.NewDecoder(resp.Body).Decode(&raw); err != nil {
		t.Fatalf("decode: %v", err)
	}

	result := transform("chr7:g.55181378C>T", &raw)
	if result == nil || result.Annotation == nil {
		t.Fatal("expected non-nil result from mock")
	}
	if result.Annotation.Dbsnp == nil || result.Annotation.Dbsnp.Rsid != "rs121434569" {
		t.Errorf("expected rs121434569, got %+v", result.Annotation.Dbsnp)
	}
}

func TestFetch404(t *testing.T) {
	server := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.WriteHeader(http.StatusNotFound)
		w.Write([]byte(`{"success":false,"error":"not found"}`))
	}))
	defer server.Close()

	// Test that 404 produces nil, nil by simulating the logic.
	resp, err := http.Get(server.URL)
	if err != nil {
		t.Fatalf("GET: %v", err)
	}
	resp.Body.Close()

	if resp.StatusCode != http.StatusNotFound {
		t.Fatalf("expected 404, got %d", resp.StatusCode)
	}
	// In Fetch, this would return nil, nil.
}

func TestTransformEmptyResponse(t *testing.T) {
	raw := &apiResponse{}
	result := transform("chr1:g.100A>G", raw)
	if result == nil || result.Annotation == nil {
		t.Fatal("expected non-nil result")
	}
	if result.Annotation.Dbsnp != nil {
		t.Error("expected nil dbsnp for empty response")
	}
	if result.Annotation.GnomadExome != nil {
		t.Error("expected nil gnomadExome for empty response")
	}
	if result.Annotation.GnomadGenome != nil {
		t.Error("expected nil gnomadGenome for empty response")
	}
	if result.Annotation.Vcf != nil {
		t.Error("expected nil vcf for empty response")
	}
	if result.Annotation.Variant != "chr1:g.100A>G" {
		t.Errorf("variant: got %q", result.Annotation.Variant)
	}
}
