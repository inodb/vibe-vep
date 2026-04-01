// Package myvariantinfo provides a client for the myvariant.info API,
// returning gnomAD population frequencies and dbSNP rsid data in
// genome-nexus compatible format.
package myvariantinfo

import (
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"strings"
	"time"
)

// Client queries the myvariant.info API.
type Client struct {
	httpClient *http.Client
}

// NewClient creates a new myvariant.info client with a 10-second timeout.
func NewClient() *Client {
	return &Client{
		httpClient: &http.Client{Timeout: 10 * time.Second},
	}
}

// MyVariantInfoResponse wraps the annotation in genome-nexus format.
type MyVariantInfoResponse struct {
	Annotation *MyVariantInfo `json:"annotation,omitempty"`
}

// MyVariantInfo holds the transformed myvariant.info data.
type MyVariantInfo struct {
	Dbsnp        *Dbsnp                 `json:"dbsnp,omitempty"`
	GnomadExome  *Gnomad                `json:"gnomadExome,omitempty"`
	GnomadGenome *Gnomad                `json:"gnomadGenome,omitempty"`
	Vcf          *Vcf                   `json:"vcf,omitempty"`
	Variant      string                 `json:"variant,omitempty"`
	Query        string                 `json:"query,omitempty"`
	Hgvs         string                 `json:"hgvs,omitempty"`
}

// Dbsnp holds dbSNP rsid.
type Dbsnp struct {
	Rsid string `json:"rsid,omitempty"`
}

// Gnomad holds gnomAD frequency data with the nested allele count/frequency/number/homozygotes maps.
type Gnomad struct {
	AlleleCount     map[string]interface{} `json:"alleleCount,omitempty"`
	AlleleFrequency map[string]interface{} `json:"alleleFrequency,omitempty"`
	AlleleNumber    map[string]interface{} `json:"alleleNumber,omitempty"`
	Homozygotes     map[string]interface{} `json:"homozygotes,omitempty"`
}

// Vcf holds VCF-level variant info.
type Vcf struct {
	Ref      string `json:"ref,omitempty"`
	Alt      string `json:"alt,omitempty"`
	Position string `json:"position,omitempty"`
}

// apiResponse represents the raw JSON structure from myvariant.info.
type apiResponse struct {
	Dbsnp        *apiDbsnp       `json:"dbsnp,omitempty"`
	GnomadExome  *apiGnomad      `json:"gnomad_exome,omitempty"`
	GnomadGenome *apiGnomad      `json:"gnomad_genome,omitempty"`
	Vcf          *apiVcf         `json:"vcf,omitempty"`
}

type apiDbsnp struct {
	Rsid string `json:"rsid,omitempty"`
}

type apiGnomad struct {
	AF  map[string]interface{} `json:"af,omitempty"`
	AC  map[string]interface{} `json:"ac,omitempty"`
	AN  map[string]interface{} `json:"an,omitempty"`
	Hom map[string]interface{} `json:"hom,omitempty"`
}

type apiVcf struct {
	Ref      string `json:"ref,omitempty"`
	Alt      string `json:"alt,omitempty"`
	Position string `json:"position,omitempty"`
}

// Fetch queries myvariant.info for the given HGVSg notation and returns
// the transformed response. Returns nil (not error) if the variant is not found (404).
func (c *Client) Fetch(hgvsg, assembly string) (*MyVariantInfoResponse, error) {
	hgvsg = ensureChrPrefix(hgvsg)

	url := buildURL(hgvsg, assembly)

	resp, err := c.httpClient.Get(url)
	if err != nil {
		return nil, fmt.Errorf("myvariant.info request: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode == http.StatusNotFound {
		return nil, nil
	}
	if resp.StatusCode != http.StatusOK {
		body, _ := io.ReadAll(resp.Body)
		return nil, fmt.Errorf("myvariant.info returned %d: %s", resp.StatusCode, string(body))
	}

	var raw apiResponse
	if err := json.NewDecoder(resp.Body).Decode(&raw); err != nil {
		return nil, fmt.Errorf("myvariant.info decode: %w", err)
	}

	return transform(hgvsg, &raw), nil
}

// ensureChrPrefix prepends "chr" to the HGVSg if it does not already have it.
func ensureChrPrefix(hgvsg string) string {
	if strings.HasPrefix(hgvsg, "chr") {
		return hgvsg
	}
	return "chr" + hgvsg
}

// buildURL constructs the myvariant.info API URL.
func buildURL(hgvsg, assembly string) string {
	assembly = strings.ToUpper(assembly)
	if assembly == "GRCH38" || assembly == "HG38" {
		return "https://myvariant.info/v1/variant/" + hgvsg + "?assembly=hg38"
	}
	// GRCh37 is the default for myvariant.info.
	return "https://myvariant.info/v1/variant/" + hgvsg
}

// transform converts the raw API response into genome-nexus format.
func transform(hgvsg string, raw *apiResponse) *MyVariantInfoResponse {
	info := &MyVariantInfo{
		Variant: hgvsg,
		Query:   hgvsg,
		Hgvs:    hgvsg,
	}

	if raw.Dbsnp != nil && raw.Dbsnp.Rsid != "" {
		info.Dbsnp = &Dbsnp{Rsid: raw.Dbsnp.Rsid}
	}

	if raw.GnomadExome != nil {
		info.GnomadExome = transformGnomad(raw.GnomadExome)
	}
	if raw.GnomadGenome != nil {
		info.GnomadGenome = transformGnomad(raw.GnomadGenome)
	}

	if raw.Vcf != nil {
		info.Vcf = &Vcf{
			Ref:      raw.Vcf.Ref,
			Alt:      raw.Vcf.Alt,
			Position: raw.Vcf.Position,
		}
	}

	return &MyVariantInfoResponse{Annotation: info}
}

// transformGnomad maps the raw gnomAD structure to the genome-nexus field names.
func transformGnomad(raw *apiGnomad) *Gnomad {
	if raw == nil {
		return nil
	}
	g := &Gnomad{}
	if len(raw.AC) > 0 {
		g.AlleleCount = raw.AC
	}
	if len(raw.AF) > 0 {
		g.AlleleFrequency = raw.AF
	}
	if len(raw.AN) > 0 {
		g.AlleleNumber = raw.AN
	}
	if len(raw.Hom) > 0 {
		g.Homozygotes = raw.Hom
	}
	return g
}
