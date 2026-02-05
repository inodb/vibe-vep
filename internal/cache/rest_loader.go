// Package cache provides VEP cache loading functionality.
package cache

import (
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"sort"
	"time"
)

// RESTLoader loads transcript data from Ensembl REST API.
// This is useful when the local VEP cache is unavailable or in an unsupported format.
type RESTLoader struct {
	baseURL    string
	assembly   string
	httpClient *http.Client
}

// NewRESTLoader creates a new REST API loader.
// assembly should be "GRCh37" or "GRCh38".
func NewRESTLoader(assembly string) *RESTLoader {
	baseURL := "https://rest.ensembl.org"
	if assembly == "GRCh37" {
		baseURL = "https://grch37.rest.ensembl.org"
	}

	return &RESTLoader{
		baseURL:  baseURL,
		assembly: assembly,
		httpClient: &http.Client{
			Timeout: 30 * time.Second,
		},
	}
}

// Load loads all transcripts for a chromosome (not practical for REST API).
func (l *RESTLoader) Load(c *Cache, chrom string) error {
	return fmt.Errorf("REST loader does not support loading entire chromosomes; use LoadRegion instead")
}

// LoadAll loads all transcripts (not supported for REST API).
func (l *RESTLoader) LoadAll(c *Cache) error {
	return fmt.Errorf("REST loader does not support loading all transcripts; use LoadRegion instead")
}

// LoadRegion loads transcripts overlapping a specific genomic region from the REST API.
func (l *RESTLoader) LoadRegion(c *Cache, chrom string, start, end int64) error {
	// Check if we already have transcripts for this region
	existing := c.Lookup(chrom, start)
	if len(existing) > 0 {
		return nil // Already loaded
	}

	url := fmt.Sprintf("%s/overlap/region/human/%s:%d-%d?feature=transcript;content-type=application/json",
		l.baseURL, chrom, start, end)

	resp, err := l.httpClient.Get(url)
	if err != nil {
		return fmt.Errorf("REST API request failed: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		body, _ := io.ReadAll(resp.Body)
		return fmt.Errorf("REST API error %d: %s", resp.StatusCode, string(body))
	}

	var rawTranscripts []restTranscript
	if err := json.NewDecoder(resp.Body).Decode(&rawTranscripts); err != nil {
		return fmt.Errorf("decode REST response: %w", err)
	}

	// Convert and add transcripts to cache
	for _, rt := range rawTranscripts {
		t := rt.toTranscript()
		if t != nil {
			// Fetch additional details (CDS, protein sequence) if it's a coding transcript
			if t.Biotype == "protein_coding" {
				l.fetchTranscriptDetails(t)
			}
			c.AddTranscript(t)
		}
	}

	return nil
}

// restTranscript represents the JSON response from Ensembl REST API overlap endpoint.
type restTranscript struct {
	ID              string `json:"transcript_id"`
	GeneID          string `json:"gene_id"`
	ExternalName    string `json:"external_name"`
	Start           int64  `json:"start"`
	End             int64  `json:"end"`
	Strand          int    `json:"strand"`
	Biotype         string `json:"biotype"`
	IsCanonical     int    `json:"is_canonical"`
	SeqRegionName   string `json:"seq_region_name"`
	AssemblyName    string `json:"assembly_name"`
	Source          string `json:"source"`
	TranscriptIDVer string `json:"id"`
}

func (rt *restTranscript) toTranscript() *Transcript {
	if rt.ID == "" {
		return nil
	}

	return &Transcript{
		ID:          rt.ID,
		GeneID:      rt.GeneID,
		GeneName:    rt.ExternalName,
		Chrom:       rt.SeqRegionName,
		Start:       rt.Start,
		End:         rt.End,
		Strand:      int8(rt.Strand),
		Biotype:     rt.Biotype,
		IsCanonical: rt.IsCanonical == 1,
	}
}

// fetchTranscriptDetails fetches CDS, exons, and protein sequence for a transcript.
func (l *RESTLoader) fetchTranscriptDetails(t *Transcript) {
	// Fetch full transcript info including exons
	lookupURL := fmt.Sprintf("%s/lookup/id/%s?expand=1;content-type=application/json",
		l.baseURL, t.ID)

	if resp, err := l.httpClient.Get(lookupURL); err == nil {
		defer resp.Body.Close()
		if resp.StatusCode == http.StatusOK {
			var lookup struct {
				Translation *struct {
					Start  int64  `json:"start"`
					End    int64  `json:"end"`
					Length int    `json:"length"`
				} `json:"Translation"`
				Exon []struct {
					Start int64 `json:"start"`
					End   int64 `json:"end"`
					Rank  int   `json:"rank"`
				} `json:"Exon"`
			}
			if err := json.NewDecoder(resp.Body).Decode(&lookup); err == nil {
				// Set CDS boundaries from translation
				if lookup.Translation != nil {
					t.CDSStart = lookup.Translation.Start
					t.CDSEnd = lookup.Translation.End
				}

				// Convert exons
				if len(lookup.Exon) > 0 {
					t.Exons = make([]Exon, len(lookup.Exon))
					for i, e := range lookup.Exon {
						t.Exons[i] = Exon{
							Number: e.Rank,
							Start:  e.Start,
							End:    e.End,
							Frame:  -1, // Will be calculated
						}
						// Calculate CDS portion if overlaps
						if t.CDSStart > 0 && t.CDSEnd > 0 {
							if e.End >= t.CDSStart && e.Start <= t.CDSEnd {
								t.Exons[i].CDSStart = max(e.Start, t.CDSStart)
								t.Exons[i].CDSEnd = min(e.End, t.CDSEnd)
							}
						}
					}
					// Sort exons by genomic position (start coordinate)
					// This is required for correct CDS position calculation
					sort.Slice(t.Exons, func(i, j int) bool {
						return t.Exons[i].Start < t.Exons[j].Start
					})
				}
			}
		}
	}

	// Fetch CDS sequence
	cdsURL := fmt.Sprintf("%s/sequence/id/%s?type=cds;content-type=application/json",
		l.baseURL, t.ID)

	if resp, err := l.httpClient.Get(cdsURL); err == nil {
		defer resp.Body.Close()
		if resp.StatusCode == http.StatusOK {
			var seqResp struct {
				Seq string `json:"seq"`
			}
			if json.NewDecoder(resp.Body).Decode(&seqResp) == nil {
				t.CDSSequence = seqResp.Seq
			}
		}
	}

	// Fetch protein sequence
	protURL := fmt.Sprintf("%s/sequence/id/%s?type=protein;content-type=application/json",
		l.baseURL, t.ID)

	if resp, err := l.httpClient.Get(protURL); err == nil {
		defer resp.Body.Close()
		if resp.StatusCode == http.StatusOK {
			var seqResp struct {
				Seq string `json:"seq"`
			}
			if json.NewDecoder(resp.Body).Decode(&seqResp) == nil {
				t.ProteinSequence = seqResp.Seq
			}
		}
	}
}
