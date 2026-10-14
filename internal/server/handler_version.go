package server

import (
	"encoding/json"
	"net/http"
	"sort"
)

// A genome-nexus client asks /version before it asks anything else, and reads
// back what produced the annotation it is about to receive. The shape below is
// genome-nexus's, so a client reading it needs no special case for this
// server; the contents are this server's own, naming what it actually has
// loaded rather than what genome-nexus would have had.

type versionInfo struct {
	Version string `json:"version"`
	Static  bool   `json:"static"`
}

type sourceInfo struct {
	ID      string `json:"id"`
	Name    string `json:"name"`
	Version string `json:"version"`
	// Type is "mirrored" for data held locally and "external" for data
	// fetched from a service as it is needed.
	Type        string `json:"type"`
	Description string `json:"description,omitempty"`
	URL         string `json:"url,omitempty"`
}

func (s *Server) handleVersion(w http.ResponseWriter, r *http.Request) {
	s.mu.RLock()
	defer s.mu.RUnlock()

	// The annotator is the server, so its version answers for both.
	response := map[string]any{
		"genomeNexus": map[string]any{
			"server": versionInfo{Version: s.version, Static: true},
		},
		"vep": map[string]any{
			"server": versionInfo{Version: s.version, Static: true},
		},
	}

	assemblies := make([]string, 0, len(s.assemblies))
	for _, ctx := range s.assemblies {
		assemblies = append(assemblies, ctx.assembly)
	}
	sort.Strings(assemblies)

	sources := []sourceInfo{{
		ID:          "vibe-vep",
		Name:        "vibe-vep",
		Version:     s.version,
		Type:        "mirrored",
		Description: "Variant effect prediction from GENCODE transcript sequences.",
		URL:         "https://github.com/inodb/vibe-vep",
	}}
	for _, assembly := range assemblies {
		sources = append(sources, sourceInfo{
			ID:          "gencode",
			Name:        "GENCODE",
			Version:     assembly,
			Type:        "mirrored",
			Description: "Gene annotation the consequences are predicted against.",
			URL:         "https://www.gencodegenes.org/",
		})
	}
	if len(s.overrideTranscripts) > 0 {
		sources = append(sources, sourceInfo{
			ID:          "hgnc",
			Name:        "HGNC",
			Version:     "ensembl_biomart_canonical_transcripts_per_hgnc",
			Type:        "mirrored",
			Description: "Approved gene symbols and the canonical transcript each isoform override names.",
			URL:         "https://www.genenames.org/",
		})
	}
	if len(s.cancerGenes) > 0 {
		sources = append(sources, sourceInfo{
			ID:          "oncokb",
			Name:        "OncoKB",
			Version:     "cancerGeneList.tsv",
			Type:        "mirrored",
			Description: "Cancer gene list, used to decide which gene a variant is reported against.",
			URL:         "https://www.oncokb.org/",
		})
	}
	if s.revue != nil {
		sources = append(sources, sourceInfo{
			ID:          "revue",
			Name:        "reVUE",
			Version:     "VUEs.json",
			Type:        "mirrored",
			Description: "Published protein effects for variants whose real effect the sequence does not show.",
			URL:         "https://www.cancerrevue.org/",
		})
	}
	response["annotationSourcesInfo"] = sources

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(response)
}
