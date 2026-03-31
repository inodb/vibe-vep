package server

import (
	"encoding/json"
	"net/http"
	"strings"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/datasource/pfam"
)

// ensemblTranscriptResponse is the JSON response for canonical transcript / transcript lookup endpoints.
type ensemblTranscriptResponse struct {
	TranscriptID  string              `json:"transcriptId"`
	GeneID        string              `json:"geneId"`
	ProteinID     string              `json:"proteinId,omitempty"`
	HugoSymbols   []string            `json:"hugoSymbols"`
	ProteinLength int                 `json:"proteinLength"`
	PfamDomains   []pfamDomainJSON    `json:"pfamDomains"`
	Exons         []exonJSON          `json:"exons"`
	UTRs          []utrJSON           `json:"utrs"`
	UniprotID     string              `json:"uniprotId,omitempty"`
	EntrezGeneID  string              `json:"entrezGeneId,omitempty"`
}

type pfamDomainJSON struct {
	PfamDomainID    string `json:"pfamDomainId"`
	PfamDomainStart int    `json:"pfamDomainStart"`
	PfamDomainEnd   int    `json:"pfamDomainEnd"`
}

type exonJSON struct {
	ExonID   string `json:"exonId"`
	ExonStart int64  `json:"exonStart"`
	ExonEnd   int64  `json:"exonEnd"`
	Rank      int    `json:"rank"`
	Strand    int8   `json:"strand"`
	Version   int    `json:"version"`
}

type utrJSON struct {
	Type   string `json:"type"`
	Start  int64  `json:"start"`
	End    int64  `json:"end"`
	Strand int8   `json:"strand"`
}

// pfamDomainResponse is the JSON response for PFAM domain lookup endpoints.
type pfamDomainResponse struct {
	PfamAccession string `json:"pfamAccession"`
	Name          string `json:"name"`
	Description   string `json:"description"`
}

// handleEnsemblCanonicalByHugo handles GET /genome-nexus/{assembly}/ensembl/canonical-transcript/hgnc/{hugoSymbol}
func (s *Server) handleEnsemblCanonicalByHugo(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	hugoSymbol := r.PathValue("hugoSymbol")
	if hugoSymbol == "" {
		writeError(w, http.StatusBadRequest, "missing hugo symbol")
		return
	}

	isoformSource := r.URL.Query().Get("isoformOverrideSource")

	transcripts := ctx.cache.FindTranscriptsByGene(hugoSymbol)
	if len(transcripts) == 0 {
		writeError(w, http.StatusNotFound, "no transcripts found for gene "+hugoSymbol)
		return
	}

	tx := pickCanonical(transcripts, isoformSource)
	if tx == nil {
		writeError(w, http.StatusNotFound, "no canonical transcript found for gene "+hugoSymbol)
		return
	}

	resp := buildTranscriptResponse(tx, ctx.pfam)
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(resp)
}

// handleEnsemblTranscriptByID handles GET /genome-nexus/{assembly}/ensembl/transcript/{transcriptId}
func (s *Server) handleEnsemblTranscriptByID(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	txID := r.PathValue("transcriptId")
	if txID == "" {
		writeError(w, http.StatusBadRequest, "missing transcript ID")
		return
	}

	tx := ctx.cache.GetTranscript(txID)
	if tx == nil {
		// Try without version suffix
		if idx := strings.IndexByte(txID, '.'); idx >= 0 {
			tx = ctx.cache.GetTranscript(txID[:idx])
		}
	}
	if tx == nil {
		writeError(w, http.StatusNotFound, "transcript not found: "+txID)
		return
	}

	resp := buildTranscriptResponse(tx, ctx.pfam)
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(resp)
}

// handlePfamDomainPost handles POST /genome-nexus/{assembly}/pfam/domain
// Body: JSON array of PFAM accession strings.
func (s *Server) handlePfamDomainPost(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	if ctx.pfam == nil {
		writeError(w, http.StatusServiceUnavailable, "PFAM data not loaded")
		return
	}

	var accessions []string
	if err := json.NewDecoder(r.Body).Decode(&accessions); err != nil {
		writeError(w, http.StatusBadRequest, "invalid JSON body: "+err.Error())
		return
	}

	results := make([]pfamDomainResponse, 0, len(accessions))
	for _, acc := range accessions {
		if d, ok := ctx.pfam.LookupDomain(acc); ok {
			results = append(results, pfamDomainResponse{
				PfamAccession: d.Accession,
				Name:          d.Name,
				Description:   d.Description,
			})
		}
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(results)
}

// handlePfamDomainGet handles GET /genome-nexus/{assembly}/pfam/domain/{pfamAccession}
func (s *Server) handlePfamDomainGet(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	if ctx.pfam == nil {
		writeError(w, http.StatusServiceUnavailable, "PFAM data not loaded")
		return
	}

	acc := r.PathValue("pfamAccession")
	d, ok := ctx.pfam.LookupDomain(acc)
	if !ok {
		writeError(w, http.StatusNotFound, "PFAM domain not found: "+acc)
		return
	}

	resp := pfamDomainResponse{
		PfamAccession: d.Accession,
		Name:          d.Name,
		Description:   d.Description,
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(resp)
}

// pickCanonical selects the canonical transcript based on isoform override source.
func pickCanonical(transcripts []*cache.Transcript, isoformSource string) *cache.Transcript {
	useMSK := strings.EqualFold(isoformSource, "mskcc")

	// First pass: look for the matching canonical flag.
	for _, tx := range transcripts {
		if useMSK && tx.IsCanonicalMSK {
			return tx
		}
		if !useMSK && tx.IsCanonicalEnsembl {
			return tx
		}
	}

	// Fallback: try the other canonical flag.
	for _, tx := range transcripts {
		if tx.IsCanonicalEnsembl || tx.IsCanonicalMSK {
			return tx
		}
	}

	// Last resort: pick the first protein-coding transcript.
	for _, tx := range transcripts {
		if tx.IsProteinCoding() {
			return tx
		}
	}

	// Nothing suitable found.
	if len(transcripts) > 0 {
		return transcripts[0]
	}
	return nil
}

// buildTranscriptResponse constructs the JSON response for a transcript.
func buildTranscriptResponse(tx *cache.Transcript, pfamStore *pfam.Store) ensemblTranscriptResponse {
	resp := ensemblTranscriptResponse{
		TranscriptID:  tx.ID,
		GeneID:        tx.GeneID,
		ProteinID:     tx.ProteinID,
		HugoSymbols:   []string{tx.GeneName},
		ProteinLength: tx.ProteinLength,
		EntrezGeneID:  tx.EntrezGeneID,
	}

	// PFAM domains
	if pfamStore != nil {
		ranges := pfamStore.LookupTranscript(tx.ID)
		resp.PfamDomains = make([]pfamDomainJSON, 0, len(ranges))
		for _, dr := range ranges {
			resp.PfamDomains = append(resp.PfamDomains, pfamDomainJSON{
				PfamDomainID:    dr.PfamDomainID,
				PfamDomainStart: dr.Start,
				PfamDomainEnd:   dr.End,
			})
		}
	}
	if resp.PfamDomains == nil {
		resp.PfamDomains = []pfamDomainJSON{}
	}

	// Exons
	resp.Exons = make([]exonJSON, 0, len(tx.Exons))
	for _, e := range tx.Exons {
		resp.Exons = append(resp.Exons, exonJSON{
			ExonID:    "", // Exon IDs not available in current GENCODE parsing
			ExonStart: e.Start,
			ExonEnd:   e.End,
			Rank:      e.Number,
			Strand:    tx.Strand,
			Version:   1,
		})
	}

	// UTRs: derive from CDS boundaries and transcript boundaries.
	resp.UTRs = buildUTRs(tx)

	return resp
}

// buildUTRs derives UTR regions from the transcript CDS and exon boundaries.
func buildUTRs(tx *cache.Transcript) []utrJSON {
	if !tx.IsProteinCoding() || len(tx.Exons) == 0 {
		return []utrJSON{}
	}

	var utrs []utrJSON

	// 5' UTR: region before CDS start (in transcript orientation).
	// 3' UTR: region after CDS end (in transcript orientation).
	if tx.Strand == 1 {
		// Forward strand: 5'UTR is tx.Start..CDSStart-1, 3'UTR is CDSEnd+1..tx.End
		if tx.CDSStart > tx.Start {
			utrs = append(utrs, utrJSON{
				Type:   "five_prime_UTR",
				Start:  tx.Start,
				End:    tx.CDSStart - 1,
				Strand: tx.Strand,
			})
		}
		if tx.CDSEnd < tx.End {
			utrs = append(utrs, utrJSON{
				Type:   "three_prime_UTR",
				Start:  tx.CDSEnd + 1,
				End:    tx.End,
				Strand: tx.Strand,
			})
		}
	} else {
		// Reverse strand: 5'UTR is CDSEnd+1..tx.End, 3'UTR is tx.Start..CDSStart-1
		if tx.CDSEnd < tx.End {
			utrs = append(utrs, utrJSON{
				Type:   "five_prime_UTR",
				Start:  tx.CDSEnd + 1,
				End:    tx.End,
				Strand: tx.Strand,
			})
		}
		if tx.CDSStart > tx.Start {
			utrs = append(utrs, utrJSON{
				Type:   "three_prime_UTR",
				Start:  tx.Start,
				End:    tx.CDSStart - 1,
				Strand: tx.Strand,
			})
		}
	}

	return utrs
}
