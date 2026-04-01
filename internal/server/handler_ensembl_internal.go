package server

import (
	"encoding/json"
	"net/http"
	"strings"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/datasource/pfam"
	"github.com/inodb/vibe-vep/internal/datasource/ptm"
	"github.com/inodb/vibe-vep/internal/datasource/uniprot"
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

	resp := buildTranscriptResponse(tx, ctx.pfam, ctx.uniprot)
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
		// Input is likely unversioned (e.g. "ENST00000357654") but cache
		// stores versioned IDs (e.g. "ENST00000357654.9"). Try prefix match.
		tx = ctx.cache.GetTranscriptByPrefix(txID)
	}
	if tx == nil {
		writeError(w, http.StatusNotFound, "transcript not found: "+txID)
		return
	}

	resp := buildTranscriptResponse(tx, ctx.pfam, ctx.uniprot)
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
func buildTranscriptResponse(tx *cache.Transcript, pfamStore *pfam.Store, uniprotStore ...*uniprot.Store) ensemblTranscriptResponse {
	resp := ensemblTranscriptResponse{
		TranscriptID:  stripTxVersion(tx.ID),
		GeneID:        tx.GeneID,
		ProteinID:     tx.ProteinID,
		HugoSymbols:   []string{tx.GeneName},
		ProteinLength: tx.ProteinLength,
		EntrezGeneID:  tx.EntrezGeneID,
	}

	// UniProt ID
	if len(uniprotStore) > 0 && uniprotStore[0] != nil {
		resp.UniprotID = uniprotStore[0].LookupByTranscript(tx.ID)
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

// buildUTRs derives per-exon UTR regions from the transcript CDS and exon boundaries.
// UTRs are the non-coding exonic portions before/after the CDS.
func buildUTRs(tx *cache.Transcript) []utrJSON {
	if !tx.IsProteinCoding() || len(tx.Exons) == 0 {
		return []utrJSON{}
	}

	var utrs []utrJSON
	strand := tx.Strand

	for _, exon := range tx.Exons {
		// Determine if this exon has a UTR portion.
		// UTR is the part of the exon outside the CDS.
		if exon.End < tx.CDSStart || exon.Start > tx.CDSEnd {
			// Entire exon is outside CDS — it's a UTR exon.
			utrType := "three_prime_UTR"
			if strand == 1 && exon.End < tx.CDSStart {
				utrType = "five_prime_UTR"
			} else if strand == -1 && exon.Start > tx.CDSEnd {
				utrType = "five_prime_UTR"
			}
			utrs = append(utrs, utrJSON{
				Type:   utrType,
				Start:  exon.Start,
				End:    exon.End,
				Strand: strand,
			})
		} else {
			// Exon partially overlaps CDS — split into UTR + CDS portions.
			if exon.Start < tx.CDSStart {
				// Left portion is UTR.
				utrType := "five_prime_UTR"
				if strand == -1 {
					utrType = "three_prime_UTR"
				}
				utrs = append(utrs, utrJSON{
					Type:   utrType,
					Start:  exon.Start,
					End:    tx.CDSStart - 1,
					Strand: strand,
				})
			}
			if exon.End > tx.CDSEnd {
				// Right portion is UTR.
				utrType := "three_prime_UTR"
				if strand == -1 {
					utrType = "five_prime_UTR"
				}
				utrs = append(utrs, utrJSON{
					Type:   utrType,
					Start:  tx.CDSEnd + 1,
					End:    exon.End,
					Strand: strand,
				})
			}
		}
	}

	return utrs
}

// handleEnsemblTranscriptPost handles POST /genome-nexus/{assembly}/ensembl/transcript
// Body: {"transcriptIds":["ENST..."], "geneIds":[], "hugoSymbols":[], "proteinIds":[]}
// Returns an array of transcript responses matching any of the filter criteria.
func (s *Server) handleEnsemblTranscriptPost(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	var filter struct {
		TranscriptIDs []string `json:"transcriptIds"`
		GeneIDs       []string `json:"geneIds"`
		HugoSymbols   []string `json:"hugoSymbols"`
		ProteinIDs    []string `json:"proteinIds"`
	}
	if err := json.NewDecoder(r.Body).Decode(&filter); err != nil {
		writeError(w, http.StatusBadRequest, "invalid JSON: "+err.Error())
		return
	}

	pfamStore := ctx.pfam
	uniprotStore := ctx.uniprot
	var results []ensemblTranscriptResponse

	// Lookup by transcript IDs
	for _, txID := range filter.TranscriptIDs {
		tx := ctx.cache.GetTranscript(txID)
		if tx == nil {
			tx = ctx.cache.GetTranscriptByPrefix(txID)
		}
		if tx != nil {
			results = append(results, buildTranscriptResponse(tx, pfamStore, uniprotStore))
		}
	}

	// Lookup by hugo symbols (all transcripts for the gene).
	for _, symbol := range filter.HugoSymbols {
		for _, chrom := range ctx.cache.Chromosomes() {
			for _, t := range ctx.cache.FindTranscriptsByChrom(chrom) {
				if t.GeneName == symbol {
					results = append(results, buildTranscriptResponse(t, pfamStore, uniprotStore))
				}
			}
		}
	}

	w.Header().Set("Content-Type", "application/json")
	writeJSON(w, http.StatusOK, results)
}

// ensemblGeneResponse is the JSON response for canonical gene lookups.
type ensemblGeneResponse struct {
	GeneID       string `json:"geneId"`
	HugoSymbol   string `json:"hugoSymbol"`
	EntrezGeneID string `json:"entrezGeneId,omitempty"`
	UniprotID    string `json:"uniprotId,omitempty"`
}

// handleEnsemblCanonicalGeneByHugo handles GET /genome-nexus/{assembly}/ensembl/canonical-gene/hgnc/{hugoSymbol}
func (s *Server) handleEnsemblCanonicalGeneByHugo(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	symbol := r.PathValue("hugoSymbol")
	tx := findCanonicalByGene(ctx.cache, symbol)
	if tx == nil {
		writeError(w, http.StatusNotFound, "gene not found: "+symbol)
		return
	}

	uniprotID := ""
	if ctx.uniprot != nil {
		uniprotID = ctx.uniprot.LookupByTranscript(tx.ID)
	}
	writeJSON(w, http.StatusOK, ensemblGeneResponse{
		GeneID:       tx.GeneID,
		HugoSymbol:   tx.GeneName,
		EntrezGeneID: tx.EntrezGeneID,
		UniprotID:    uniprotID,
	})
}

// handleEnsemblCanonicalGeneByEntrez handles GET /genome-nexus/{assembly}/ensembl/canonical-gene/entrez/{entrezGeneId}
func (s *Server) handleEnsemblCanonicalGeneByEntrez(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	entrezID := r.PathValue("entrezGeneId")
	var found *cache.Transcript
	for _, chrom := range ctx.cache.Chromosomes() {
		for _, t := range ctx.cache.FindTranscriptsByChrom(chrom) {
			if t.EntrezGeneID == entrezID && t.IsCanonicalMSK {
				found = t
				break
			}
		}
		if found != nil {
			break
		}
	}
	if found == nil {
		writeError(w, http.StatusNotFound, "gene not found for entrez ID: "+entrezID)
		return
	}

	entrezUniprotID := ""
	if ctx.uniprot != nil {
		entrezUniprotID = ctx.uniprot.LookupByTranscript(found.ID)
	}
	writeJSON(w, http.StatusOK, ensemblGeneResponse{
		GeneID:       found.GeneID,
		HugoSymbol:   found.GeneName,
		EntrezGeneID: found.EntrezGeneID,
		UniprotID:    entrezUniprotID,
	})
}

// handleCancerHotspotsTranscript handles GET /genome-nexus/{assembly}/cancer_hotspots/transcript/{transcriptId}
// Returns hotspot positions on the protein for the given transcript.
func (s *Server) handleCancerHotspotsTranscript(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	// TODO: look up hotspots from the hotspot store by transcript ID
	// For now return empty array (matches GN behavior for most transcripts).
	w.Header().Set("Content-Type", "application/json")
	w.Write([]byte("[]\n"))
}

// findCanonicalByGene finds the MSK canonical transcript for a gene symbol.
func findCanonicalByGene(c *cache.Cache, geneName string) *cache.Transcript {
	for _, chrom := range c.Chromosomes() {
		for _, t := range c.FindTranscriptsByChrom(chrom) {
			if t.GeneName == geneName && t.IsCanonicalMSK {
				return t
			}
		}
	}
	return nil
}

// handlePtmExperimentalGet handles GET /genome-nexus/{assembly}/ptm/experimental?ensemblTranscriptId={id}
func (s *Server) handlePtmExperimentalGet(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	txID := r.URL.Query().Get("ensemblTranscriptId")
	if txID == "" {
		w.Header().Set("Content-Type", "application/json")
		w.Write([]byte("[]\n"))
		return
	}

	ptms := lookupPtms(ctx.ptm, txID)
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(ptms)
}

// handlePtmExperimentalPost handles POST /genome-nexus/{assembly}/ptm/experimental
// Body: JSON object with ensemblTranscriptId field, or array of transcript IDs.
func (s *Server) handlePtmExperimentalPost(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	var body struct {
		EnsemblTranscriptIds []string `json:"ensemblTranscriptIds"`
		EnsemblTranscriptId  string   `json:"ensemblTranscriptId"`
	}
	if err := json.NewDecoder(r.Body).Decode(&body); err != nil {
		writeError(w, http.StatusBadRequest, "invalid JSON: "+err.Error())
		return
	}

	var allPtms []ptm.PTM
	if body.EnsemblTranscriptId != "" {
		allPtms = append(allPtms, lookupPtms(ctx.ptm, body.EnsemblTranscriptId)...)
	}
	for _, txID := range body.EnsemblTranscriptIds {
		allPtms = append(allPtms, lookupPtms(ctx.ptm, txID)...)
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(allPtms)
}

// lookupPtms returns PTMs for a transcript, falling back to an empty slice.
func lookupPtms(store *ptm.Store, transcriptID string) []ptm.PTM {
	if store == nil {
		return []ptm.PTM{}
	}
	result := store.LookupByTranscript(transcriptID)
	if result == nil {
		return []ptm.PTM{}
	}
	return result
}

// stripTxVersion removes the version suffix from a transcript ID (e.g. "ENST00000357654.9" → "ENST00000357654").
func stripTxVersion(id string) string {
	if i := strings.IndexByte(id, '.'); i >= 0 {
		return id[:i]
	}
	return id
}
