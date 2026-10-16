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
	TranscriptID  string           `json:"transcriptId"`
	GeneID        string           `json:"geneId"`
	ProteinID     string           `json:"proteinId,omitempty"`
	HugoSymbols   []string         `json:"hugoSymbols"`
	ProteinLength int              `json:"proteinLength"`
	PfamDomains   []pfamDomainJSON `json:"pfamDomains"`
	Exons         []exonJSON       `json:"exons"`
	UTRs          []utrJSON        `json:"utrs"`
	UniprotID     string           `json:"uniprotId,omitempty"`
	EntrezGeneID  string           `json:"entrezGeneId,omitempty"`
}

type pfamDomainJSON struct {
	PfamDomainID    string `json:"pfamDomainId"`
	PfamDomainStart int    `json:"pfamDomainStart"`
	PfamDomainEnd   int    `json:"pfamDomainEnd"`
}

type exonJSON struct {
	ExonID    string `json:"exonId"`
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

	tx := pickCanonical(transcripts, isoformSource, s.OverrideTranscripts(isoformSource))
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

// pickCanonical selects the canonical transcript for an isoform override
// source.
//
// The override set is consulted first and is keyed by transcript id across all
// genes, not per gene. That matters wherever the GENCODE release still carries
// an older gene symbol than the override table -- GENCODE calls NSD2 "WHSC1",
// so a per-gene lookup finds nothing and the gene silently loses its override.
func pickCanonical(transcripts []*cache.Transcript, isoformSource string, overrides map[string]bool) *cache.Transcript {
	for _, tx := range transcripts {
		if overrides[stripTxVersion(tx.ID)] {
			return tx
		}
	}

	// No override knows any of these transcripts: fall back to the flags the
	// transcripts carry.
	useMSK := strings.EqualFold(cache.NormalizeIsoformSource(isoformSource), cache.IsoformSourceMSKCC)
	for _, tx := range transcripts {
		if useMSK && tx.IsCanonicalMSK {
			return tx
		}
		if !useMSK && tx.IsCanonicalEnsembl {
			return tx
		}
	}

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

	// The four lists are alternatives, not a union: genome-nexus takes the
	// first non-empty one and ignores the rest, in this order. Sending both
	// transcriptIds and hugoSymbols returns the transcripts, not both sets.
	results := make([]ensemblTranscriptResponse, 0)
	switch {
	case len(filter.TranscriptIDs) > 0:
		for _, txID := range filter.TranscriptIDs {
			tx := ctx.cache.GetTranscript(txID)
			if tx == nil {
				tx = ctx.cache.GetTranscriptByPrefix(txID)
			}
			if tx != nil {
				results = append(results, buildTranscriptResponse(tx, ctx.pfam, ctx.uniprot))
			}
		}
	case len(filter.GeneIDs) > 0:
		for _, geneID := range filter.GeneIDs {
			results = appendTranscripts(results, ctx, ctx.cache.FindTranscriptsByGeneID(geneID))
		}
	case len(filter.ProteinIDs) > 0:
		for _, proteinID := range filter.ProteinIDs {
			results = appendTranscripts(results, ctx, ctx.cache.FindTranscriptsByProteinID(proteinID))
		}
	case len(filter.HugoSymbols) > 0:
		for _, symbol := range filter.HugoSymbols {
			results = appendTranscripts(results, ctx, ctx.cache.FindTranscriptsByGene(symbol))
		}
	}

	writeJSON(w, http.StatusOK, results)
}

// handleEnsemblTranscriptGet handles GET {prefix}/ensembl/transcript, the
// single-valued form of the filter above.
//
// Here the parameters are conjunctive rather than alternative: a gene and a
// protein together mean transcripts matching both. With no parameters at all
// the answer is an empty list -- genome-nexus's own documentation claims it
// returns every transcript, but its implementation does not, and the
// implementation is what clients were built against.
func (s *Server) handleEnsemblTranscriptGet(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	q := r.URL.Query()
	geneID := q.Get("geneId")
	proteinID := q.Get("proteinId")
	hugoSymbol := q.Get("hugoSymbol")

	results := make([]ensemblTranscriptResponse, 0)
	if geneID == "" && proteinID == "" && hugoSymbol == "" {
		writeJSON(w, http.StatusOK, results)
		return
	}

	// Start from whichever identifier is most selective, then narrow.
	var candidates []*cache.Transcript
	switch {
	case geneID != "":
		candidates = ctx.cache.FindTranscriptsByGeneID(geneID)
	case proteinID != "":
		candidates = ctx.cache.FindTranscriptsByProteinID(proteinID)
	default:
		candidates = ctx.cache.FindTranscriptsByGene(hugoSymbol)
	}

	for _, tx := range candidates {
		if proteinID != "" && stripTxVersion(tx.ProteinID) != stripTxVersion(proteinID) {
			continue
		}
		if hugoSymbol != "" && tx.GeneName != hugoSymbol {
			continue
		}
		results = append(results, buildTranscriptResponse(tx, ctx.pfam, ctx.uniprot))
	}

	writeJSON(w, http.StatusOK, results)
}

// appendTranscripts renders a set of transcripts into the response slice.
func appendTranscripts(dst []ensemblTranscriptResponse, ctx *assemblyContext, txs []*cache.Transcript) []ensemblTranscriptResponse {
	for _, tx := range txs {
		dst = append(dst, buildTranscriptResponse(tx, ctx.pfam, ctx.uniprot))
	}
	return dst
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
	source := r.URL.Query().Get("isoformOverrideSource")
	tx := pickCanonical(ctx.cache.FindTranscriptsByGene(symbol), source, s.OverrideTranscripts(source))
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
	source := r.URL.Query().Get("isoformOverrideSource")
	found := pickCanonical(ctx.cache.FindTranscriptsByEntrezID(entrezID), source, s.OverrideTranscripts(source))
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

// decodeStringList reads the plain JSON array of identifiers the batch
// endpoints take as a body.
func decodeStringList(w http.ResponseWriter, r *http.Request) ([]string, bool) {
	var ids []string
	if err := json.NewDecoder(r.Body).Decode(&ids); err != nil {
		writeError(w, http.StatusBadRequest, "invalid JSON body: "+err.Error())
		return nil, false
	}
	return ids, true
}

// geneResponseFor renders a gene from the transcript chosen to represent it.
func geneResponseFor(ctx *assemblyContext, tx *cache.Transcript) ensemblGeneResponse {
	uniprotID := ""
	if ctx.uniprot != nil {
		uniprotID = ctx.uniprot.LookupByTranscript(tx.ID)
	}
	return ensemblGeneResponse{
		GeneID:       tx.GeneID,
		HugoSymbol:   tx.GeneName,
		EntrezGeneID: tx.EntrezGeneID,
		UniprotID:    uniprotID,
	}
}

// handleEnsemblCanonicalByHugoPost handles POST
// {prefix}/ensembl/canonical-transcript/hgnc, the batch form of the lookup
// above. A symbol with no canonical transcript is left out of the answer
// rather than failing the request.
func (s *Server) handleEnsemblCanonicalByHugoPost(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}
	symbols, ok := decodeStringList(w, r)
	if !ok {
		return
	}

	isoformSource := r.URL.Query().Get("isoformOverrideSource")
	overrides := s.OverrideTranscripts(isoformSource)

	results := make([]ensemblTranscriptResponse, 0, len(symbols))
	for _, symbol := range symbols {
		tx := pickCanonical(ctx.cache.FindTranscriptsByGene(symbol), isoformSource, overrides)
		if tx == nil {
			continue
		}
		results = append(results, buildTranscriptResponse(tx, ctx.pfam, ctx.uniprot))
	}
	writeJSON(w, http.StatusOK, results)
}

// handleEnsemblCanonicalGeneByHugoPost handles POST
// {prefix}/ensembl/canonical-gene/hgnc.
func (s *Server) handleEnsemblCanonicalGeneByHugoPost(w http.ResponseWriter, r *http.Request) {
	s.canonicalGenesPost(w, r, func(ctx *assemblyContext, id string) []*cache.Transcript {
		return ctx.cache.FindTranscriptsByGene(id)
	})
}

// handleEnsemblCanonicalGeneByEntrezPost handles POST
// {prefix}/ensembl/canonical-gene/entrez.
func (s *Server) handleEnsemblCanonicalGeneByEntrezPost(w http.ResponseWriter, r *http.Request) {
	s.canonicalGenesPost(w, r, func(ctx *assemblyContext, id string) []*cache.Transcript {
		return ctx.cache.FindTranscriptsByEntrezID(id)
	})
}

// canonicalGenesPost answers a batch of gene lookups, differing only in how an
// identifier reaches its transcripts.
func (s *Server) canonicalGenesPost(w http.ResponseWriter, r *http.Request, find func(*assemblyContext, string) []*cache.Transcript) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}
	ids, ok := decodeStringList(w, r)
	if !ok {
		return
	}

	source := r.URL.Query().Get("isoformOverrideSource")
	overrides := s.OverrideTranscripts(source)

	results := make([]ensemblGeneResponse, 0, len(ids))
	for _, id := range ids {
		tx := pickCanonical(find(ctx, id), source, overrides)
		if tx == nil {
			continue
		}
		results = append(results, geneResponseFor(ctx, tx))
	}
	writeJSON(w, http.StatusOK, results)
}
