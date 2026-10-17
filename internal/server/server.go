// Package server provides an HTTP server for variant annotation.
package server

import (
	"encoding/json"
	"net/http"
	"strings"
	"sync"

	"go.uber.org/zap"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/datasource/ensembltranscript"
	"github.com/inodb/vibe-vep/internal/datasource/hotspots"
	"github.com/inodb/vibe-vep/internal/datasource/mutationassessor"
	"github.com/inodb/vibe-vep/internal/datasource/myvariantinfo"
	"github.com/inodb/vibe-vep/internal/datasource/pfam"
	"github.com/inodb/vibe-vep/internal/datasource/ptm"
	"github.com/inodb/vibe-vep/internal/datasource/revue"
	"github.com/inodb/vibe-vep/internal/datasource/uniprot"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// assemblyContext holds the annotator and sources for one assembly.
type assemblyContext struct {
	annotator *annotate.Annotator
	sources   []annotate.AnnotationSource
	cache     *cache.Cache
	pfam      *pfam.Store
	ptm       *ptm.Store
	uniprot   *uniprot.Store
	hotspots  *hotspots.Store
	mutAssess *mutationassessor.Store
	// ensemblTx holds the Ensembl transcript records the /ensembl endpoints
	// report. They cannot be derived from the GENCODE cache, so without it
	// those endpoints answer from the cache and omit what it cannot supply.
	ensemblTx *ensembltranscript.Store
	assembly  string // normalized: "GRCh38"
}

// Server is the HTTP annotation server.
type Server struct {
	mu              sync.RWMutex
	assemblies      map[string]*assemblyContext // keyed by lowercase: "grch38"
	logger          *zap.Logger
	version         string
	myVariantClient *myvariantinfo.Client

	// cancerGenes and curatedGenes are the OncoKB cancer gene list and the
	// curated subset of it. genome-nexus uses them to decide which gene a
	// variant falling in two of them is reported against, so without them
	// those variants can land on the wrong gene.
	cancerGenes  map[string]bool
	curatedGenes map[string]bool

	// overrideTranscripts holds, per isoform override source, the transcript
	// ids that source names canonical. Keyed by lowercase source name.
	overrideTranscripts map[string]map[string]bool

	// officialGeneSymbols maps a retired gene symbol to the approved one, so
	// a gene the GENCODE release still calls WHSC1 is reported as NSD2.
	officialGeneSymbols map[string]string

	// geneMetadata holds what the biomart HGNC table says about each gene
	// beyond its canonical transcript: its HGNC id, which VEP reports on
	// every transcript consequence, and the other names it answers to.
	geneMetadata map[string]cache.GeneMetadata

	// revue holds the published protein effects that replace the predicted
	// one for variants whose real effect the sequence does not show.
	revue *revue.Store
}

// New creates a new Server.
func New(logger *zap.Logger, version string) *Server {
	return &Server{
		assemblies:      make(map[string]*assemblyContext),
		logger:          logger,
		version:         version,
		myVariantClient: myvariantinfo.NewClient(),
	}
}

// SetHotspotsStore records the cancer hotspots table for an assembly, which
// the /cancer_hotspots endpoints answer from directly rather than through a
// variant annotation.
func (s *Server) SetHotspotsStore(assembly string, store *hotspots.Store) {
	s.mu.Lock()
	defer s.mu.Unlock()
	if ctx, ok := s.assemblies[strings.ToLower(assembly)]; ok {
		ctx.hotspots = store
	}
}

// SetEnsemblTranscriptStore records the Ensembl transcript records for an
// assembly, which the /ensembl endpoints report.
func (s *Server) SetEnsemblTranscriptStore(assembly string, store *ensembltranscript.Store) {
	s.mu.Lock()
	defer s.mu.Unlock()
	if ctx, ok := s.assemblies[strings.ToLower(assembly)]; ok {
		ctx.ensemblTx = store
	}
}

// SetMutationAssessorStore records the Mutation Assessor predictions for an
// assembly.
func (s *Server) SetMutationAssessorStore(assembly string, store *mutationassessor.Store) {
	s.mu.Lock()
	defer s.mu.Unlock()
	if ctx, ok := s.assemblies[strings.ToLower(assembly)]; ok {
		ctx.mutAssess = store
	}
}

// SetCancerGenes records the OncoKB cancer gene list and its curated subset,
// used to decide which of two overlapping genes a variant is reported against.
func (s *Server) SetCancerGenes(cancer, curated map[string]bool) {
	s.mu.Lock()
	defer s.mu.Unlock()
	s.cancerGenes = cancer
	s.curatedGenes = curated
}

// CancerGenes returns the OncoKB cancer gene list and its curated subset.
func (s *Server) CancerGenes() (cancer, curated map[string]bool) {
	s.mu.RLock()
	defer s.mu.RUnlock()
	return s.cancerGenes, s.curatedGenes
}

// SetOverrideTranscripts records the transcript ids an isoform override
// source names canonical.
func (s *Server) SetOverrideTranscripts(source string, transcripts map[string]bool) {
	s.mu.Lock()
	defer s.mu.Unlock()
	if s.overrideTranscripts == nil {
		s.overrideTranscripts = make(map[string]map[string]bool)
	}
	s.overrideTranscripts[cache.NormalizeIsoformSource(source)] = transcripts
}

// OverrideTranscripts returns the transcript ids the named source makes
// canonical. An unnamed source, and anything not recognised, means MSKCC,
// which is what the annotation pipeline asks for and what genome-nexus
// falls back to.
func (s *Server) OverrideTranscripts(source string) map[string]bool {
	s.mu.RLock()
	defer s.mu.RUnlock()
	return s.overrideTranscripts[cache.NormalizeIsoformSource(source)]
}

// SetOfficialGeneSymbols records the retired-to-approved gene symbol map.
func (s *Server) SetOfficialGeneSymbols(symbols map[string]string) {
	s.mu.Lock()
	defer s.mu.Unlock()
	s.officialGeneSymbols = symbols
}

// OfficialGeneSymbols returns the retired-to-approved gene symbol map.
func (s *Server) OfficialGeneSymbols() map[string]string {
	s.mu.RLock()
	defer s.mu.RUnlock()
	return s.officialGeneSymbols
}

// SetGeneMetadata records the biomart HGNC table's gene identifiers and
// alternative names, keyed by approved symbol.
func (s *Server) SetGeneMetadata(meta map[string]cache.GeneMetadata) {
	s.mu.Lock()
	defer s.mu.Unlock()
	s.geneMetadata = meta
}

// GeneMetadata returns what is known about a gene by its approved symbol,
// following the retired-symbol map first so that a GENCODE release still
// naming a gene WHSC1 finds NSD2's entry.
func (s *Server) GeneMetadata(symbol string) (cache.GeneMetadata, bool) {
	s.mu.RLock()
	defer s.mu.RUnlock()
	if meta, ok := s.geneMetadata[symbol]; ok {
		return meta, true
	}
	if official, ok := s.officialGeneSymbols[symbol]; ok {
		meta, found := s.geneMetadata[official]
		return meta, found
	}
	return cache.GeneMetadata{}, false
}

// SetRevue records the reVUE catalogue of revised protein effects.
func (s *Server) SetRevue(store *revue.Store) {
	s.mu.Lock()
	defer s.mu.Unlock()
	s.revue = store
}

// Revue returns the reVUE catalogue, or nil when none is loaded.
func (s *Server) Revue() *revue.Store {
	s.mu.RLock()
	defer s.mu.RUnlock()
	return s.revue
}

// AddAssembly registers an assembly with its annotator and sources.
func (s *Server) AddAssembly(assembly string, c *cache.Cache, ann *annotate.Annotator, sources []annotate.AnnotationSource) {
	s.mu.Lock()
	defer s.mu.Unlock()
	s.assemblies[strings.ToLower(assembly)] = &assemblyContext{
		annotator: ann,
		sources:   sources,
		cache:     c,
		assembly:  assembly,
	}
}

// SetPfamStore sets the PFAM store for the given assembly.
func (s *Server) SetPfamStore(assembly string, store *pfam.Store) {
	s.mu.Lock()
	defer s.mu.Unlock()
	if ctx, ok := s.assemblies[strings.ToLower(assembly)]; ok {
		ctx.pfam = store
	}
}

// SetPtmStore sets the PTM store for the given assembly.
func (s *Server) SetPtmStore(assembly string, store *ptm.Store) {
	s.mu.Lock()
	defer s.mu.Unlock()
	if ctx, ok := s.assemblies[strings.ToLower(assembly)]; ok {
		ctx.ptm = store
	}
}

// SetUniprotStore sets the UniProt mapping store for the given assembly.
func (s *Server) SetUniprotStore(assembly string, store *uniprot.Store) {
	s.mu.Lock()
	defer s.mu.Unlock()
	if ctx, ok := s.assemblies[strings.ToLower(assembly)]; ok {
		ctx.uniprot = store
	}
}

// getAssembly returns the assembly context for the given name (case-insensitive).
// soleAssembly returns the only loaded assembly, or nil when more than one is
// loaded. It is what lets the genome-nexus routes be served at the root: a
// genome-nexus deployment serves one assembly and its clients send no
// assembly in the path, so a single-assembly server can answer them.
func (s *Server) soleAssembly() *assemblyContext {
	s.mu.RLock()
	defer s.mu.RUnlock()
	if len(s.assemblies) != 1 {
		return nil
	}
	for _, ctx := range s.assemblies {
		return ctx
	}
	return nil
}

func (s *Server) getAssembly(name string) *assemblyContext {
	s.mu.RLock()
	defer s.mu.RUnlock()
	return s.assemblies[strings.ToLower(name)]
}

// Handler returns the HTTP handler with all routes registered.
func (s *Server) Handler() http.Handler {
	mux := http.NewServeMux()

	// Health and info endpoints. These describe the server rather than an
	// assembly, so they live at the root whatever is loaded.
	mux.HandleFunc("GET /health", s.handleHealth)
	mux.HandleFunc("GET /info", s.handleInfo)
	mux.HandleFunc("GET /version", s.handleVersion)

	// Ensembl VEP compatibility endpoints.
	mux.HandleFunc("GET /ensembl/{assembly}/vep/human/region/{region}/{allele}", s.handleEnsemblRegionGet)
	mux.HandleFunc("POST /ensembl/{assembly}/vep/human/region", s.handleEnsemblRegionPost)
	mux.HandleFunc("GET /ensembl/{assembly}/vep/human/hgvs/{notation...}", s.handleEnsemblHGVSGet)
	mux.HandleFunc("POST /ensembl/{assembly}/vep/human/hgvs", s.handleEnsemblHGVSPost)

	// Genome-nexus compatibility endpoints, under an assembly-qualified
	// prefix and, where there is only one assembly to serve, at the root as
	// well. genome-nexus serves one assembly per deployment and its clients
	// send no assembly in the path, so serving both shapes lets a client be
	// pointed at this server without being reconfigured.
	s.registerGenomeNexusRoutes(mux, "/genome-nexus/{assembly}")
	if s.soleAssembly() != nil {
		s.registerGenomeNexusRoutes(mux, "")
	}

	return corsMiddleware(mux)
}

// registerGenomeNexusRoutes registers the genome-nexus API under a path
// prefix. An empty prefix serves it at the root.
func (s *Server) registerGenomeNexusRoutes(mux *http.ServeMux, prefix string) {
	// Server-level endpoints. They are already registered at the root by
	// Handler, so only the prefixed mount adds them -- but a genome-nexus
	// client pointed at .../genome-nexus/grch37 asks for .../version there,
	// and used to get a 404.
	if prefix != "" {
		mux.HandleFunc("GET "+prefix+"/health", s.handleHealth)
		mux.HandleFunc("GET "+prefix+"/info", s.handleInfo)
		mux.HandleFunc("GET "+prefix+"/version", s.handleVersion)
	}

	mux.HandleFunc("GET "+prefix+"/annotation/genomic/{genomicLocation}", s.handleGNGenomicGet)
	mux.HandleFunc("POST "+prefix+"/annotation/genomic", s.handleGNGenomicPost)
	mux.HandleFunc("GET "+prefix+"/annotation/{variant...}", s.handleGNHGVSGet)
	mux.HandleFunc("POST "+prefix+"/annotation", s.handleGNHGVSPost)

	// Ensembl internal endpoints (transcript/PFAM lookups for frontend).
	mux.HandleFunc("GET "+prefix+"/ensembl/canonical-transcript/hgnc/{hugoSymbol}", s.handleEnsemblCanonicalByHugo)
	mux.HandleFunc("POST "+prefix+"/ensembl/canonical-transcript/hgnc", s.handleEnsemblCanonicalByHugoPost)
	mux.HandleFunc("GET "+prefix+"/ensembl/transcript/{transcriptId}", s.handleEnsemblTranscriptByID)
	mux.HandleFunc("POST "+prefix+"/ensembl/transcript", s.handleEnsemblTranscriptPost)
	mux.HandleFunc("GET "+prefix+"/ensembl/transcript", s.handleEnsemblTranscriptGet)
	mux.HandleFunc("POST "+prefix+"/pfam/domain", s.handlePfamDomainPost)
	mux.HandleFunc("GET "+prefix+"/pfam/domain/{pfamAccession}", s.handlePfamDomainGet)
	mux.HandleFunc("GET "+prefix+"/ensembl/canonical-gene/hgnc/{hugoSymbol}", s.handleEnsemblCanonicalGeneByHugo)
	mux.HandleFunc("GET "+prefix+"/ensembl/canonical-gene/entrez/{entrezGeneId}", s.handleEnsemblCanonicalGeneByEntrez)
	mux.HandleFunc("POST "+prefix+"/ensembl/canonical-gene/hgnc", s.handleEnsemblCanonicalGeneByHugoPost)
	mux.HandleFunc("POST "+prefix+"/ensembl/canonical-gene/entrez", s.handleEnsemblCanonicalGeneByEntrezPost)
	// Cancer hotspots: by transcript, by variant, by genomic location, and by
	// residue span, each in single and batch form.
	mux.HandleFunc("GET "+prefix+"/cancer_hotspots/transcript/{transcriptId}", s.handleHotspotsByTranscriptGet)
	mux.HandleFunc("POST "+prefix+"/cancer_hotspots/transcript", s.handleHotspotsByTranscriptPost)
	mux.HandleFunc("GET "+prefix+"/cancer_hotspots/hgvs/{variant...}", s.handleHotspotsByHgvsGet)
	mux.HandleFunc("POST "+prefix+"/cancer_hotspots/hgvs", s.handleHotspotsByHgvsPost)
	mux.HandleFunc("GET "+prefix+"/cancer_hotspots/genomic/{genomicLocation}", s.handleHotspotsByGenomicGet)
	mux.HandleFunc("POST "+prefix+"/cancer_hotspots/proteinLocations", s.handleHotspotsByProteinLocationsPost)

	// Mutation Assessor.
	mux.HandleFunc("GET "+prefix+"/mutation_assessor/{variant...}", s.handleMutationAssessorGet)
	mux.HandleFunc("POST "+prefix+"/mutation_assessor", s.handleMutationAssessorPost)

	// PTM endpoints.
	mux.HandleFunc("GET "+prefix+"/ptm/experimental", s.handlePtmExperimentalGet)
	mux.HandleFunc("POST "+prefix+"/ptm/experimental", s.handlePtmExperimentalPost)

	// Stub endpoints — return empty responses so a client does not hang
	// waiting for data this server does not have.
	emptyArray := func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Content-Type", "application/json")
		w.Write([]byte("[]\n"))
	}
	mux.HandleFunc("POST "+prefix+"/cancer_hotspots/genomic", s.handleHotspotsByGenomicPost)
	mux.HandleFunc("GET "+prefix+"/ensembl/xrefs", emptyArray)
}

// corsMiddleware adds CORS headers to allow cross-origin requests.
func corsMiddleware(next http.Handler) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Access-Control-Allow-Origin", "*")
		w.Header().Set("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
		w.Header().Set("Access-Control-Allow-Headers", "Content-Type")
		if r.Method == http.MethodOptions {
			w.WriteHeader(http.StatusNoContent)
			return
		}
		next.ServeHTTP(w, r)
	})
}

// handleHealth returns 200 OK for load balancer health checks.
func (s *Server) handleHealth(w http.ResponseWriter, r *http.Request) {
	w.Header().Set("Content-Type", "application/json")
	w.WriteHeader(http.StatusOK)
	json.NewEncoder(w).Encode(map[string]string{"status": "ok"})
}

// handleInfo returns server version and loaded assemblies.
func (s *Server) handleInfo(w http.ResponseWriter, r *http.Request) {
	s.mu.RLock()
	defer s.mu.RUnlock()

	type assemblyInfo struct {
		Name              string `json:"name"`
		TranscriptCount   int    `json:"transcript_count"`
		AnnotationSources int    `json:"annotation_sources"`
	}

	assemblies := make([]assemblyInfo, 0, len(s.assemblies))
	for _, ctx := range s.assemblies {
		assemblies = append(assemblies, assemblyInfo{
			Name:              ctx.assembly,
			TranscriptCount:   ctx.cache.TranscriptCount(),
			AnnotationSources: len(ctx.sources),
		})
	}

	info := map[string]any{
		"version":    s.version,
		"assemblies": assemblies,
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(info)
}

// annotateVariant runs the annotator + sources for a single variant in the given assembly context.
func (s *Server) annotateVariant(ctx *assemblyContext, v *vcf.Variant) ([]*annotate.Annotation, error) {
	anns, err := ctx.annotator.Annotate(v)
	if err != nil {
		return nil, err
	}
	for _, src := range ctx.sources {
		src.Annotate(v, anns)
	}
	return anns, nil
}

// writeJSON writes a JSON response.
func writeJSON(w http.ResponseWriter, status int, v any) {
	w.Header().Set("Content-Type", "application/json")
	w.WriteHeader(status)
	json.NewEncoder(w).Encode(v)
}

// writeError writes a JSON error response.
// writeError reports a failure the way genome-nexus does, under "message".
// Clients read that key, so reporting it under any other name is a failure
// they cannot see the reason for.
func writeError(w http.ResponseWriter, status int, msg string) {
	writeJSON(w, status, map[string]string{"message": msg})
}

// writeGNError writes a genome-nexus style error response that includes
// assembly_name and variant fields. This allows clients (like the frontend)
// to detect the genome build even when annotation fails.
// gnErrorEntry is one failed variant inside a batch response. genome-nexus
// reports a failure per item and still answers for the rest, so a single
// unparseable notation must not cost the caller the whole batch.
func gnErrorEntry(variant, assembly, msg string) json.RawMessage {
	data, err := json.Marshal(gnErrorBody(variant, assembly, msg))
	if err != nil {
		return json.RawMessage(`{"successfully_annotated":false}`)
	}
	return data
}

// gnErrorBody is the shape both the single and batch error paths report.
func gnErrorBody(variant, assembly, msg string) map[string]interface{} {
	return map[string]interface{}{
		"variant":                variant,
		"originalVariantQuery":   variant,
		"assembly_name":          assembly,
		"successfully_annotated": false,
		"errorMessage":           msg,
	}
}

func writeGNError(w http.ResponseWriter, variant, assembly, msg string) {
	writeJSON(w, http.StatusOK, gnErrorBody(variant, assembly, msg))
}

// requireAssembly extracts and validates the assembly from the URL path.
func (s *Server) requireAssembly(w http.ResponseWriter, r *http.Request) *assemblyContext {
	name := r.PathValue("assembly")
	if name == "" {
		if ctx := s.soleAssembly(); ctx != nil {
			return ctx
		}
	}
	ctx := s.getAssembly(name)
	if ctx == nil {
		loaded := make([]string, 0)
		s.mu.RLock()
		for k := range s.assemblies {
			loaded = append(loaded, k)
		}
		s.mu.RUnlock()
		writeError(w, http.StatusNotFound, "assembly "+name+" not loaded (available: "+strings.Join(loaded, ", ")+")")
		return nil
	}
	return ctx
}
