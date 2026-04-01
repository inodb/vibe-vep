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
	"github.com/inodb/vibe-vep/internal/datasource/pfam"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// assemblyContext holds the annotator and sources for one assembly.
type assemblyContext struct {
	annotator *annotate.Annotator
	sources   []annotate.AnnotationSource
	cache     *cache.Cache
	pfam      *pfam.Store
	assembly  string // normalized: "GRCh38"
}

// Server is the HTTP annotation server.
type Server struct {
	mu         sync.RWMutex
	assemblies map[string]*assemblyContext // keyed by lowercase: "grch38"
	logger     *zap.Logger
	version    string
}

// New creates a new Server.
func New(logger *zap.Logger, version string) *Server {
	return &Server{
		assemblies: make(map[string]*assemblyContext),
		logger:     logger,
		version:    version,
	}
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

// getAssembly returns the assembly context for the given name (case-insensitive).
func (s *Server) getAssembly(name string) *assemblyContext {
	s.mu.RLock()
	defer s.mu.RUnlock()
	return s.assemblies[strings.ToLower(name)]
}

// Handler returns the HTTP handler with all routes registered.
func (s *Server) Handler() http.Handler {
	mux := http.NewServeMux()

	// Health and info endpoints.
	mux.HandleFunc("GET /health", s.handleHealth)
	mux.HandleFunc("GET /info", s.handleInfo)

	// Ensembl VEP compatibility endpoints.
	mux.HandleFunc("GET /ensembl/{assembly}/vep/human/region/{region}/{allele}", s.handleEnsemblRegionGet)
	mux.HandleFunc("POST /ensembl/{assembly}/vep/human/region", s.handleEnsemblRegionPost)
	mux.HandleFunc("GET /ensembl/{assembly}/vep/human/hgvs/{notation...}", s.handleEnsemblHGVSGet)
	mux.HandleFunc("POST /ensembl/{assembly}/vep/human/hgvs", s.handleEnsemblHGVSPost)

	// Genome-nexus compatibility endpoints.
	mux.HandleFunc("GET /genome-nexus/{assembly}/annotation/genomic/{genomicLocation}", s.handleGNGenomicGet)
	mux.HandleFunc("POST /genome-nexus/{assembly}/annotation/genomic", s.handleGNGenomicPost)
	mux.HandleFunc("GET /genome-nexus/{assembly}/annotation/{variant...}", s.handleGNHGVSGet)
	mux.HandleFunc("POST /genome-nexus/{assembly}/annotation", s.handleGNHGVSPost)

	// Ensembl internal endpoints (transcript/PFAM lookups for frontend).
	mux.HandleFunc("GET /genome-nexus/{assembly}/ensembl/canonical-transcript/hgnc/{hugoSymbol}", s.handleEnsemblCanonicalByHugo)
	mux.HandleFunc("GET /genome-nexus/{assembly}/ensembl/transcript/{transcriptId}", s.handleEnsemblTranscriptByID)
	mux.HandleFunc("POST /genome-nexus/{assembly}/pfam/domain", s.handlePfamDomainPost)
	mux.HandleFunc("GET /genome-nexus/{assembly}/pfam/domain/{pfamAccession}", s.handlePfamDomainGet)

	// POST ensembl/transcript filter (used by mutation mapper to fetch transcript details).
	mux.HandleFunc("POST /genome-nexus/{assembly}/ensembl/transcript", s.handleEnsemblTranscriptPost)

	// Stub endpoints — return empty responses so the frontend mutation mapper
	// doesn't hang waiting for data we don't have yet.
	emptyArray := func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Content-Type", "application/json")
		w.Write([]byte("[]\n"))
	}
	mux.HandleFunc("GET /genome-nexus/{assembly}/ensembl/canonical-gene/hgnc/{hugoSymbol}", s.handleEnsemblCanonicalGeneByHugo)
	mux.HandleFunc("GET /genome-nexus/{assembly}/ensembl/canonical-gene/entrez/{entrezGeneId}", s.handleEnsemblCanonicalGeneByEntrez)
	mux.HandleFunc("GET /genome-nexus/{assembly}/cancer_hotspots/transcript/{transcriptId}", s.handleCancerHotspotsTranscript)

	// Stub endpoints — return empty responses for data we don't have yet.
	mux.HandleFunc("POST /genome-nexus/{assembly}/cancer_hotspots/genomic", emptyArray)
	mux.HandleFunc("GET /genome-nexus/{assembly}/ptm/experimental", emptyArray)
	mux.HandleFunc("POST /genome-nexus/{assembly}/ptm/experimental", emptyArray)
	mux.HandleFunc("GET /genome-nexus/{assembly}/ensembl/xrefs", emptyArray)

	return corsMiddleware(mux)
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
		Name             string `json:"name"`
		TranscriptCount  int    `json:"transcript_count"`
		AnnotationSources int   `json:"annotation_sources"`
	}

	assemblies := make([]assemblyInfo, 0, len(s.assemblies))
	for _, ctx := range s.assemblies {
		assemblies = append(assemblies, assemblyInfo{
			Name:             ctx.assembly,
			TranscriptCount:  ctx.cache.TranscriptCount(),
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
func writeError(w http.ResponseWriter, status int, msg string) {
	writeJSON(w, status, map[string]string{"error": msg})
}

// writeGNError writes a genome-nexus style error response that includes
// assembly_name and variant fields. This allows clients (like the frontend)
// to detect the genome build even when annotation fails.
func writeGNError(w http.ResponseWriter, variant, assembly, msg string) {
	writeJSON(w, http.StatusOK, map[string]interface{}{
		"variant":                variant,
		"originalVariantQuery":   variant,
		"assembly_name":          assembly,
		"successfully_annotated": false,
		"errorMessage":           msg,
	})
}

// requireAssembly extracts and validates the assembly from the URL path.
func (s *Server) requireAssembly(w http.ResponseWriter, r *http.Request) *assemblyContext {
	name := r.PathValue("assembly")
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
