package server

import (
	"encoding/json"
	"fmt"
	"net/http"
	"strconv"
	"strings"

	"go.uber.org/zap"

	"github.com/inodb/vibe-vep/internal/input"
	"github.com/inodb/vibe-vep/internal/output"
)

// handleGNGenomicGet handles GET /genome-nexus/{assembly}/annotation/genomic/{genomicLocation}
// genomicLocation format: "7,140753336,140753336,A,T" (comma-separated)
func (s *Server) handleGNGenomicGet(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	locStr := r.PathValue("genomicLocation")
	gl, err := parseCommaGenomicLocation(locStr)
	if err != nil {
		writeError(w, http.StatusBadRequest, err.Error())
		return
	}

	v := gl.ToVariant()
	inputLabel := gl.FormatInput()

	anns, err := s.annotateVariant(ctx, v)
	if err != nil {
		s.logger.Error("annotation error", zap.Error(err), zap.String("input", inputLabel))
		writeError(w, http.StatusInternalServerError, "annotation failed: "+err.Error())
		return
	}

	opts := parseGNMarshalOptions(r)
	data, err := output.MarshalGNAnnotation(inputLabel, v, anns, ctx.assembly, opts)
	if err != nil {
		writeError(w, http.StatusInternalServerError, "marshal error: "+err.Error())
		return
	}

	w.Header().Set("Content-Type", "application/json")
	w.Write(data)
	w.Write([]byte("\n"))
}

// handleGNGenomicPost handles POST /genome-nexus/{assembly}/annotation/genomic
// Body: JSON array of GenomicLocation objects.
func (s *Server) handleGNGenomicPost(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	var locations []input.GenomicLocation
	if err := json.NewDecoder(r.Body).Decode(&locations); err != nil {
		writeError(w, http.StatusBadRequest, "invalid JSON body: "+err.Error())
		return
	}
	if len(locations) == 0 {
		writeError(w, http.StatusBadRequest, "empty array")
		return
	}

	opts := parseGNMarshalOptions(r)
	results := make([]json.RawMessage, 0, len(locations))
	for _, gl := range locations {
		v := gl.ToVariant()
		inputLabel := gl.FormatInput()

		anns, err := s.annotateVariant(ctx, v)
		if err != nil {
			s.logger.Error("annotation error", zap.Error(err), zap.String("input", inputLabel))
			writeError(w, http.StatusInternalServerError, "annotation failed: "+err.Error())
			return
		}

		data, err := output.MarshalGNAnnotation(inputLabel, v, anns, ctx.assembly, opts)
		if err != nil {
			writeError(w, http.StatusInternalServerError, "marshal error: "+err.Error())
			return
		}
		results = append(results, data)
	}

	w.Header().Set("Content-Type", "application/json")
	writeJSONArray(w, results)
}

// handleGNHGVSGet handles GET /genome-nexus/{assembly}/annotation/{variant}
// variant is an HGVS notation.
func (s *Server) handleGNHGVSGet(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	notation := r.PathValue("variant")
	variants, err := s.resolveHGVS(ctx, notation)
	if err != nil {
		writeError(w, http.StatusBadRequest, err.Error())
		return
	}

	// GN returns a single annotation object for single-variant GET.
	if len(variants) == 0 {
		writeError(w, http.StatusNotFound, "could not resolve variant: "+notation)
		return
	}

	v := variants[0]
	anns, err := s.annotateVariant(ctx, v)
	if err != nil {
		s.logger.Error("annotation error", zap.Error(err), zap.String("input", notation))
		writeError(w, http.StatusInternalServerError, "annotation failed: "+err.Error())
		return
	}

	opts := parseGNMarshalOptions(r)
	data, err := output.MarshalGNAnnotation(notation, v, anns, ctx.assembly, opts)
	if err != nil {
		writeError(w, http.StatusInternalServerError, "marshal error: "+err.Error())
		return
	}

	w.Header().Set("Content-Type", "application/json")
	w.Write(data)
	w.Write([]byte("\n"))
}

// handleGNHGVSPost handles POST /genome-nexus/{assembly}/annotation
// Body: JSON array of HGVS notation strings.
func (s *Server) handleGNHGVSPost(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	var notations []string
	if err := json.NewDecoder(r.Body).Decode(&notations); err != nil {
		writeError(w, http.StatusBadRequest, "invalid JSON body: "+err.Error())
		return
	}
	if len(notations) == 0 {
		writeError(w, http.StatusBadRequest, "empty array")
		return
	}

	opts := parseGNMarshalOptions(r)
	results := make([]json.RawMessage, 0, len(notations))
	for _, notation := range notations {
		variants, err := s.resolveHGVS(ctx, notation)
		if err != nil {
			writeError(w, http.StatusBadRequest, fmt.Sprintf("notation %q: %s", notation, err.Error()))
			return
		}

		for _, v := range variants {
			anns, err := s.annotateVariant(ctx, v)
			if err != nil {
				s.logger.Error("annotation error", zap.Error(err), zap.String("input", notation))
				writeError(w, http.StatusInternalServerError, "annotation failed: "+err.Error())
				return
			}

			data, err := output.MarshalGNAnnotation(notation, v, anns, ctx.assembly, opts)
			if err != nil {
				writeError(w, http.StatusInternalServerError, "marshal error: "+err.Error())
				return
			}
			results = append(results, data)
		}
	}

	w.Header().Set("Content-Type", "application/json")
	writeJSONArray(w, results)
}

// parseGNMarshalOptions extracts GN marshal options from query parameters.
func parseGNMarshalOptions(r *http.Request) output.GNMarshalOptions {
	fields := r.URL.Query().Get("fields")
	return output.GNMarshalOptions{
		IncludeAnnotationSummary: strings.Contains(fields, "annotation_summary"),
		IncludeClinVar:           strings.Contains(fields, "clinvar"),
		IncludeHotspots:          strings.Contains(fields, "hotspots"),
	}
}

// parseCommaGenomicLocation parses "7,140753336,140753336,A,T" into a GenomicLocation.
func parseCommaGenomicLocation(s string) (input.GenomicLocation, error) {
	parts := strings.Split(s, ",")
	if len(parts) != 5 {
		return input.GenomicLocation{}, fmt.Errorf("invalid genomic location %q (expected chrom,start,end,ref,alt)", s)
	}

	start, err := strconv.ParseInt(parts[1], 10, 64)
	if err != nil {
		return input.GenomicLocation{}, fmt.Errorf("invalid start in genomic location %q", s)
	}
	end, err := strconv.ParseInt(parts[2], 10, 64)
	if err != nil {
		return input.GenomicLocation{}, fmt.Errorf("invalid end in genomic location %q", s)
	}

	return input.GenomicLocation{
		Chromosome:     parts[0],
		Start:          start,
		End:            end,
		ReferenceAllele: parts[3],
		VariantAllele:  parts[4],
	}, nil
}

