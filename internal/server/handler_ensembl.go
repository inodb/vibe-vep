package server

import (
	"encoding/json"
	"fmt"
	"net/http"
	"strconv"
	"strings"

	"go.uber.org/zap"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/output"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// handleEnsemblRegionGet handles GET /ensembl/{assembly}/vep/human/region/{region}/{allele}
// Region format: "7:140753336-140753336:1" (chrom:start-end:strand)
func (s *Server) handleEnsemblRegionGet(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	region := r.PathValue("region")
	allele := r.PathValue("allele")

	v, err := parseEnsemblRegion(region, allele)
	if err != nil {
		writeError(w, http.StatusBadRequest, err.Error())
		return
	}

	input := region + "/" + allele
	anns, err := s.annotateVariant(ctx, v)
	if err != nil {
		s.logger.Error("annotation error", zap.Error(err), zap.String("input", input))
		writeError(w, http.StatusInternalServerError, "annotation failed: "+err.Error())
		return
	}

	data, err := output.MarshalVEPAnnotation(input, v, anns, ctx.assembly)
	if err != nil {
		writeError(w, http.StatusInternalServerError, "marshal error: "+err.Error())
		return
	}

	w.Header().Set("Content-Type", "application/json")
	// VEP REST API returns a JSON array with one element for single-variant GET.
	w.Write([]byte("["))
	w.Write(data)
	w.Write([]byte("]\n"))
}

// handleEnsemblRegionPost handles POST /ensembl/{assembly}/vep/human/region
// Body: {"variants": ["7:140753336-140753336:1/T", ...]}
func (s *Server) handleEnsemblRegionPost(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	var body struct {
		Variants []string `json:"variants"`
	}
	if err := json.NewDecoder(r.Body).Decode(&body); err != nil {
		writeError(w, http.StatusBadRequest, "invalid JSON body: "+err.Error())
		return
	}
	if len(body.Variants) == 0 {
		writeError(w, http.StatusBadRequest, "variants array is empty")
		return
	}

	results := make([]json.RawMessage, 0, len(body.Variants))
	for _, input := range body.Variants {
		// Format: "7:140753336-140753336:1/T" → region="7:140753336-140753336:1", allele="T"
		lastSlash := strings.LastIndex(input, "/")
		if lastSlash < 0 {
			writeError(w, http.StatusBadRequest, fmt.Sprintf("invalid variant format %q (expected region/allele)", input))
			return
		}
		region := input[:lastSlash]
		allele := input[lastSlash+1:]

		v, err := parseEnsemblRegion(region, allele)
		if err != nil {
			writeError(w, http.StatusBadRequest, fmt.Sprintf("variant %q: %s", input, err.Error()))
			return
		}

		anns, err := s.annotateVariant(ctx, v)
		if err != nil {
			s.logger.Error("annotation error", zap.Error(err), zap.String("input", input))
			writeError(w, http.StatusInternalServerError, "annotation failed: "+err.Error())
			return
		}

		data, err := output.MarshalVEPAnnotation(input, v, anns, ctx.assembly)
		if err != nil {
			writeError(w, http.StatusInternalServerError, "marshal error: "+err.Error())
			return
		}
		results = append(results, data)
	}

	w.Header().Set("Content-Type", "application/json")
	writeJSONArray(w, results)
}

// handleEnsemblHGVSGet handles GET /ensembl/{assembly}/vep/human/hgvs/{notation}
func (s *Server) handleEnsemblHGVSGet(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	notation := r.PathValue("notation")
	variants, err := s.resolveHGVS(ctx, notation)
	if err != nil {
		writeError(w, http.StatusBadRequest, err.Error())
		return
	}

	results := make([]json.RawMessage, 0, len(variants))
	for _, v := range variants {
		anns, err := s.annotateVariant(ctx, v)
		if err != nil {
			s.logger.Error("annotation error", zap.Error(err), zap.String("input", notation))
			writeError(w, http.StatusInternalServerError, "annotation failed: "+err.Error())
			return
		}

		data, err := output.MarshalVEPAnnotation(notation, v, anns, ctx.assembly)
		if err != nil {
			writeError(w, http.StatusInternalServerError, "marshal error: "+err.Error())
			return
		}
		results = append(results, data)
	}

	w.Header().Set("Content-Type", "application/json")
	writeJSONArray(w, results)
}

// handleEnsemblHGVSPost handles POST /ensembl/{assembly}/vep/human/hgvs
// Body: {"hgvs_notations": ["ENST00000311936:c.35G>T", ...]}
func (s *Server) handleEnsemblHGVSPost(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	var body struct {
		HGVSNotations []string `json:"hgvs_notations"`
	}
	if err := json.NewDecoder(r.Body).Decode(&body); err != nil {
		writeError(w, http.StatusBadRequest, "invalid JSON body: "+err.Error())
		return
	}
	if len(body.HGVSNotations) == 0 {
		writeError(w, http.StatusBadRequest, "hgvs_notations array is empty")
		return
	}

	results := make([]json.RawMessage, 0, len(body.HGVSNotations))
	for _, notation := range body.HGVSNotations {
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

			data, err := output.MarshalVEPAnnotation(notation, v, anns, ctx.assembly)
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

// parseEnsemblRegion parses "7:140753336-140753336:1" into a vcf.Variant.
// The format is chrom:start-end:strand. For SNPs start==end.
// The ref allele is inferred from context or left empty (VEP convention for region queries).
func parseEnsemblRegion(region, allele string) (*vcf.Variant, error) {
	// Split by ":"
	parts := strings.Split(region, ":")
	if len(parts) < 2 || len(parts) > 3 {
		return nil, fmt.Errorf("invalid region format %q (expected chrom:start-end or chrom:start-end:strand)", region)
	}

	chrom := strings.TrimPrefix(parts[0], "chr")

	// Parse start-end
	rangeParts := strings.Split(parts[1], "-")
	start, err := strconv.ParseInt(rangeParts[0], 10, 64)
	if err != nil {
		return nil, fmt.Errorf("invalid start position in region %q", region)
	}

	end := start
	if len(rangeParts) == 2 {
		end, err = strconv.ParseInt(rangeParts[1], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("invalid end position in region %q", region)
		}
	}

	// For SNPs (start==end), ref is inferred as single base at that position.
	// We leave ref empty since the annotator works with position + alt.
	// For deletions (end > start), the deleted length gives the ref length.
	ref := ""
	if end > start {
		// Deletion or multi-base: ref is a placeholder of the right length.
		ref = strings.Repeat("N", int(end-start+1))
	}

	return &vcf.Variant{
		Chrom: chrom,
		Pos:   start,
		Ref:   ref,
		Alt:   allele,
	}, nil
}

// resolveHGVS resolves an HGVS notation string to genomic variants.
func (s *Server) resolveHGVS(ctx *assemblyContext, notation string) ([]*vcf.Variant, error) {
	spec, err := annotate.ParseVariantSpec(notation)
	if err != nil {
		return nil, err
	}

	switch spec.Type {
	case annotate.SpecGenomic:
		return []*vcf.Variant{{Chrom: spec.Chrom, Pos: spec.Pos, Ref: spec.Ref, Alt: spec.Alt}}, nil
	case annotate.SpecHGVSc:
		variants, _, err := annotate.ReverseMapHGVScWithWarning(ctx.cache, spec.TranscriptID, spec.CDSChange)
		return variants, err
	case annotate.SpecHGVSg:
		return annotate.ResolveHGVSg(ctx.cache, spec.Chrom, spec.GenomicChange)
	case annotate.SpecProtein:
		return annotate.ReverseMapProteinChange(ctx.cache, spec.GeneName, spec.RefAA, spec.Position, spec.AltAA)
	default:
		return nil, fmt.Errorf("unsupported variant spec type for %q", notation)
	}
}

// writeJSONArray writes a JSON array from pre-marshaled elements.
func writeJSONArray(w http.ResponseWriter, items []json.RawMessage) {
	w.Write([]byte("["))
	for i, item := range items {
		if i > 0 {
			w.Write([]byte(","))
		}
		w.Write(item)
	}
	w.Write([]byte("]\n"))
}
