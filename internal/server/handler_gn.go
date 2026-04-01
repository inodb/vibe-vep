package server

import (
	"encoding/json"
	"fmt"
	"net/http"
	"strconv"
	"strings"

	"go.uber.org/zap"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/datasource/myvariantinfo"
	"github.com/inodb/vibe-vep/internal/input"
	"github.com/inodb/vibe-vep/internal/output"
	"github.com/inodb/vibe-vep/internal/vcf"
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
	s.enrichMyVariantInfo(&opts, v, ctx.assembly, anns)
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

	baseOpts := parseGNMarshalOptions(r)
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

		opts := baseOpts
		s.enrichMyVariantInfo(&opts, v, ctx.assembly, anns)
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
		// Return a partial response with assembly_name (like genome-nexus does)
		// so that clients can still detect the genome build.
		writeGNError(w, notation, ctx.assembly, err.Error())
		return
	}

	// GN returns a single annotation object for single-variant GET.
	if len(variants) == 0 {
		writeGNError(w, notation, ctx.assembly, "could not resolve variant: "+notation)
		return
	}

	v := variants[0]
	anns, err := s.annotateVariant(ctx, v)
	if err != nil {
		s.logger.Error("annotation error", zap.Error(err), zap.String("input", notation))
		writeGNError(w, notation, ctx.assembly, "annotation failed: "+err.Error())
		return
	}

	opts := parseGNMarshalOptions(r)
	s.enrichMyVariantInfo(&opts, v, ctx.assembly, anns)
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

	baseOpts := parseGNMarshalOptions(r)
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

			opts := baseOpts
			s.enrichMyVariantInfo(&opts, v, ctx.assembly, anns)
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

// enrichMyVariantInfo populates opts.MyVariantInfoData from local annotation extras
// (gnomAD, dbSNP) when available, falling back to myvariant.info API when local data
// is not available.
func (s *Server) enrichMyVariantInfo(opts *output.GNMarshalOptions, v *vcf.Variant, assembly string, anns []*annotate.Annotation) {
	if !opts.IncludeMyVariantInfo {
		return
	}

	ref := v.Ref
	if ref == "" {
		ref = "-"
	}
	alt := v.Alt
	if alt == "" {
		alt = "-"
	}

	// Try local data first (from genomic index extras).
	src := firstAnnotationWithExtras(anns)
	if src != nil {
		gnomadAF := src.GetExtraKey("gnomad.af")
		dbsnpID := src.GetExtraKey("dbsnp.id")
		if gnomadAF != "" || dbsnpID != "" {
			mvi := &output.GNMyVariantInfo{
				Variant: fmt.Sprintf("chr%s:g.%d%s>%s", v.Chrom, v.Pos, ref, alt),
				Query:   fmt.Sprintf("chr%s:g.%d%s>%s", v.Chrom, v.Pos, ref, alt),
				Hgvs:    fmt.Sprintf("chr%s:g.%d%s>%s", v.Chrom, v.Pos, ref, alt),
				Vcf: &output.GNVcf{
					Ref:      ref,
					Alt:      alt,
					Position: fmt.Sprintf("%d", v.Pos),
				},
			}
			if dbsnpID != "" {
				mvi.Dbsnp = &output.GNDbsnp{Rsid: dbsnpID}
			}
			if gnomadAF != "" {
				gnomadAC := src.GetExtraKey("gnomad.ac")
				gnomadAN := src.GetExtraKey("gnomad.an")
				gnomadNhomalt := src.GetExtraKey("gnomad.nhomalt")
				mvi.GnomadExome = buildLocalGnomad(gnomadAF, gnomadAC, gnomadAN, gnomadNhomalt)
			}
			opts.MyVariantInfoData = &output.GNMyVariantInfoAnnotation{Annotation: mvi}
			return
		}
	}

	// Fall back to myvariant.info API.
	if s.myVariantClient == nil {
		return
	}
	hgvsg := fmt.Sprintf("%s:g.%d%s>%s", v.Chrom, v.Pos, ref, alt)
	resp, err := s.myVariantClient.Fetch(hgvsg, assembly)
	if err != nil {
		s.logger.Warn("myvariant.info fetch failed", zap.Error(err), zap.String("hgvsg", hgvsg))
		return
	}
	if resp == nil {
		return
	}
	opts.MyVariantInfoData = &output.GNMyVariantInfoAnnotation{
		Annotation: convertMyVariantInfo(resp.Annotation),
	}
}

// buildLocalGnomad constructs a GNGnomad from local genomic index values.
// Only overall frequency is available locally (not per-population).
func buildLocalGnomad(af, ac, an, nhomalt string) *output.GNGnomad {
	g := &output.GNGnomad{
		AlleleFrequency: map[string]interface{}{},
		AlleleCount:     map[string]interface{}{},
		AlleleNumber:    map[string]interface{}{},
		Homozygotes:     map[string]interface{}{},
	}
	if v, err := strconv.ParseFloat(af, 64); err == nil {
		g.AlleleFrequency["af"] = v
	}
	if v, err := strconv.Atoi(ac); err == nil {
		g.AlleleCount["ac"] = v
	}
	if v, err := strconv.Atoi(an); err == nil {
		g.AlleleNumber["an"] = v
	}
	if v, err := strconv.Atoi(nhomalt); err == nil {
		g.Homozygotes["hom"] = v
	}
	return g
}

// firstAnnotationWithExtras returns the first annotation with extras populated.
func firstAnnotationWithExtras(anns []*annotate.Annotation) *annotate.Annotation {
	for _, ann := range anns {
		if len(ann.Extra) > 0 {
			return ann
		}
	}
	if len(anns) > 0 {
		return anns[0]
	}
	return nil
}

// convertMyVariantInfo converts from the myvariantinfo package types to the output package types.
func convertMyVariantInfo(src *myvariantinfo.MyVariantInfo) *output.GNMyVariantInfo {
	if src == nil {
		return nil
	}
	dst := &output.GNMyVariantInfo{
		Variant: src.Variant,
		Query:   src.Query,
		Hgvs:    src.Hgvs,
	}
	if src.Dbsnp != nil {
		dst.Dbsnp = &output.GNDbsnp{Rsid: src.Dbsnp.Rsid}
	}
	if src.GnomadExome != nil {
		dst.GnomadExome = &output.GNGnomad{
			AlleleCount:     src.GnomadExome.AlleleCount,
			AlleleFrequency: src.GnomadExome.AlleleFrequency,
			AlleleNumber:    src.GnomadExome.AlleleNumber,
			Homozygotes:     src.GnomadExome.Homozygotes,
		}
	}
	if src.GnomadGenome != nil {
		dst.GnomadGenome = &output.GNGnomad{
			AlleleCount:     src.GnomadGenome.AlleleCount,
			AlleleFrequency: src.GnomadGenome.AlleleFrequency,
			AlleleNumber:    src.GnomadGenome.AlleleNumber,
			Homozygotes:     src.GnomadGenome.Homozygotes,
		}
	}
	if src.Vcf != nil {
		dst.Vcf = &output.GNVcf{
			Ref:      src.Vcf.Ref,
			Alt:      src.Vcf.Alt,
			Position: src.Vcf.Position,
		}
	}
	return dst
}

// parseGNMarshalOptions extracts GN marshal options from query parameters.
// Supports both comma-separated (?fields=a,b) and repeated (?fields=a&fields=b).
func parseGNMarshalOptions(r *http.Request) output.GNMarshalOptions {
	fields := strings.Join(r.URL.Query()["fields"], ",")
	return output.GNMarshalOptions{
		IncludeAnnotationSummary: strings.Contains(fields, "annotation_summary"),
		IncludeClinVar:           strings.Contains(fields, "clinvar"),
		IncludeHotspots:          strings.Contains(fields, "hotspots"),
		IncludeSignal:            strings.Contains(fields, "signal"),
		IncludeMyVariantInfo:     strings.Contains(fields, "my_variant_info"),
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

