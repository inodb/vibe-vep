package server

import (
	"encoding/json"
	"net/http"
	"regexp"
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/datasource/hotspots"
	"github.com/inodb/vibe-vep/internal/input"
	"github.com/inodb/vibe-vep/internal/output"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// The /cancer_hotspots endpoints answer from the hotspots table directly.
// Asked about a transcript they report every hotspot on it; asked about a
// variant they report only the hotspots that variant could actually be,
// which means agreeing on both the residue and the kind of mutation.

// gnHotspot is one hotspot in genome-nexus's wire shape. It is shared with
// the "hotspots" field of an annotation, which reports the same records.
type gnHotspot = output.GNHotspotEntry

// proteinLocation is a span of residues on a transcript, with the kind of
// mutation being asked about.
type proteinLocation struct {
	TranscriptID string `json:"transcriptId"`
	Start        int    `json:"start"`
	End          int    `json:"end"`
	MutationType string `json:"mutationType,omitempty"`
}

// aggregatedHotspots is what the batch endpoints answer with: the question
// repeated back alongside its hotspots, so a caller can match them up.
type aggregatedHotspots struct {
	GenomicLocation     *input.GenomicLocation `json:"genomicLocation,omitempty"`
	Variant             string                 `json:"variant,omitempty"`
	TranscriptID        string                 `json:"transcriptId,omitempty"`
	TranscriptIDVersion string                 `json:"transcriptIdVersion,omitempty"`
	ProteinLocation     *proteinLocation       `json:"proteinLocation,omitempty"`
	Hotspots            []gnHotspot            `json:"hotspots"`
}

// hotspotVersion resolves the ?version= parameter. genome-nexus defaults it to
// v3 for real -- unlike most of its parameters, whose documented defaults are
// Swagger-only.
func hotspotVersion(r *http.Request) string {
	v := strings.ToLower(strings.TrimSpace(r.URL.Query().Get("version")))
	if v == "" {
		return "v3"
	}
	return v
}

// versionMatches is genome-nexus's deliberately asymmetric rule: v3 accepts
// every row, while v2 accepts only rows that say v2 or say nothing.
func versionMatches(requested, rowVersion string) bool {
	if requested == "v3" {
		return true
	}
	return rowVersion == "" || strings.EqualFold(rowVersion, requested)
}

// renderHotspots converts stored hotspots to the wire shape, dropping those
// the requested version excludes.
func renderHotspots(spots []hotspots.Hotspot, version string) []gnHotspot {
	out := make([]gnHotspot, 0, len(spots))
	for _, h := range spots {
		if !versionMatches(version, h.Version) {
			continue
		}
		out = append(out, gnHotspot{
			HugoSymbol:          h.HugoSymbol,
			TranscriptID:        h.TranscriptID,
			TranscriptIDVersion: h.TranscriptIDVersion,
			Residue:             h.Residue,
			Version:             h.Version,
			TumorCount:          h.TumorCount,
			Type:                h.Type,
			MissenseCount:       h.MissenseCount,
			TruncatingCount:     h.TruncatingCount,
			InframeCount:        h.InframeCount,
			SpliceCount:         h.SpliceCount,
		})
	}
	return out
}

// digitsRe pulls the residue numbers out of a hotspot residue such as "Q61"
// or "V600_K601".
var digitsRe = regexp.MustCompile(`\d+`)

// residueRange returns the first and last residue numbers a hotspot covers.
func residueRange(residue string) (int, int, bool) {
	found := digitsRe.FindAllString(residue, -1)
	if len(found) == 0 {
		return 0, 0, false
	}
	start, err := strconv.Atoi(found[0])
	if err != nil {
		return 0, 0, false
	}
	end := start
	if len(found) > 1 {
		if e, err := strconv.Atoi(found[len(found)-1]); err == nil {
			end = e
		}
	}
	return start, end, true
}

// overlaps reports whether two residue spans touch.
func overlaps(aStart, aEnd, bStart, bEnd int) bool {
	return aStart <= bEnd && bStart <= aEnd
}

// hotspotTypeMatches reports whether a hotspot of this type can be the kind of
// mutation these classifications describe. A type nobody recognises matches
// anything, which is how genome-nexus treats it.
func hotspotTypeMatches(hotspotType string, classifications []string) bool {
	want := func(substr string) bool {
		for _, c := range classifications {
			if strings.Contains(strings.ToLower(c), substr) {
				return true
			}
		}
		return false
	}

	switch {
	case strings.Contains(hotspotType, "single residue"), strings.Contains(hotspotType, "3d"):
		return want("missense")
	case strings.Contains(hotspotType, "splice site"):
		return want("splice")
	case strings.Contains(hotspotType, "in-frame"):
		return want("inframe") || want("in_frame")
	default:
		return true
	}
}

// hotspotsForVariant reports the hotspots a variant could be: on one of its
// transcripts, covering the residues it changes, and of a kind its
// consequence terms allow.
func (s *Server) hotspotsForVariant(ctx *assemblyContext, v *vcf.Variant, version string) []gnHotspot {
	if ctx.hotspots == nil {
		return []gnHotspot{}
	}
	anns, err := s.annotateVariant(ctx, v)
	if err != nil {
		return []gnHotspot{}
	}

	ref, alt := v.Ref, v.Alt
	if ref == "" {
		ref = "-"
	}
	if alt == "" {
		alt = "-"
	}
	variantType := output.GNVariantType(ref, alt)
	inframe := output.GNIsInframe(ref, alt)

	seen := make(map[string]bool)
	matched := make([]gnHotspot, 0)
	for _, ann := range anns {
		for _, h := range hotspotsForAnnotation(ctx, ann, variantType, inframe, version) {
			key := h.TranscriptID + "|" + h.Residue + "|" + h.Type
			if seen[key] {
				continue
			}
			seen[key] = true
			matched = append(matched, h)
		}
	}
	return matched
}

// hotspotsForAnnotation reports the hotspots one transcript's view of a
// variant could be: covering the residues it changes, and of a kind its
// consequence terms allow.
func hotspotsForAnnotation(ctx *assemblyContext, ann *annotate.Annotation, variantType string, inframe bool, version string) []gnHotspot {
	if ctx.hotspots == nil {
		return nil
	}
	start, end := ann.ProteinStart, ann.ProteinEnd
	if start == 0 {
		start = ann.ProteinPosition
	}
	if end == 0 {
		end = start
	}
	if start == 0 {
		return nil
	}

	classes := output.GNVariantClassifications(
		strings.Split(ann.Consequence, ","), variantType, inframe)

	var matched []gnHotspot
	for _, h := range renderHotspots(ctx.hotspots.LookupAll(ann.TranscriptID), version) {
		hStart, hEnd, ok := residueRange(h.Residue)
		if !ok || !overlaps(hStart, hEnd, int(start), int(end)) {
			continue
		}
		if !hotspotTypeMatches(h.Type, classes) {
			continue
		}
		matched = append(matched, h)
	}
	return matched
}

// handleHotspotsByTranscriptGet handles GET
// {prefix}/cancer_hotspots/transcript/{transcriptId}.
func (s *Server) handleHotspotsByTranscriptGet(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}
	writeJSON(w, http.StatusOK, s.transcriptHotspots(ctx, r.PathValue("transcriptId"), hotspotVersion(r)))
}

// handleHotspotsByTranscriptPost handles POST
// {prefix}/cancer_hotspots/transcript.
func (s *Server) handleHotspotsByTranscriptPost(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}
	ids, ok := decodeStringList(w, r)
	if !ok {
		return
	}

	version := hotspotVersion(r)
	out := make([]aggregatedHotspots, 0, len(ids))
	for _, id := range ids {
		out = append(out, aggregatedHotspots{
			TranscriptID: stripTxVersion(id),
			Hotspots:     s.transcriptHotspots(ctx, id, version),
		})
	}
	writeJSON(w, http.StatusOK, out)
}

// transcriptHotspots reports every hotspot on a transcript.
func (s *Server) transcriptHotspots(ctx *assemblyContext, transcriptID, version string) []gnHotspot {
	if ctx.hotspots == nil {
		return []gnHotspot{}
	}
	return renderHotspots(ctx.hotspots.LookupAll(transcriptID), version)
}

// handleHotspotsByHgvsGet handles GET {prefix}/cancer_hotspots/hgvs/{variant}.
func (s *Server) handleHotspotsByHgvsGet(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}
	variants, err := s.resolveHGVS(ctx, r.PathValue("variant"))
	if err != nil || len(variants) == 0 {
		writeJSON(w, http.StatusOK, []gnHotspot{})
		return
	}
	writeJSON(w, http.StatusOK, s.hotspotsForVariant(ctx, variants[0], hotspotVersion(r)))
}

// handleHotspotsByHgvsPost handles POST {prefix}/cancer_hotspots/hgvs.
func (s *Server) handleHotspotsByHgvsPost(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}
	notations, ok := decodeStringList(w, r)
	if !ok {
		return
	}

	version := hotspotVersion(r)
	out := make([]aggregatedHotspots, 0, len(notations))
	for _, notation := range notations {
		entry := aggregatedHotspots{Variant: notation, Hotspots: []gnHotspot{}}
		if variants, err := s.resolveHGVS(ctx, notation); err == nil && len(variants) > 0 {
			entry.Hotspots = s.hotspotsForVariant(ctx, variants[0], version)
		}
		out = append(out, entry)
	}
	writeJSON(w, http.StatusOK, out)
}

// handleHotspotsByGenomicGet handles GET
// {prefix}/cancer_hotspots/genomic/{genomicLocation}.
func (s *Server) handleHotspotsByGenomicGet(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}
	gl, err := parseCommaGenomicLocation(r.PathValue("genomicLocation"))
	if err != nil {
		writeError(w, http.StatusBadRequest, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, s.hotspotsForVariant(ctx, gl.ToVariant(), hotspotVersion(r)))
}

// handleHotspotsByGenomicPost handles POST {prefix}/cancer_hotspots/genomic.
func (s *Server) handleHotspotsByGenomicPost(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}
	var locations []input.GenomicLocation
	if err := json.NewDecoder(r.Body).Decode(&locations); err != nil {
		writeError(w, http.StatusBadRequest, "invalid JSON body: "+err.Error())
		return
	}

	version := hotspotVersion(r)
	out := make([]aggregatedHotspots, 0, len(locations))
	for i := range locations {
		gl := locations[i]
		out = append(out, aggregatedHotspots{
			GenomicLocation: &gl,
			Hotspots:        s.hotspotsForVariant(ctx, gl.ToVariant(), version),
		})
	}
	writeJSON(w, http.StatusOK, out)
}

// handleHotspotsByProteinLocationsPost handles POST
// {prefix}/cancer_hotspots/proteinLocations, which asks about residue spans
// directly rather than about a variant.
func (s *Server) handleHotspotsByProteinLocationsPost(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}
	var locations []proteinLocation
	if err := json.NewDecoder(r.Body).Decode(&locations); err != nil {
		writeError(w, http.StatusBadRequest, "invalid JSON body: "+err.Error())
		return
	}

	version := hotspotVersion(r)
	out := make([]aggregatedHotspots, 0, len(locations))
	for i := range locations {
		loc := locations[i]
		matched := make([]gnHotspot, 0)
		for _, h := range s.transcriptHotspots(ctx, loc.TranscriptID, version) {
			hStart, hEnd, ok := residueRange(h.Residue)
			if !ok || !overlaps(hStart, hEnd, loc.Start, loc.End) {
				continue
			}
			if loc.MutationType != "" && !hotspotTypeMatches(h.Type, []string{loc.MutationType}) {
				continue
			}
			matched = append(matched, h)
		}
		out = append(out, aggregatedHotspots{
			TranscriptID:    stripTxVersion(loc.TranscriptID),
			ProteinLocation: &loc,
			Hotspots:        matched,
		})
	}
	writeJSON(w, http.StatusOK, out)
}
