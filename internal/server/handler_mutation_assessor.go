package server

import (
	"net/http"

	"github.com/inodb/vibe-vep/internal/output"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Mutation Assessor predicts the functional impact of an amino acid
// substitution, so its predictions are keyed by the protein -- a UniProt
// accession -- and the substitution, not by genomic coordinate. Answering for
// a variant therefore means annotating it first, then asking what protein
// change it makes.

// gnMutationAssessor is one prediction, in genome-nexus's wire shape.
type gnMutationAssessor struct {
	UniprotID                  string   `json:"uniprotId,omitempty"`
	SV                         int      `json:"sv,omitempty"`
	HGVSpShort                 string   `json:"hgvspShort,omitempty"`
	FunctionalImpactScore      *float64 `json:"functionalImpactScore,omitempty"`
	FunctionalImpactPrediction string   `json:"functionalImpactPrediction,omitempty"`
	MSA                        string   `json:"msa,omitempty"`
	MAV                        int      `json:"mav,omitempty"`
}

// handleMutationAssessorGet handles GET {prefix}/mutation_assessor/{variant}.
func (s *Server) handleMutationAssessorGet(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}

	variants, err := s.resolveHGVS(ctx, r.PathValue("variant"))
	if err != nil || len(variants) == 0 {
		writeJSON(w, http.StatusOK, gnMutationAssessor{})
		return
	}
	pred, _ := s.mutationAssessorFor(ctx, variants[0])
	writeJSON(w, http.StatusOK, pred)
}

// handleMutationAssessorPost handles POST {prefix}/mutation_assessor.
func (s *Server) handleMutationAssessorPost(w http.ResponseWriter, r *http.Request) {
	ctx := s.requireAssembly(w, r)
	if ctx == nil {
		return
	}
	notations, ok := decodeStringList(w, r)
	if !ok {
		return
	}

	out := make([]gnMutationAssessor, 0, len(notations))
	for _, notation := range notations {
		variants, err := s.resolveHGVS(ctx, notation)
		if err != nil || len(variants) == 0 {
			out = append(out, gnMutationAssessor{})
			continue
		}
		pred, _ := s.mutationAssessorFor(ctx, variants[0])
		out = append(out, pred)
	}
	writeJSON(w, http.StatusOK, out)
}

// mutationAssessorFor looks up the prediction for a variant.
//
// Only a single-base substitution has an amino acid substitution to predict
// for, which is the same restriction genome-nexus applies.
//
// The UniProt accession is taken from the canonical transcript, the same one
// the protein change comes from. genome-nexus takes the accession from the
// first transcript consequence instead, while still taking the protein change
// from the canonical -- an inconsistency that only agrees with itself when the
// canonical happens to sort first.
func (s *Server) mutationAssessorFor(ctx *assemblyContext, v *vcf.Variant) (gnMutationAssessor, bool) {
	if ctx.mutAssess == nil || ctx.uniprot == nil {
		return gnMutationAssessor{}, false
	}
	if len(v.Ref) != 1 || len(v.Alt) != 1 {
		return gnMutationAssessor{}, false
	}

	anns, err := s.annotateVariant(ctx, v)
	if err != nil || len(anns) == 0 {
		return gnMutationAssessor{}, false
	}

	opts := output.GNMarshalOptions{
		IncludeAnnotationSummary: true,
		OverrideTranscripts:      s.OverrideTranscripts(""),
	}
	cancer, curated := s.CancerGenes()
	opts.CancerGenes, opts.CuratedGenes = cancer, curated

	data, err := output.MarshalGNAnnotation("", v, anns, ctx.assembly, opts)
	if err != nil {
		return gnMutationAssessor{}, false
	}
	parsed, err := output.ParseGNAnnotation(data)
	if err != nil || parsed.AnnotationSummary == nil ||
		parsed.AnnotationSummary.TranscriptConsequenceSummary == nil {
		return gnMutationAssessor{}, false
	}

	summary := parsed.AnnotationSummary.TranscriptConsequenceSummary
	uniprotID := ctx.uniprot.LookupByTranscript(summary.TranscriptID)
	if uniprotID == "" || summary.HGVSpShort == "" {
		return gnMutationAssessor{}, false
	}

	entry, ok := ctx.mutAssess.Lookup(uniprotID, summary.HGVSpShort)
	if !ok {
		return gnMutationAssessor{}, false
	}

	score := entry.FunctionalImpactScore
	return gnMutationAssessor{
		UniprotID:                  entry.UniprotID,
		SV:                         entry.SV,
		HGVSpShort:                 entry.HGVSpShort,
		FunctionalImpactScore:      &score,
		FunctionalImpactPrediction: entry.FunctionalImpactPrediction,
		MSA:                        entry.MSA,
		MAV:                        entry.MAV,
	}, true
}
