package output

import (
	"encoding/json"
	"sort"
	"strconv"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// MarshalVEPAnnotation builds a VEPVariantAnnotation from a variant and its annotations,
// then marshals it to JSON. This is the standalone equivalent of JSONLWriter.marshalVEP().
func MarshalVEPAnnotation(input string, v *vcf.Variant, anns []*annotate.Annotation, assembly string) ([]byte, error) {
	ref, alt := alleleStrings(v)

	start, end := vepCoordinates(v)

	result := VEPVariantAnnotation{
		Input:         input,
		ID:            annotate.FormatVariantID(v.Chrom, v.Pos, v.Ref, v.Alt),
		SeqRegionName: v.Chrom,
		Start:         start,
		End:           end,
		AlleleString:  ref + "/" + alt,
		Strand:        1,
		AssemblyName:  assembly,
	}

	result.MostSevereConsequence = mostSevereConsequence(anns)

	for _, ann := range anns {
		pStart, pEnd := proteinRange(ann)
		tc := VEPTranscriptConsequence{
			TranscriptID:        stripVersion(ann.TranscriptID),
			GeneID:              ann.GeneID,
			GeneSymbol:          ann.GeneName,
			GeneSymbolSource:    "HGNC",
			Biotype:             ann.Biotype,
			ConsequenceTerms:    splitConsequence(ann.Consequence),
			Impact:              ann.Impact,
			VariantAllele:       ann.Allele,
			AminoAcids:          aminoAcidsField(ann),
			Codons:              ann.CodonChange,
			ProteinStart:        pStart,
			ProteinEnd:          pEnd,
			CDSStart:            ann.CDSPosition,
			CDSEnd:              ann.CDSPosition,
			CDNAStart:           ann.CDNAPosition,
			CDNAEnd:             ann.CDNAPosition,
			HGVSc:               prependTranscriptID(ann.TranscriptID, ann.HGVSc),
			HGVSp:               ann.HGVSp,
			Exon:                ann.ExonNumber,
			Intron:              ann.IntronNumber,
			RefSeqTranscriptIDs: ann.RefSeqIDs,
			HGVSOffset:          ann.HGVSOffset,
		}

		// SIFT/PolyPhen from annotation source extras.
		if s := ann.GetExtraKey("sift.score"); s != "" {
			if f, err := strconv.ParseFloat(s, 64); err == nil {
				tc.SIFTScore = &f
			}
		}
		tc.SIFTPrediction = ann.GetExtraKey("sift.prediction")
		if s := ann.GetExtraKey("polyphen.score"); s != "" {
			if f, err := strconv.ParseFloat(s, 64); err == nil {
				tc.PolyPhenScore = &f
			}
		}
		tc.PolyPhenPrediction = ann.GetExtraKey("polyphen.prediction")

		result.TranscriptConsequences = append(result.TranscriptConsequences, tc)
	}

	sortVEPTranscripts(result.TranscriptConsequences)

	return json.Marshal(result)
}

// proteinRange returns the protein_start/protein_end VEP reports for an
// annotation. These are the UNSHIFTED positions of the change, which differ
// from the 3'-shifted position the HGVSp is written at: VEP reports PBRM1 as
// protein_start 1169 alongside its own p.I1170Sfs*23. Falls back to the HGVS
// position when the unshifted one was not computed (non-coding paths).
func proteinRange(ann *annotate.Annotation) (start, end int64) {
	// Withholding this wherever the HGVSp is withheld was tried and measured
	// worse. VEP does give a protein position to an intronic indel -- CYLD
	// c.2099+1del is p.R700fs -- and genome-nexus builds that short form from
	// the position and the amino acid together, so dropping only the position
	// left it rendering p.Rnullfs on 112 variants to recover 11.
	if annotate.ProteinPositionWithheld(ann.Consequence) {
		return 0, 0
	}
	start = ann.ProteinStart
	if start == 0 {
		start = ann.ProteinPosition
	}
	end = ann.ProteinEnd
	if end == 0 {
		end = start
	}
	return start, end
}

// sortVEPTranscripts orders transcript consequences by ascending transcript ID,
// which is the order Ensembl VEP emits them in (verified across 22 real VEP
// fixtures, with no exceptions).
//
// The order is load-bearing, not cosmetic. genome-nexus breaks a tie between
// equally-severe canonical candidates with bestCandidates.get(0) -- i.e. the
// first one in this list -- so emitting a different order silently changes the
// transcript picked for the MAF. CDKN2A is the clearest case: ENST00000304494
// and ENST00000579755 are both protein_coding frameshift_variant and both get
// flagged canonical, so whichever is listed first wins.
func sortVEPTranscripts(tcs []VEPTranscriptConsequence) {
	sort.SliceStable(tcs, func(i, j int) bool {
		return tcs[i].TranscriptID < tcs[j].TranscriptID
	})
}

// mostSevereConsequence returns the term VEP reports as
// most_severe_consequence: the highest ranked term across every transcript.
//
// Ranking by impact alone and keeping the first annotation to reach the best
// one is not the same thing, because a whole class of terms shares an impact.
// CDKN2A 9:21971192 reaches seventeen transcripts, eleven of them HIGH; three
// say stop_gained and the rest frameshift_variant, and stop_gained outranks
// frameshift_variant. Taking the first HIGH annotation in transcript order
// left frameshift_variant standing as the most severe, and genome-nexus
// selects the reported transcript against this value.
func mostSevereConsequence(anns []*annotate.Annotation) string {
	best := ""
	for _, ann := range anns {
		for rest := ann.Consequence; rest != ""; {
			term := rest
			if i := strings.IndexByte(rest, ','); i >= 0 {
				term, rest = rest[:i], rest[i+1:]
			} else {
				rest = ""
			}
			if best == "" || annotate.MoreSevereTerm(term, best) {
				best = term
			}
		}
	}
	return best
}
