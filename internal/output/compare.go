// Package output provides output formatting for annotations.
package output

import (
	"fmt"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/maf"
)

// transcriptBaseID strips the version suffix (e.g. ".4") from a transcript ID.
func transcriptBaseID(id string) string {
	if idx := strings.LastIndexByte(id, '.'); idx >= 0 {
		return id[:idx]
	}
	return id
}

// SelectBestAnnotation picks the best VEP annotation to compare against a MAF entry.
func SelectBestAnnotation(mafAnn *maf.MAFAnnotation, vepAnns []*annotate.Annotation) *annotate.Annotation {
	var bestAnn *annotate.Annotation

	// Pass 1: transcript ID match ignoring version suffix.
	if mafAnn.TranscriptID != "" {
		mafBase := transcriptBaseID(mafAnn.TranscriptID)
		for _, ann := range vepAnns {
			if transcriptBaseID(ann.TranscriptID) == mafBase {
				if isCodingConsequence(mafAnn.Consequence) && !isProteinCodingBiotype(ann.Biotype) {
					break
				}
				return ann
			}
		}
	}

	// Pass 2: same gene, prefer protein-coding > canonical > higher impact
	var sameGene *annotate.Annotation
	if mafAnn.HugoSymbol != "" {
		for _, ann := range vepAnns {
			if ann.GeneName == mafAnn.HugoSymbol {
				if sameGene == nil || AnnotationBetter(ann, sameGene) {
					sameGene = ann
				}
			}
		}
	}
	if sameGene != nil {
		return sameGene
	}

	// Pass 3: prefer protein-coding > canonical > higher impact
	for _, ann := range vepAnns {
		if bestAnn == nil || AnnotationBetter(ann, bestAnn) {
			bestAnn = ann
		}
	}
	return bestAnn
}

// PickBestAnnotation selects the best annotation without MAF context.
func PickBestAnnotation(anns []*annotate.Annotation) *annotate.Annotation {
	if len(anns) == 0 {
		return nil
	}
	best := anns[0]
	for _, ann := range anns[1:] {
		if AnnotationBetter(ann, best) {
			best = ann
		}
	}
	return best
}

// PickMostSevere selects the highest-impact annotation.
func PickMostSevere(anns []*annotate.Annotation) *annotate.Annotation {
	if len(anns) == 0 {
		return nil
	}
	best := anns[0]
	for _, ann := range anns[1:] {
		annImpact := annotate.ImpactRank(ann.Impact)
		bestImpact := annotate.ImpactRank(best.Impact)
		if annImpact > bestImpact {
			best = ann
		} else if annImpact == bestImpact {
			if ann.IsCanonicalMSK && !best.IsCanonicalMSK {
				best = ann
			} else if ann.IsCanonicalMSK == best.IsCanonicalMSK {
				if isProteinCodingBiotype(ann.Biotype) && !isProteinCodingBiotype(best.Biotype) {
					best = ann
				}
			}
		}
	}
	return best
}

// AnnotationBetter returns true if ann is a better pick than current for comparison.
// Priority order: protein-coding biotype > canonical > impact > has HGVSp.
func AnnotationBetter(ann, current *annotate.Annotation) bool {
	annCoding := isProteinCodingBiotype(ann.Biotype)
	curCoding := isProteinCodingBiotype(current.Biotype)
	if annCoding != curCoding {
		return annCoding
	}
	if ann.IsCanonicalMSK != current.IsCanonicalMSK {
		return ann.IsCanonicalMSK
	}
	annImpact := annotate.ImpactRank(ann.Impact)
	curImpact := annotate.ImpactRank(current.Impact)
	if annImpact != curImpact {
		return annImpact > curImpact
	}
	if (ann.HGVSp != "") != (current.HGVSp != "") {
		return ann.HGVSp != ""
	}
	return false
}

// isProteinCodingBiotype returns true if the biotype has coding potential.
func isProteinCodingBiotype(biotype string) bool {
	switch biotype {
	case "protein_coding", "nonsense_mediated_decay", "non_stop_decay",
		"IG_V_gene", "IG_D_gene", "IG_J_gene", "IG_C_gene",
		"TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene",
		"protein_coding_LoF":
		return true
	}
	return false
}

// FormatAllEffects packs all annotations into a semicolon-delimited string.
// Format per entry: gene,consequence,HGVSp_Short,transcript,HGVSc,impact,canonical_msk
func FormatAllEffects(anns []*annotate.Annotation) string {
	if len(anns) == 0 {
		return ""
	}
	var b strings.Builder
	for i, ann := range anns {
		if i > 0 {
			b.WriteByte(';')
		}
		canonMSK := ""
		if ann.IsCanonicalMSK {
			canonMSK = "YES"
		}
		b.WriteString(ann.GeneName)
		b.WriteByte(',')
		b.WriteString(ann.Consequence)
		b.WriteByte(',')
		b.WriteString(HGVSpToShort(ann.HGVSp))
		b.WriteByte(',')
		b.WriteString(ann.TranscriptID)
		b.WriteByte(',')
		b.WriteString(ann.HGVSc)
		b.WriteByte(',')
		b.WriteString(ann.Impact)
		b.WriteByte(',')
		b.WriteString(canonMSK)
	}
	return b.String()
}

// ValidOutputColumns returns the set of excludable output column names.
func ValidOutputColumns() map[string]bool {
	m := make(map[string]bool, len(annotate.CoreColumns)+1)
	for _, col := range annotate.CoreColumns {
		m[col.Name] = true
	}
	m["all_effects"] = true
	return m
}

// ValidateExcludeColumns returns an error if any column name is not a valid excludable column.
func ValidateExcludeColumns(cols []string) error {
	valid := ValidOutputColumns()
	for _, col := range cols {
		if !valid[col] {
			var names []string
			for _, c := range annotate.CoreColumns {
				names = append(names, c.Name)
			}
			names = append(names, "all_effects")
			return fmt.Errorf("unknown column %q; valid columns: %s", col, strings.Join(names, ", "))
		}
	}
	return nil
}

// HGVSpToShort converts 3-letter HGVSp notation to single-letter.
func HGVSpToShort(hgvsp string) string {
	if len(hgvsp) == 0 {
		return hgvsp
	}
	var b strings.Builder
	b.Grow(len(hgvsp))
	for i := 0; i < len(hgvsp); {
		if i+3 <= len(hgvsp) {
			if single, ok := annotate.AminoAcidThreeToSingle[hgvsp[i:i+3]]; ok {
				b.WriteByte(single)
				i += 3
				continue
			}
		}
		b.WriteByte(hgvsp[i])
		i++
	}
	return b.String()
}
