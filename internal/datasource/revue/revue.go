// Package revue reads the reVUE catalogue of variants whose real protein
// effect is not the one VEP predicts.
//
// A handful of variants sit at a splice site and look, to any annotator
// working from the sequence alone, like a disturbed junction with no protein
// change worth naming. Published experiments show what they actually do: KIT
// 4:55593576_55593606del shortens exon 11 and deletes nine residues, which is
// why it is actionable in GIST, and p.X550_splice says none of that. reVUE
// records those findings, and genome-nexus replaces the predicted protein
// change with the revised one wherever the catalogue has an entry.
package revue

import (
	"encoding/json"
	"fmt"
	"os"
	"strings"
)

// Entry is one revised protein effect.
type Entry struct {
	// HugoGeneSymbol and TranscriptID name what the entry applies to.
	HugoGeneSymbol string
	TranscriptID   string

	// ProteinEffect is the revised protein change, in the short form, and
	// VariantClassification the MAF classification that goes with it.
	ProteinEffect         string
	VariantClassification string

	// Confirmed marks an entry supported by published evidence, as opposed to
	// one still under review.
	Confirmed bool
}

// Store answers whether a variant on a transcript has a revised effect.
type Store struct {
	byTranscriptAndLocation map[string]Entry
}

// Count returns the number of revised effects held.
func (s *Store) Count() int {
	if s == nil {
		return 0
	}
	return len(s.byTranscriptAndLocation)
}

// Lookup returns the revised effect for a variant on a transcript, if the
// catalogue has one. genomicLocation is the comma-separated
// chromosome,start,end,reference,variant form.
func (s *Store) Lookup(transcriptID, genomicLocation string) (Entry, bool) {
	if s == nil {
		return Entry{}, false
	}
	e, ok := s.byTranscriptAndLocation[key(transcriptID, genomicLocation)]
	return e, ok
}

func key(transcriptID, genomicLocation string) string {
	return transcriptID + "-" + strings.ReplaceAll(genomicLocation, " ", "")
}

// jsonGene is one gene's record in the catalogue file.
type jsonGene struct {
	HugoGeneSymbol       string       `json:"hugoGeneSymbol"`
	RevisedProteinEffect []jsonEffect `json:"revisedProteinEffects"`
}

type jsonEffect struct {
	TranscriptID string `json:"transcriptId"`
	// GenomicLocation is "chromosome,start,end,reference,variant".
	GenomicLocation string `json:"genomicLocation"`
	ProteinEffect   string `json:"revisedProteinEffect"`
	// StandardVariantClassification is the MAF spelling; VariantClassification
	// is reVUE's own, more specific name for the mechanism.
	StandardVariantClassification string `json:"revisedStandardVariantClassification"`
	VariantClassification         string `json:"revisedVariantClassification"`
	Confirmed                     bool   `json:"confirmed"`
}

// Load reads a reVUE VUEs.json catalogue.
func Load(path string) (*Store, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("open reVUE catalogue: %w", err)
	}
	return Parse(data)
}

// Parse reads a reVUE catalogue from JSON.
func Parse(data []byte) (*Store, error) {
	var genes []jsonGene
	if err := json.Unmarshal(data, &genes); err != nil {
		return nil, fmt.Errorf("parse reVUE catalogue: %w", err)
	}

	store := &Store{byTranscriptAndLocation: make(map[string]Entry)}
	for _, gene := range genes {
		for _, effect := range gene.RevisedProteinEffect {
			if effect.TranscriptID == "" || effect.GenomicLocation == "" {
				continue
			}
			class := effect.StandardVariantClassification
			if class == "" {
				class = StandardVariantClassification(effect.VariantClassification)
			}
			store.byTranscriptAndLocation[key(effect.TranscriptID, effect.GenomicLocation)] = Entry{
				HugoGeneSymbol:        gene.HugoGeneSymbol,
				TranscriptID:          effect.TranscriptID,
				ProteinEffect:         effect.ProteinEffect,
				VariantClassification: class,
				Confirmed:             effect.Confirmed,
			}
		}
	}
	return store, nil
}

// StandardVariantClassification translates reVUE's name for a mechanism into
// the MAF classification for it. A mechanism with no MAF spelling returns
// nothing, which leaves the predicted classification in place.
func StandardVariantClassification(revised string) string {
	switch revised {
	case "Splice_Exon_Skip_In_Frame", "Splice_Exon_Skip_Non_Start",
		"Splice_Exon_Shortening_In_Frame":
		return "In_Frame_Del"
	case "Splice_Exon_Skip_Out_Of_Frame", "Splice_Exon_Skip_Frame_Shift",
		"Splice_Exon_Shortening_Out_Of_Frame", "Splice_Exon_Shortening_Frame_Shift":
		return "Frame_Shift_Del"
	case "Splice_Exon_Extension_In_Frame", "Splice_Exon_Extension_Nonsense",
		"Splice_Intron_Retention_In_Frame":
		return "In_Frame_Ins"
	case "Splice_Exon_Extension_Out_Of_Frame", "Splice_Exon_Extension_Frame_Shift",
		"Splice_Intron_Retention_Out_Of_Frame":
		return "Frame_Shift_Ins"
	default:
		return ""
	}
}
