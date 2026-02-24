// Package annotate provides variant effect prediction functionality.
package annotate

import (
	"fmt"
	"io"
	"log"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// TranscriptLookup defines the interface for finding transcripts at a position.
type TranscriptLookup interface {
	FindTranscripts(chrom string, pos int64) []*cache.Transcript
}

// Annotator annotates variants with consequence predictions.
type Annotator struct {
	cache         TranscriptLookup
	canonicalOnly bool
	warnings      io.Writer
}

// NewAnnotator creates a new annotator with the given cache.
func NewAnnotator(c TranscriptLookup) *Annotator {
	return &Annotator{
		cache: c,
	}
}

// SetCanonicalOnly configures whether to only report canonical transcript annotations.
func (a *Annotator) SetCanonicalOnly(canonical bool) {
	a.canonicalOnly = canonical
}

// SetWarnings sets the writer for warning messages.
func (a *Annotator) SetWarnings(w io.Writer) {
	a.warnings = w
}

// Annotate annotates a single variant and returns all annotations.
func (a *Annotator) Annotate(v *vcf.Variant) ([]*Annotation, error) {
	// Normalize chromosome
	chrom := v.NormalizeChrom()

	// Find overlapping transcripts
	transcripts := a.cache.FindTranscripts(chrom, v.Pos)

	if len(transcripts) == 0 {
		// Intergenic variant
		ann := &Annotation{
			VariantID:   FormatVariantID(v.Chrom, v.Pos, v.Ref, v.Alt),
			Consequence: ConsequenceIntergenicVariant,
			Impact:      GetImpact(ConsequenceIntergenicVariant),
			Allele:      v.Alt,
		}
		return []*Annotation{ann}, nil
	}

	var annotations []*Annotation

	for _, t := range transcripts {
		// Skip non-canonical if canonicalOnly is set
		if a.canonicalOnly && !t.IsCanonical {
			continue
		}

		result := PredictConsequence(v, t)

		// Append biotype-specific modifier terms per VEP convention
		consequence := result.Consequence
		if t.Biotype == "nonsense_mediated_decay" {
			consequence += ",NMD_transcript_variant"
		}

		ann := &Annotation{
			VariantID:       FormatVariantID(v.Chrom, v.Pos, v.Ref, v.Alt),
			TranscriptID:    t.ID,
			GeneName:        t.GeneName,
			GeneID:          t.GeneID,
			Consequence:     consequence,
			Impact:          result.Impact,
			CDSPosition:     result.CDSPosition,
			ProteinPosition: result.ProteinPosition,
			AminoAcidChange: result.AminoAcidChange,
			CodonChange:     result.CodonChange,
			IsCanonical:     t.IsCanonical,
			Allele:          v.Alt,
			Biotype:         t.Biotype,
			ExonNumber:      result.ExonNumber,
			IntronNumber:    result.IntronNumber,
			CDNAPosition:    result.CDNAPosition,
			HGVSp:           result.HGVSp,
		}

		annotations = append(annotations, ann)
	}

	// If no annotations after filtering, add intergenic
	if len(annotations) == 0 {
		ann := &Annotation{
			VariantID:   FormatVariantID(v.Chrom, v.Pos, v.Ref, v.Alt),
			Consequence: ConsequenceIntergenicVariant,
			Impact:      GetImpact(ConsequenceIntergenicVariant),
			Allele:      v.Alt,
		}
		return []*Annotation{ann}, nil
	}

	return annotations, nil
}

// AnnotateAll annotates all variants from a parser.
// The parser can be any type that implements vcf.VariantParser (VCF, MAF, etc.).
func (a *Annotator) AnnotateAll(parser vcf.VariantParser, writer AnnotationWriter) error {
	variantCount := 0

	for {
		v, err := parser.Next()
		if err != nil {
			return fmt.Errorf("read variant: %w", err)
		}
		if v == nil {
			break
		}

		// Split multi-allelic variants
		variants := vcf.SplitMultiAllelic(v)

		for _, variant := range variants {
			annotations, err := a.Annotate(variant)
			if err != nil {
				if a.warnings != nil {
					log.New(a.warnings, "", 0).Printf("Warning: failed to annotate %s:%d: %v\n",
						variant.Chrom, variant.Pos, err)
				}
				continue
			}

			for _, ann := range annotations {
				if err := writer.Write(variant, ann); err != nil {
					return fmt.Errorf("write annotation: %w", err)
				}
			}
		}

		variantCount++
	}

	if variantCount == 0 && a.warnings != nil {
		log.New(a.warnings, "", 0).Println("Info: 0 variants processed")
	}

	return writer.Flush()
}

// AnnotationWriter defines the interface for writing annotations.
type AnnotationWriter interface {
	WriteHeader() error
	Write(v *vcf.Variant, ann *Annotation) error
	Flush() error
}
