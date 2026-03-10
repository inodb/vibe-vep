package main

import (
	"fmt"
	"os"
	"runtime"
	"strings"
	"text/tabwriter"
	"time"

	"go.uber.org/zap"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/duckdb"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/output"
	"github.com/inodb/vibe-vep/internal/vcf"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

func newAnnotateCmd(verbose *bool) *cobra.Command {
	cmd := &cobra.Command{
		Use:   "annotate",
		Short: "Annotate variants",
		Long:  "Annotate variants in a MAF file, VCF file, or by variant specification.",
	}

	cmd.AddCommand(newAnnotateMAFCmd(verbose))
	cmd.AddCommand(newAnnotateVCFCmd(verbose))
	cmd.AddCommand(newAnnotateVariantCmd(verbose))

	return cmd
}

func newAnnotateMAFCmd(verbose *bool) *cobra.Command {
	var (
		assembly      string
		outputFile    string
		canonicalOnly bool
		saveResults   bool
		pick          bool
		mostSevere    bool
		replace       bool
	)

	cmd := &cobra.Command{
		Use:   "maf <file>",
		Short: "Annotate variants in a MAF file",
		Long: `Annotate variants in a MAF file with consequence predictions.

By default, annotations are appended as vibe.* namespaced columns.
With --replace, core columns (Hugo_Symbol, Consequence, Variant_Classification,
Transcript_ID, HGVSc, HGVSp, HGVSp_Short) are overwritten in-place.`,
		Example: `  vibe-vep annotate maf input.maf
  vibe-vep annotate maf -o output.maf input.maf
  vibe-vep annotate maf --replace -o annotated.maf input.maf
  vibe-vep annotate maf --save-results data_mutations.txt`,
		Args: cobra.ExactArgs(1),
		PreRunE: func(cmd *cobra.Command, args []string) error {
			return viper.BindPFlags(cmd.Flags())
		},
		RunE: func(cmd *cobra.Command, args []string) error {
			if viper.GetBool("pick") && viper.GetBool("most-severe") {
				return fmt.Errorf("--pick and --most-severe are mutually exclusive")
			}
			logger, err := newLogger(*verbose)
			if err != nil {
				return fmt.Errorf("creating logger: %w", err)
			}
			defer logger.Sync()
			return runAnnotateMAF(logger, args[0],
				viper.GetString("assembly"),
				viper.GetString("output"),
				viper.GetBool("canonical"),
				viper.GetBool("save-results"),
				viper.GetBool("no-cache"),
				viper.GetBool("clear-cache"),
				viper.GetBool("most-severe"),
				viper.GetBool("replace"),
			)
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")
	cmd.Flags().StringVarP(&outputFile, "output", "o", "", "Output file (default: stdout)")
	cmd.Flags().BoolVar(&canonicalOnly, "canonical", false, "Only report canonical transcript annotations")
	cmd.Flags().BoolVar(&saveResults, "save-results", false, "Save annotation results to DuckDB for later lookup")
	cmd.Flags().BoolVar(&pick, "pick", false, "One annotation per variant (best transcript)")
	cmd.Flags().BoolVar(&mostSevere, "most-severe", false, "One annotation per variant (highest impact)")
	cmd.Flags().BoolVar(&replace, "replace", false, "Overwrite core MAF columns in-place instead of appending vibe.* columns")
	addCacheFlags(cmd)

	return cmd
}

func newAnnotateVCFCmd(verbose *bool) *cobra.Command {
	var (
		assembly      string
		outputFile    string
		canonicalOnly bool
		saveResults   bool
		pick          bool
		mostSevere    bool
	)

	cmd := &cobra.Command{
		Use:   "vcf <file>",
		Short: "Annotate variants in a VCF file",
		Long:  "Annotate variants in a VCF file with consequence predictions.",
		Example: `  vibe-vep annotate vcf input.vcf
  vibe-vep annotate vcf -o output.vcf input.vcf
  vibe-vep annotate vcf --pick input.vcf
  cat input.vcf | vibe-vep annotate vcf -`,
		Args: cobra.ExactArgs(1),
		PreRunE: func(cmd *cobra.Command, args []string) error {
			return viper.BindPFlags(cmd.Flags())
		},
		RunE: func(cmd *cobra.Command, args []string) error {
			if viper.GetBool("pick") && viper.GetBool("most-severe") {
				return fmt.Errorf("--pick and --most-severe are mutually exclusive")
			}
			logger, err := newLogger(*verbose)
			if err != nil {
				return fmt.Errorf("creating logger: %w", err)
			}
			defer logger.Sync()
			return runAnnotateVCF(logger, args[0],
				viper.GetString("assembly"),
				viper.GetString("output"),
				viper.GetBool("canonical"),
				viper.GetBool("save-results"),
				viper.GetBool("no-cache"),
				viper.GetBool("clear-cache"),
				viper.GetBool("pick"),
				viper.GetBool("most-severe"),
			)
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")
	cmd.Flags().StringVarP(&outputFile, "output", "o", "", "Output file (default: stdout)")
	cmd.Flags().BoolVar(&canonicalOnly, "canonical", false, "Only report canonical transcript annotations")
	cmd.Flags().BoolVar(&saveResults, "save-results", false, "Save annotation results to DuckDB for later lookup")
	cmd.Flags().BoolVar(&pick, "pick", false, "One annotation per variant (best transcript)")
	cmd.Flags().BoolVar(&mostSevere, "most-severe", false, "One annotation per variant (highest impact)")
	addCacheFlags(cmd)

	return cmd
}

func newAnnotateVariantCmd(verbose *bool) *cobra.Command {
	var (
		assembly string
		specType string
	)

	cmd := &cobra.Command{
		Use:   "variant <spec>",
		Short: "Annotate a single variant",
		Long: `Annotate a single variant by genomic coordinates, protein change, or HGVSc notation.

Supported formats:
  Genomic:  12:25245350:C:A  or  chr12:25245350:C>A  or  12-25245350-C-A
  Protein:  KRAS G12C  or  KRAS p.G12C  or  KRAS p.Gly12Cys
  HGVSc:    KRAS c.35G>T  or  ENST00000311936:c.35G>T`,
		Example: `  vibe-vep annotate variant 12:25245350:C:A
  vibe-vep annotate variant KRAS G12C
  vibe-vep annotate variant KRAS c.35G>T
  vibe-vep annotate variant ENST00000311936:c.35G>T`,
		Args: cobra.MinimumNArgs(1),
		PreRunE: func(cmd *cobra.Command, args []string) error {
			return viper.BindPFlags(cmd.Flags())
		},
		RunE: func(cmd *cobra.Command, args []string) error {
			logger, err := newLogger(*verbose)
			if err != nil {
				return fmt.Errorf("creating logger: %w", err)
			}
			defer logger.Sync()
			return runAnnotateVariant(logger, strings.Join(args, " "),
				viper.GetString("assembly"),
				viper.GetString("type"),
				viper.GetBool("no-cache"),
				viper.GetBool("clear-cache"),
			)
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")
	cmd.Flags().StringVar(&specType, "type", "", "Force variant type: genomic, protein, or hgvsc (auto-detected if not specified)")
	addCacheFlags(cmd)

	return cmd
}

func runAnnotateMAF(logger *zap.Logger, inputPath, assembly, outputFile string, canonicalOnly, saveResults, noCache, clearCache, mostSevere, replace bool) error {
	parser, err := maf.NewParser(inputPath)
	if err != nil {
		if os.IsNotExist(err) {
			return fmt.Errorf("%w (check that the file path is correct)", err)
		}
		return err
	}
	defer parser.Close()

	cr, err := loadCache(logger, assembly, noCache, clearCache)
	if err != nil {
		return err
	}
	if cr.store != nil {
		defer cr.store.Close()
	}
	defer cr.closeSources()

	ann := annotate.NewAnnotator(cr.cache)
	ann.SetCanonicalOnly(canonicalOnly)
	ann.SetLogger(logger)

	var out *os.File
	if outputFile == "" {
		out = os.Stdout
	} else {
		out, err = os.Create(outputFile)
		if err != nil {
			return fmt.Errorf("creating output file: %w", err)
		}
		defer out.Close()
	}

	var variantResults []duckdb.VariantResult
	var collectResults *[]duckdb.VariantResult
	if saveResults && cr.store != nil {
		collectResults = &variantResults
	}

	if err := runMAFOutput(logger, parser, ann, out, cr.sources, collectResults, mostSevere, replace); err != nil {
		return err
	}

	// Write new variant results to DuckDB
	if len(variantResults) > 0 {
		start := time.Now()
		if err := cr.store.WriteVariantResults(variantResults); err != nil {
			logger.Warn("could not write variant results to cache", zap.Error(err))
		} else {
			logger.Info("wrote variant results to cache",
				zap.Int("results", len(variantResults)),
				zap.Duration("elapsed", time.Since(start)))
		}
	}
	return nil
}

func runAnnotateVCF(logger *zap.Logger, inputPath, assembly, outputFile string, canonicalOnly, saveResults, noCache, clearCache, pick, mostSevere bool) error {
	parser, err := vcf.NewParser(inputPath)
	if err != nil {
		if os.IsNotExist(err) {
			return fmt.Errorf("%w (check that the file path is correct)", err)
		}
		return err
	}
	defer parser.Close()

	cr, err := loadCache(logger, assembly, noCache, clearCache)
	if err != nil {
		return err
	}
	if cr.store != nil {
		defer cr.store.Close()
	}
	defer cr.closeSources()

	ann := annotate.NewAnnotator(cr.cache)
	ann.SetCanonicalOnly(canonicalOnly)
	ann.SetLogger(logger)

	var out *os.File
	if outputFile == "" {
		out = os.Stdout
	} else {
		out, err = os.Create(outputFile)
		if err != nil {
			return fmt.Errorf("creating output file: %w", err)
		}
		defer out.Close()
	}

	writer := output.NewVCFWriter(out, parser.Header())
	writer.SetSources(cr.sources)
	if err := writer.WriteHeader(); err != nil {
		return fmt.Errorf("writing header: %w", err)
	}

	if len(cr.sources) > 0 || pick || mostSevere {
		for {
			v, err := parser.Next()
			if err != nil {
				return fmt.Errorf("reading variant: %w", err)
			}
			if v == nil {
				break
			}
			anns, err := ann.Annotate(v)
			if err != nil {
				logger.Warn("annotation failed", zap.Error(err))
				continue
			}
			for _, src := range cr.sources {
				src.Annotate(v, anns)
			}

			// Apply pick/most-severe filtering
			if pick && len(anns) > 1 {
				anns = []*annotate.Annotation{output.PickBestAnnotation(anns)}
			} else if mostSevere && len(anns) > 1 {
				anns = []*annotate.Annotation{output.PickMostSevere(anns)}
			}

			for _, a := range anns {
				if err := writer.Write(v, a); err != nil {
					return fmt.Errorf("writing annotation: %w", err)
				}
			}
		}
		return writer.Flush()
	}

	return ann.AnnotateAll(parser, writer)
}

func runAnnotateVariant(logger *zap.Logger, specInput, assembly, specType string, noCache, clearCache bool) error {
	// Parse variant specification
	spec, err := annotate.ParseVariantSpec(specInput)
	if err != nil {
		return err
	}

	// Override type if --type flag is set
	if specType != "" {
		switch strings.ToLower(specType) {
		case "genomic":
			spec.Type = annotate.SpecGenomic
		case "protein":
			spec.Type = annotate.SpecProtein
		case "hgvsc":
			spec.Type = annotate.SpecHGVSc
		default:
			return fmt.Errorf("unknown variant type %q (use genomic, protein, or hgvsc)", specType)
		}
	}

	// Load transcript cache
	cr, err := loadCache(logger, assembly, noCache, clearCache)
	if err != nil {
		return err
	}
	if cr.store != nil {
		defer cr.store.Close()
	}
	defer cr.closeSources()

	ann := annotate.NewAnnotator(cr.cache)
	ann.SetLogger(logger)

	// Convert spec to genomic variant(s)
	var variants []*vcf.Variant
	switch spec.Type {
	case annotate.SpecGenomic:
		variants = []*vcf.Variant{{
			Chrom: spec.Chrom,
			Pos:   spec.Pos,
			Ref:   spec.Ref,
			Alt:   spec.Alt,
		}}

	case annotate.SpecProtein:
		variants, err = annotate.ReverseMapProteinChange(cr.cache, spec.GeneName, spec.RefAA, spec.Position, spec.AltAA)
		if err != nil {
			return err
		}
		fmt.Fprintf(os.Stderr, "Query: %s %c%d%c\n\n", spec.GeneName, spec.RefAA, spec.Position, spec.AltAA)
		fmt.Fprintf(os.Stderr, "Found %d genomic variant(s):\n", len(variants))
		for _, v := range variants {
			fmt.Fprintf(os.Stderr, "  %s:%d %s>%s\n", v.Chrom, v.Pos, v.Ref, v.Alt)
		}
		fmt.Fprintln(os.Stderr)

	case annotate.SpecHGVSc:
		variants, err = annotate.ReverseMapHGVSc(cr.cache, spec.TranscriptID, spec.CDSChange)
		if err != nil {
			return err
		}
		fmt.Fprintf(os.Stderr, "Query: %s c.%s\n\n", spec.TranscriptID, spec.CDSChange)
		fmt.Fprintf(os.Stderr, "Found %d genomic variant(s):\n", len(variants))
		for _, v := range variants {
			fmt.Fprintf(os.Stderr, "  %s:%d %s>%s\n", v.Chrom, v.Pos, v.Ref, v.Alt)
		}
		fmt.Fprintln(os.Stderr)
	}

	hasSources := len(cr.sources) > 0

	// Annotate each variant and display
	w := tabwriter.NewWriter(os.Stdout, 0, 4, 2, ' ', 0)

	for _, v := range variants {
		fmt.Fprintf(os.Stdout, "Variant: %s:%d %s>%s\n\n", v.Chrom, v.Pos, v.Ref, v.Alt)

		anns, err := ann.Annotate(v)
		if err != nil {
			logger.Warn("annotation failed", zap.Error(err))
			continue
		}

		for _, src := range cr.sources {
			src.Annotate(v, anns)
		}

		if hasSources {
			// Build dynamic header from sources
			header := "Gene\tTranscript\tCanon\tConsequence\tImpact\tHGVSc\tHGVSp"
			for _, src := range cr.sources {
				name := src.Name()
				for _, col := range src.Columns() {
					if name == "" {
						header += "\t" + col.Name
					} else {
						header += "\t" + name + "_" + col.Name
					}
				}
			}
			fmt.Fprintln(w, header)
			for _, a := range anns {
				canon := "no"
				if a.IsCanonical {
					canon = "yes"
				}
				line := fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s",
					a.GeneName, a.TranscriptID, canon,
					a.Consequence, a.Impact, a.HGVSc, a.HGVSp)
				for _, src := range cr.sources {
					name := src.Name()
					for _, col := range src.Columns() {
						if name == "" {
							line += "\t" + a.GetExtraKey(col.Name)
						} else {
							line += "\t" + a.GetExtra(name, col.Name)
						}
					}
				}
				fmt.Fprintln(w, line)
			}
		} else {
			fmt.Fprintln(w, "Gene\tTranscript\tCanon\tConsequence\tImpact\tHGVSc\tHGVSp")
			for _, a := range anns {
				canon := "no"
				if a.IsCanonical {
					canon = "yes"
				}
				fmt.Fprintf(w, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
					a.GeneName, a.TranscriptID, canon,
					a.Consequence, a.Impact, a.HGVSc, a.HGVSp)
			}
		}
		w.Flush()
	}

	// DuckDB lookup for previously seen results
	if cr.store != nil && len(variants) > 0 {
		for _, v := range variants {
			chrom := v.NormalizeChrom()
			prev, err := cr.store.LookupVariant(chrom, v.Pos, v.Ref, v.Alt)
			if err != nil {
				logger.Debug("DuckDB lookup failed", zap.Error(err))
				continue
			}
			if len(prev) > 0 {
				fmt.Fprintf(os.Stdout, "\nPreviously seen (%d cached annotations):\n", len(prev))
				pw := tabwriter.NewWriter(os.Stdout, 0, 4, 2, ' ', 0)
				fmt.Fprintln(pw, "Gene\tTranscript\tCanon\tConsequence\tImpact\tHGVSc\tHGVSp")
				for _, a := range prev {
					canon := "no"
					if a.IsCanonical {
						canon = "yes"
					}
					fmt.Fprintf(pw, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
						a.GeneName, a.TranscriptID, canon,
						a.Consequence, a.Impact, a.HGVSc, a.HGVSp)
				}
				pw.Flush()
			}
		}
	}

	return nil
}

// runMAFOutput runs MAF annotation mode, preserving all original columns.
func runMAFOutput(logger *zap.Logger, parser *maf.Parser, ann *annotate.Annotator, out *os.File, sources []annotate.AnnotationSource, newResults *[]duckdb.VariantResult, mostSevere, replace bool) error {
	mafWriter := output.NewMAFWriter(out, parser.Header(), parser.Columns())
	mafWriter.SetSources(sources)
	mafWriter.SetReplace(replace)

	if err := mafWriter.WriteHeader(); err != nil {
		return fmt.Errorf("writing header: %w", err)
	}

	// Parse variants in a goroutine, send to worker pool.
	items := make(chan annotate.WorkItem, 2*runtime.NumCPU())
	var parseErr error
	go func() {
		defer close(items)
		seq := 0
		for {
			v, mafAnn, err := parser.NextWithAnnotation()
			if err != nil {
				parseErr = fmt.Errorf("reading variant: %w", err)
				return
			}
			if v == nil {
				return
			}
			items <- annotate.WorkItem{Seq: seq, Variant: v, Extra: mafAnn}
			seq++
		}
	}()

	results := ann.ParallelAnnotate(items, 0)

	progress := func(n int) {
		logger.Info("progress", zap.Int("variants_processed", n))
	}

	if err := annotate.OrderedCollectWithProgress(results, 2*time.Second, progress, func(r annotate.WorkResult) error {
		mafAnn := r.Extra.(*maf.MAFAnnotation)
		if r.Err != nil {
			logger.Warn("failed to annotate variant",
				zap.String("chrom", r.Variant.Chrom),
				zap.Int64("pos", r.Variant.Pos),
				zap.Error(r.Err))
			return mafWriter.WriteRow(mafAnn.RawFields, nil, r.Variant)
		}

		// Collect all results for DuckDB persistence
		if newResults != nil && len(r.Anns) > 0 {
			chrom := r.Variant.NormalizeChrom()
			for _, a := range r.Anns {
				*newResults = append(*newResults, duckdb.VariantResult{
					Chrom: chrom, Pos: r.Variant.Pos, Ref: r.Variant.Ref, Alt: r.Variant.Alt, Ann: a,
				})
			}
		}

		var best *annotate.Annotation
		if mostSevere {
			best = output.PickMostSevere(r.Anns)
		} else {
			best = output.SelectBestAnnotation(mafAnn, r.Anns)
		}
		// Enrich best annotation with annotation sources
		if best != nil {
			for _, src := range sources {
				src.Annotate(r.Variant, []*annotate.Annotation{best})
			}
		}
		return mafWriter.WriteRow(mafAnn.RawFields, best, r.Variant)
	}); err != nil {
		return err
	}

	if parseErr != nil {
		return parseErr
	}

	return mafWriter.Flush()
}
