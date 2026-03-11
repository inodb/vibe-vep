package main

import (
	"fmt"
	"os"
	"runtime"
	"time"

	"go.uber.org/zap"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/output"
	"github.com/inodb/vibe-vep/internal/vcf"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

func newConvertCmd(verbose *bool) *cobra.Command {
	cmd := &cobra.Command{
		Use:   "convert",
		Short: "Convert between variant file formats",
		Long:  "Convert between variant file formats (e.g., VCF to MAF).",
	}

	cmd.AddCommand(newVCF2MAFCmd(verbose))

	return cmd
}

func newVCF2MAFCmd(verbose *bool) *cobra.Command {
	var (
		assembly      string
		outputFile    string
		canonicalOnly bool
		saveResults   bool
	)

	cmd := &cobra.Command{
		Use:   "vcf2maf <input.vcf>",
		Short: "Convert VCF to MAF format",
		Long:  "Convert a VCF file to MAF format with consequence annotations.",
		Example: `  vibe-vep convert vcf2maf input.vcf
  vibe-vep convert vcf2maf -o output.maf input.vcf
  vibe-vep convert vcf2maf --assembly GRCh37 input.vcf`,
		Args: cobra.ExactArgs(1),
		PreRunE: func(cmd *cobra.Command, args []string) error {
			return viper.BindPFlags(cmd.Flags())
		},
		RunE: func(cmd *cobra.Command, args []string) error {
			logger, err := newLogger(*verbose)
			if err != nil {
				return fmt.Errorf("creating logger: %w", err)
			}
			defer logger.Sync()
			return runConvertVCF2MAF(logger, args[0],
				viper.GetString("assembly"),
				viper.GetString("output"),
				viper.GetBool("canonical"),
				viper.GetBool("save-results"),
				viper.GetBool("no-cache"),
				viper.GetBool("clear-cache"),
			)
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")
	cmd.Flags().StringVarP(&outputFile, "output", "o", "", "Output file (default: stdout)")
	cmd.Flags().BoolVar(&canonicalOnly, "canonical", false, "Only report canonical transcript annotations")
	cmd.Flags().BoolVar(&saveResults, "save-results", false, "Save annotation results to DuckDB for later lookup")
	addCacheFlags(cmd)

	return cmd
}

func runConvertVCF2MAF(logger *zap.Logger, inputPath, assembly, outputFile string, canonicalOnly, saveResults, noCache, clearCache bool) error {
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

	// Determine tumor sample barcode
	tumorSampleID := "TUMOR"
	if names := parser.SampleNames(); len(names) > 0 {
		tumorSampleID = names[0]
	}

	writer := output.NewVCF2MAFWriter(out, assembly, tumorSampleID)
	writer.SetSources(cr.sources)
	if err := writer.WriteHeader(); err != nil {
		return fmt.Errorf("writing header: %w", err)
	}

	// Process variants
	items := make(chan annotate.WorkItem, 2*runtime.NumCPU())
	var parseErr error
	go func() {
		defer close(items)
		seq := 0
		for {
			v, err := parser.Next()
			if err != nil {
				parseErr = fmt.Errorf("reading variant: %w", err)
				return
			}
			if v == nil {
				return
			}
			// Split multi-allelic variants
			variants := vcf.SplitMultiAllelic(v)
			for _, variant := range variants {
				items <- annotate.WorkItem{Seq: seq, Variant: variant}
				seq++
			}
		}
	}()

	results := ann.ParallelAnnotate(items, 0)

	progress := func(n int) {
		logger.Info("progress", zap.Int("variants_processed", n))
	}

	if err := annotate.OrderedCollectWithProgress(results, 2*time.Second, progress, func(r annotate.WorkResult) error {
		if r.Err != nil {
			logger.Warn("failed to annotate variant",
				zap.String("chrom", r.Variant.Chrom),
				zap.Int64("pos", r.Variant.Pos),
				zap.Error(r.Err))
			return writer.WriteRow(r.Variant, nil, nil)
		}

		for _, src := range cr.sources {
			src.Annotate(r.Variant, r.Anns)
		}

		best := output.PickBestAnnotation(r.Anns)
		return writer.WriteRow(r.Variant, best, r.Anns)
	}); err != nil {
		return err
	}

	if parseErr != nil {
		return parseErr
	}

	return writer.Flush()
}
