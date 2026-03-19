package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"

	"go.uber.org/zap"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/input"
	"github.com/inodb/vibe-vep/internal/output"
	"github.com/inodb/vibe-vep/internal/vcf"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

func newAnnotateStreamCmd(verbose *bool) *cobra.Command {
	var (
		assembly     string
		inputFormat  string
		outputFormat string
	)

	cmd := &cobra.Command{
		Use:   "stream",
		Short: "Stream-annotate variants from stdin",
		Long: `Read variants from stdin (one JSON per line) and write annotations to stdout.

Designed for use as a long-running subprocess — loads transcript cache and
annotation sources once, then processes variants as they arrive. Each input
line produces exactly one output line.

Supported input formats:
  genome-nexus-genomic-location-jsonl  Genome Nexus GenomicLocation JSON
    {"chromosome":"7","start":140453136,"end":140453136,"referenceAllele":"A","variantAllele":"T"}

Supported output formats:
  ensembl-vep-jsonl   VEP-compatible JSON (default)
  vibe-vep-jsonl      vibe-vep native JSON with all annotation source extras`,
		Example: `  # Stream annotation (default: genome-nexus input, VEP output)
  echo '{"chromosome":"12","start":25245350,"end":25245350,"referenceAllele":"C","variantAllele":"A"}' \
    | vibe-vep annotate stream

  # Use vibe-vep native output format
  echo '{"chromosome":"7","start":140453136,"end":140453136,"referenceAllele":"A","variantAllele":"T"}' \
    | vibe-vep annotate stream --output-format vibe-vep-jsonl`,
		Args: cobra.NoArgs,
		PreRunE: func(cmd *cobra.Command, args []string) error {
			return viper.BindPFlags(cmd.Flags())
		},
		RunE: func(cmd *cobra.Command, args []string) error {
			logger, err := newLogger(*verbose)
			if err != nil {
				return fmt.Errorf("creating logger: %w", err)
			}
			defer logger.Sync()
			return runAnnotateStream(logger,
				viper.GetString("assembly"),
				viper.GetString("input-format"),
				viper.GetString("output-format"),
				viper.GetBool("no-cache"),
				viper.GetBool("clear-cache"),
			)
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")
	cmd.Flags().StringVar(&inputFormat, "input-format", "genome-nexus-genomic-location-jsonl", "Input format")
	cmd.Flags().StringVar(&outputFormat, "output-format", "ensembl-vep-jsonl", "Output format: ensembl-vep-jsonl or vibe-vep-jsonl")
	addCacheFlags(cmd)

	return cmd
}

func runAnnotateStream(logger *zap.Logger, assembly, inputFmt, outputFmt string, noCache, clearCache bool) error {
	// Validate formats.
	switch inputFmt {
	case "genome-nexus-genomic-location-jsonl":
		// ok
	default:
		return fmt.Errorf("unsupported input format %q (use: genome-nexus-genomic-location-jsonl)", inputFmt)
	}
	switch outputFmt {
	case "ensembl-vep-jsonl", "vibe-vep-jsonl":
		// ok
	default:
		return fmt.Errorf("unsupported output format %q (use: ensembl-vep-jsonl or vibe-vep-jsonl)", outputFmt)
	}

	assembly, err := normalizeAssembly(assembly)
	if err != nil {
		return err
	}

	// Load cache + annotation sources once.
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

	writer := output.NewJSONLWriter(os.Stdout, outputFmt, assembly)

	logger.Info("stream mode ready",
		zap.String("input_format", inputFmt),
		zap.String("output_format", outputFmt),
		zap.String("assembly", assembly),
		zap.Int("sources", len(cr.sources)))

	// Read stdin line by line.
	scanner := bufio.NewScanner(os.Stdin)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	for scanner.Scan() {
		line := scanner.Bytes()
		if len(line) == 0 {
			continue
		}

		// Trim whitespace.
		lineStr := strings.TrimSpace(string(line))
		if lineStr == "" {
			continue
		}

		// Resolve input to genomic variant(s).
		var variants []*vcf.Variant
		var inputLabel string

		if lineStr[0] == '{' {
			// JSON input: GenomicLocation
			gl, err := input.ParseGenomicLocation([]byte(lineStr))
			if err != nil {
				logger.Warn("parse error", zap.Error(err), zap.String("input", lineStr))
				writer.WriteError(lineStr, err.Error())
				continue
			}
			variants = []*vcf.Variant{gl.ToVariant()}
			inputLabel = gl.FormatInput()
		} else {
			// Plain string: try as variant spec (genomic coords, HGVSc, HGVSg, protein)
			spec, err := annotate.ParseVariantSpec(lineStr)
			if err != nil {
				logger.Warn("parse error", zap.Error(err), zap.String("input", lineStr))
				writer.WriteError(lineStr, err.Error())
				continue
			}
			inputLabel = lineStr

			switch spec.Type {
			case annotate.SpecGenomic:
				variants = []*vcf.Variant{{Chrom: spec.Chrom, Pos: spec.Pos, Ref: spec.Ref, Alt: spec.Alt}}
			case annotate.SpecHGVSc:
				variants, err = annotate.ReverseMapHGVSc(cr.cache, spec.TranscriptID, spec.CDSChange)
			case annotate.SpecHGVSg:
				variants, err = annotate.ResolveHGVSg(cr.cache, spec.Chrom, spec.GenomicChange)
			case annotate.SpecProtein:
				variants, err = annotate.ReverseMapProteinChange(cr.cache, spec.GeneName, spec.RefAA, spec.Position, spec.AltAA)
			}
			if err != nil {
				logger.Warn("resolve error", zap.Error(err), zap.String("input", lineStr))
				writer.WriteError(lineStr, err.Error())
				continue
			}
		}

		for _, v := range variants {
			writer.SetInput(inputLabel)

			anns, err := ann.Annotate(v)
			if err != nil {
				logger.Warn("annotation error", zap.Error(err), zap.String("input", lineStr))
				writer.WriteError(lineStr, err.Error())
				continue
			}

			// Apply annotation sources.
			for _, src := range cr.sources {
				src.Annotate(v, anns)
			}

			// Write all annotations for this variant.
			for _, a := range anns {
				if err := writer.Write(v, a); err != nil {
					return fmt.Errorf("write annotation: %w", err)
				}
			}
			if err := writer.Flush(); err != nil {
				return fmt.Errorf("flush: %w", err)
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return fmt.Errorf("read stdin: %w", err)
	}

	return nil
}
