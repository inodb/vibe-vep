// Package main provides the vibe-vep command-line tool.
package main

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"go.uber.org/zap"
	"go.uber.org/zap/zapcore"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/datasource/oncokb"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/output"
	"github.com/inodb/vibe-vep/internal/vcf"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

// Exit codes
const (
	ExitSuccess = 0
	ExitError   = 1
)

// Version information (set at build time)
var (
	version = "dev"
	commit  = "none"
	date    = "unknown"
)

func newRootCmd() *cobra.Command {
	var (
		verbose    bool
		configFile string
	)

	rootCmd := &cobra.Command{
		Use:   "vibe-vep",
		Short: "Variant Effect Predictor",
		Long:  "vibe-vep - Variant Effect Predictor using GENCODE annotations",
		Version: fmt.Sprintf("%s (%s) built %s", version, commit, date),
		PersistentPreRunE: func(cmd *cobra.Command, args []string) error {
			return initConfig(configFile)
		},
	}

	rootCmd.PersistentFlags().BoolVar(&verbose, "verbose", false, "Enable debug-level logging")
	rootCmd.PersistentFlags().StringVar(&configFile, "config", "", "Config file (default: $HOME/.vibe-vep.yaml)")

	rootCmd.AddCommand(newAnnotateCmd(&verbose))
	rootCmd.AddCommand(newCompareCmd(&verbose))
	rootCmd.AddCommand(newDownloadCmd(&verbose))

	return rootCmd
}

// initConfig reads in config file and ENV variables if set.
func initConfig(configFile string) error {
	if configFile != "" {
		viper.SetConfigFile(configFile)
	} else {
		home, err := os.UserHomeDir()
		if err == nil {
			viper.AddConfigPath(home)
		}
		viper.AddConfigPath(".")
		viper.SetConfigName(".vibe-vep")
		viper.SetConfigType("yaml")
	}

	viper.SetEnvPrefix("VIBE_VEP")
	viper.AutomaticEnv()
	viper.SetEnvKeyReplacer(strings.NewReplacer("-", "_"))

	// Read config file if it exists (not an error if missing)
	if err := viper.ReadInConfig(); err != nil {
		if _, ok := err.(viper.ConfigFileNotFoundError); !ok {
			return fmt.Errorf("reading config file: %w", err)
		}
	}
	return nil
}

// newLogger creates a zap logger for the CLI. In verbose mode it logs at DEBUG
// level; otherwise at INFO level.
func newLogger(verbose bool) (*zap.Logger, error) {
	cfg := zap.NewDevelopmentConfig()
	cfg.EncoderConfig.EncodeLevel = zapcore.CapitalColorLevelEncoder
	cfg.DisableStacktrace = true
	if !verbose {
		cfg.Level.SetLevel(zap.InfoLevel)
	}
	return cfg.Build()
}

func main() {
	if err := newRootCmd().Execute(); err != nil {
		os.Exit(ExitError)
	}
}

func newAnnotateCmd(verbose *bool) *cobra.Command {
	var (
		assembly      string
		outputFormat  string
		outputFile    string
		canonicalOnly bool
		inputFormat   string
	)

	cmd := &cobra.Command{
		Use:   "annotate <input-file>",
		Short: "Annotate variants in a VCF or MAF file",
		Long:  "Annotate variants in a VCF or MAF file with consequence predictions.",
		Example: `  vibe-vep annotate input.vcf
  vibe-vep annotate input.maf
  vibe-vep annotate -f vcf -o output.vcf input.vcf
  cat input.vcf | vibe-vep annotate -`,
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
			return runAnnotate(logger, args[0],
				viper.GetString("assembly"),
				viper.GetString("output-format"),
				viper.GetString("output"),
				viper.GetBool("canonical"),
				viper.GetString("input-format"),
			)
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")
	cmd.Flags().StringVarP(&outputFormat, "output-format", "f", "", "Output format: vcf, maf (default: auto-detect from input)")
	cmd.Flags().StringVarP(&outputFile, "output", "o", "", "Output file (default: stdout)")
	cmd.Flags().BoolVar(&canonicalOnly, "canonical", false, "Only report canonical transcript annotations")
	cmd.Flags().StringVar(&inputFormat, "input-format", "", "Input format: vcf, maf (auto-detected if not specified)")

	return cmd
}

func newCompareCmd(verbose *bool) *cobra.Command {
	var (
		assembly string
		columns  string
		all      bool
	)

	cmd := &cobra.Command{
		Use:   "compare <maf-file>",
		Short: "Compare MAF annotations against VEP predictions",
		Long:  "Compare MAF annotations against VEP predictions with categorized mismatch analysis.",
		Example: `  vibe-vep compare data_mutations.txt
  vibe-vep compare --columns consequence,hgvsp data_mutations.txt
  vibe-vep compare --all data_mutations.txt`,
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

			// Parse columns
			colMap := make(map[string]bool)
			for _, c := range strings.Split(viper.GetString("columns"), ",") {
				c = strings.TrimSpace(strings.ToLower(c))
				if c != "" {
					colMap[c] = true
				}
			}

			return runCompare(logger, args[0],
				viper.GetString("assembly"),
				colMap,
				viper.GetBool("all"),
			)
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")
	cmd.Flags().StringVar(&columns, "columns", "consequence,hgvsp,hgvsc", "Columns to compare (comma-separated)")
	cmd.Flags().BoolVar(&all, "all", false, "Show all rows, not just non-matches")

	return cmd
}

func runAnnotate(logger *zap.Logger, inputPath, assembly, outputFormat, outputFile string, canonicalOnly bool, inputFormat string) error {
	// Detect input format if not specified
	detectedFormat := inputFormat
	if detectedFormat == "" {
		detectedFormat = detectInputFormat(inputPath)
	}

	// Create appropriate parser
	var parser vcf.VariantParser
	var err error

	switch detectedFormat {
	case "maf":
		parser, err = maf.NewParser(inputPath)
	case "vcf":
		parser, err = vcf.NewParser(inputPath)
	default:
		return fmt.Errorf("unknown input format %q (use --input-format to specify vcf or maf)", detectedFormat)
	}

	if err != nil {
		if os.IsNotExist(err) {
			return fmt.Errorf("%w (check that the file path is correct)", err)
		}
		return err
	}
	defer parser.Close()

	// Load GENCODE cache
	gtfPath, fastaPath, canonicalPath, found := FindGENCODEFiles(assembly)
	if !found {
		return fmt.Errorf("no GENCODE cache found for %s\nHint: Download GENCODE annotations with: vibe-vep download --assembly %s", assembly, assembly)
	}

	logger.Info("using GENCODE cache",
		zap.String("assembly", assembly),
		zap.String("gtf", gtfPath),
		zap.String("fasta", fastaPath))

	c := cache.New()
	loader := cache.NewGENCODELoader(gtfPath, fastaPath)

	if canonicalPath != "" {
		logger.Info("loading canonical overrides", zap.String("path", canonicalPath))
		overrides, err := cache.LoadCanonicalOverrides(canonicalPath)
		if err != nil {
			logger.Warn("could not load canonical overrides", zap.Error(err))
		} else {
			loader.SetCanonicalOverrides(overrides)
			logger.Info("loaded canonical overrides", zap.Int("count", len(overrides)))
		}
	}

	if err := loader.Load(c); err != nil {
		return fmt.Errorf("loading GENCODE cache: %w", err)
	}
	logger.Info("loaded transcripts", zap.Int("count", c.TranscriptCount()))
	transcriptCache := c

	ann := annotate.NewAnnotator(transcriptCache)
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

	// Auto-detect output format
	if outputFormat == "" {
		if detectedFormat == "maf" {
			outputFormat = "maf"
		} else {
			outputFormat = "vcf"
		}
	}

	// Load cancer gene list from config if available
	var cgl oncokb.CancerGeneList
	if cglPath := viper.GetString("oncokb.cancer-gene-list"); cglPath != "" {
		var err error
		cgl, err = oncokb.LoadCancerGeneList(cglPath)
		if err != nil {
			logger.Warn("could not load cancer gene list", zap.String("path", cglPath), zap.Error(err))
		} else {
			logger.Info("loaded cancer gene list", zap.Int("genes", len(cgl)))
		}
	}

	// MAF output
	if outputFormat == "maf" {
		if detectedFormat != "maf" {
			return fmt.Errorf("MAF output format requires MAF input")
		}
		mafParser, ok := parser.(*maf.Parser)
		if !ok {
			return fmt.Errorf("MAF output requires MAF parser")
		}
		return runMAFOutput(logger, mafParser, ann, out, cgl)
	}

	// VCF output
	if outputFormat == "vcf" {
		if detectedFormat != "vcf" {
			return fmt.Errorf("VCF output format requires VCF input (no VCF headers to preserve)")
		}
		vcfParser, ok := parser.(*vcf.Parser)
		if !ok {
			return fmt.Errorf("VCF output requires VCF parser")
		}
		writer := output.NewVCFWriter(out, vcfParser.Header())
		if err := writer.WriteHeader(); err != nil {
			return fmt.Errorf("writing header: %w", err)
		}
		if err := ann.AnnotateAll(parser, writer); err != nil {
			return err
		}
		return nil
	}

	return fmt.Errorf("unknown output format %q", outputFormat)
}

func runCompare(logger *zap.Logger, inputPath, assembly string, columns map[string]bool, showAll bool) error {
	parser, err := maf.NewParser(inputPath)
	if err != nil {
		if os.IsNotExist(err) {
			return fmt.Errorf("%w (check that the file path is correct)", err)
		}
		return err
	}
	defer parser.Close()

	// Load GENCODE cache
	gtfPath, fastaPath, canonicalPath, found := FindGENCODEFiles(assembly)
	if !found {
		return fmt.Errorf("no GENCODE cache found for %s\nHint: Download GENCODE annotations with: vibe-vep download --assembly %s", assembly, assembly)
	}

	logger.Info("using GENCODE cache",
		zap.String("assembly", assembly),
		zap.String("gtf", gtfPath),
		zap.String("fasta", fastaPath))

	c := cache.New()
	loader := cache.NewGENCODELoader(gtfPath, fastaPath)

	if canonicalPath != "" {
		logger.Info("loading canonical overrides", zap.String("path", canonicalPath))
		overrides, err := cache.LoadCanonicalOverrides(canonicalPath)
		if err != nil {
			logger.Warn("could not load canonical overrides", zap.Error(err))
		} else {
			loader.SetCanonicalOverrides(overrides)
			logger.Info("loaded canonical overrides", zap.Int("count", len(overrides)))
		}
	}

	if err := loader.Load(c); err != nil {
		return fmt.Errorf("loading GENCODE cache: %w", err)
	}
	logger.Info("loaded transcripts", zap.Int("count", c.TranscriptCount()))

	ann := annotate.NewAnnotator(c)
	ann.SetLogger(logger)

	cmpWriter := output.NewCompareWriter(os.Stdout, columns, showAll)
	if err := cmpWriter.WriteHeader(); err != nil {
		return fmt.Errorf("writing header: %w", err)
	}

	for {
		v, mafAnn, err := parser.NextWithAnnotation()
		if err != nil {
			return fmt.Errorf("reading variant: %w", err)
		}
		if v == nil {
			break
		}

		vepAnns, err := ann.Annotate(v)
		if err != nil {
			logger.Warn("failed to annotate variant",
				zap.String("chrom", v.Chrom),
				zap.Int64("pos", v.Pos),
				zap.Error(err))
			continue
		}

		if err := cmpWriter.WriteComparison(v, mafAnn, vepAnns); err != nil {
			return fmt.Errorf("writing comparison: %w", err)
		}
	}

	if err := cmpWriter.Flush(); err != nil {
		return fmt.Errorf("flushing output: %w", err)
	}

	cmpWriter.WriteSummary(os.Stderr)

	return nil
}

// runMAFOutput runs MAF annotation mode, preserving all original columns.
func runMAFOutput(logger *zap.Logger, parser *maf.Parser, ann *annotate.Annotator, out *os.File, cgl oncokb.CancerGeneList) error {
	mafWriter := output.NewMAFWriter(out, parser.Header(), parser.Columns())

	if cgl != nil {
		mafWriter.AddExtraColumn("Gene_Type")
	}

	if err := mafWriter.WriteHeader(); err != nil {
		return fmt.Errorf("writing header: %w", err)
	}

	for {
		v, mafAnn, err := parser.NextWithAnnotation()
		if err != nil {
			return fmt.Errorf("reading variant: %w", err)
		}
		if v == nil {
			break
		}

		vepAnns, err := ann.Annotate(v)
		if err != nil {
			logger.Warn("failed to annotate variant",
				zap.String("chrom", v.Chrom),
				zap.Int64("pos", v.Pos),
				zap.Error(err))
			if writeErr := mafWriter.WriteRow(mafAnn.RawFields, nil, v); writeErr != nil {
				return fmt.Errorf("writing row: %w", writeErr)
			}
			continue
		}

		best := output.SelectBestAnnotation(mafAnn, vepAnns)
		// Enrich annotation with gene type from cancer gene list
		if best != nil && cgl != nil {
			if ga, ok := cgl[best.GeneName]; ok {
				best.GeneType = ga.GeneType
			}
		}
		if err := mafWriter.WriteRow(mafAnn.RawFields, best, v); err != nil {
			return fmt.Errorf("writing row: %w", err)
		}
	}

	return mafWriter.Flush()
}

// detectInputFormat detects the input file format based on extension or content.
func detectInputFormat(path string) string {
	lowerPath := strings.ToLower(path)

	if strings.HasSuffix(lowerPath, ".gz") {
		lowerPath = lowerPath[:len(lowerPath)-3]
	}

	if strings.HasSuffix(lowerPath, ".vcf") {
		return "vcf"
	}
	if strings.HasSuffix(lowerPath, ".maf") {
		return "maf"
	}

	baseName := filepath.Base(lowerPath)
	if baseName == "data_mutations.txt" || baseName == "data_mutations_extended.txt" {
		return "maf"
	}

	if path == "-" {
		return "vcf"
	}

	file, err := os.Open(path)
	if err != nil {
		return "vcf"
	}
	defer file.Close()

	buf := make([]byte, 512)
	n, err := file.Read(buf)
	if err != nil || n == 0 {
		return "vcf"
	}

	content := string(buf[:n])

	if strings.HasPrefix(content, "##fileformat=VCF") || strings.HasPrefix(content, "#CHROM") {
		return "vcf"
	}

	if strings.Contains(content, "Hugo_Symbol") && strings.Contains(content, "Chromosome") {
		return "maf"
	}

	return "vcf"
}
