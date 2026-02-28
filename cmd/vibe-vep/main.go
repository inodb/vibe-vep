// Package main provides the vibe-vep command-line tool.
package main

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"time"

	"go.uber.org/zap"
	"go.uber.org/zap/zapcore"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/datasource/oncokb"
	"github.com/inodb/vibe-vep/internal/duckdb"
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
	rootCmd.AddCommand(newPrepareCmd(&verbose))

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

func addCacheFlags(cmd *cobra.Command) {
	cmd.Flags().Bool("no-cache", false, "Skip transcript cache, always load from GTF/FASTA")
	cmd.Flags().Bool("clear-cache", false, "Clear and rebuild transcript and variant caches")
	cmd.Flags().Bool("no-variant-cache", false, "Disable variant result caching")
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
				viper.GetBool("no-cache"),
				viper.GetBool("clear-cache"),
				viper.GetBool("no-variant-cache"),
			)
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")
	cmd.Flags().StringVarP(&outputFormat, "output-format", "f", "", "Output format: vcf, maf (default: auto-detect from input)")
	cmd.Flags().StringVarP(&outputFile, "output", "o", "", "Output file (default: stdout)")
	cmd.Flags().BoolVar(&canonicalOnly, "canonical", false, "Only report canonical transcript annotations")
	cmd.Flags().StringVar(&inputFormat, "input-format", "", "Input format: vcf, maf (auto-detected if not specified)")
	addCacheFlags(cmd)

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
				viper.GetBool("no-cache"),
				viper.GetBool("clear-cache"),
			)
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")
	cmd.Flags().StringVar(&columns, "columns", "consequence,hgvsp,hgvsc", "Columns to compare (comma-separated)")
	cmd.Flags().BoolVar(&all, "all", false, "Show all rows, not just non-matches")
	addCacheFlags(cmd)

	return cmd
}

// cacheResult holds the loaded transcript cache and optional DuckDB variant store.
type cacheResult struct {
	cache *cache.Cache
	store *duckdb.Store // variant cache (DuckDB), nil if --no-cache
}

// loadCache loads transcripts using gob transcript cache, and opens DuckDB for variant cache.
func loadCache(logger *zap.Logger, assembly string, noCache, clearCache bool) (*cacheResult, error) {
	gtfPath, fastaPath, canonicalPath, found := FindGENCODEFiles(assembly)
	if !found {
		return nil, fmt.Errorf("no GENCODE cache found for %s\nHint: Download GENCODE annotations with: vibe-vep download --assembly %s", assembly, assembly)
	}

	logger.Info("using GENCODE cache",
		zap.String("assembly", assembly),
		zap.String("gtf", gtfPath),
		zap.String("fasta", fastaPath))

	c := cache.New()
	cacheDir := DefaultGENCODEPath(assembly)

	// Fingerprint source files for cache validation
	gtfFP, err1 := duckdb.StatFile(gtfPath)
	fastaFP, err2 := duckdb.StatFile(fastaPath)
	canonicalFP := duckdb.FileFingerprint{}
	if canonicalPath != "" {
		canonicalFP, _ = duckdb.StatFile(canonicalPath)
	}

	// --- Transcript cache (gob) ---
	transcriptsLoaded := false
	tc := duckdb.NewTranscriptCache(cacheDir)

	if noCache || clearCache {
		if clearCache {
			tc.Clear()
			logger.Info("cleared transcript cache")
		}
	} else if err1 == nil && err2 == nil && tc.Valid(gtfFP, fastaFP, canonicalFP) {
		start := time.Now()
		if err := tc.Load(c); err != nil {
			logger.Warn("transcript cache load failed, falling back to GTF/FASTA", zap.Error(err))
		} else {
			logger.Info("loaded transcript cache",
				zap.Int("count", c.TranscriptCount()),
				zap.Duration("elapsed", time.Since(start)))
			transcriptsLoaded = true
		}
	}

	if !transcriptsLoaded {
		// Load from GTF/FASTA
		if err := loadFromGTFFASTA(logger, c, gtfPath, fastaPath, canonicalPath); err != nil {
			return nil, err
		}

		// Write transcript cache for next time
		if !noCache && err1 == nil && err2 == nil {
			start := time.Now()
			if err := tc.Write(c, gtfFP, fastaFP, canonicalFP); err != nil {
				logger.Warn("could not write transcript cache", zap.Error(err))
			} else {
				logger.Info("wrote transcript cache",
					zap.Int("count", c.TranscriptCount()),
					zap.Duration("elapsed", time.Since(start)))
			}
		}
	}

	// --- Variant cache (DuckDB) ---
	if noCache {
		return &cacheResult{cache: c}, nil
	}

	dbPath := filepath.Join(cacheDir, "variant_cache.duckdb")
	store, err := duckdb.Open(dbPath)
	if err != nil {
		logger.Warn("could not open variant cache", zap.Error(err))
		return &cacheResult{cache: c}, nil
	}

	// Clear variant cache when transcripts changed (annotations depend on transcript data)
	if clearCache || !transcriptsLoaded {
		if err := store.ClearVariantResults(); err != nil {
			logger.Warn("could not clear variant cache", zap.Error(err))
		} else if clearCache {
			logger.Info("cleared variant cache")
		}
	}

	return &cacheResult{cache: c, store: store}, nil
}

// loadFromGTFFASTA loads transcripts from GENCODE GTF and FASTA files.
func loadFromGTFFASTA(logger *zap.Logger, c *cache.Cache, gtfPath, fastaPath, canonicalPath string) error {
	start := time.Now()
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
	logger.Info("loaded transcripts from GTF/FASTA",
		zap.Int("count", c.TranscriptCount()),
		zap.Duration("elapsed", time.Since(start)))
	return nil
}

func newPrepareCmd(verbose *bool) *cobra.Command {
	var assembly string

	cmd := &cobra.Command{
		Use:   "prepare",
		Short: "Build transcript cache for fast startup",
		Long:  "Load GENCODE GTF/FASTA and build the transcript cache so subsequent annotate/compare runs start instantly.",
		Example: `  vibe-vep prepare
  vibe-vep prepare --assembly GRCh37`,
		PreRunE: func(cmd *cobra.Command, args []string) error {
			return viper.BindPFlags(cmd.Flags())
		},
		RunE: func(cmd *cobra.Command, args []string) error {
			logger, err := newLogger(*verbose)
			if err != nil {
				return fmt.Errorf("creating logger: %w", err)
			}
			defer logger.Sync()

			cr, err := loadCache(logger, viper.GetString("assembly"), false, true)
			if err != nil {
				return err
			}
			if cr.store != nil {
				cr.store.Close()
			}
			logger.Info("transcript cache ready",
				zap.Int("transcripts", cr.cache.TranscriptCount()))
			return nil
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")

	return cmd
}

func runAnnotate(logger *zap.Logger, inputPath, assembly, outputFormat, outputFile string, canonicalOnly bool, inputFormat string, noCache, clearCache, noVariantCache bool) error {
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

	// Load transcript cache
	cr, err := loadCache(logger, assembly, noCache, clearCache)
	if err != nil {
		return err
	}
	if cr.store != nil {
		defer cr.store.Close()
	}

	ann := annotate.NewAnnotator(cr.cache)
	ann.SetCanonicalOnly(canonicalOnly)
	ann.SetLogger(logger)

	// Load variant cache
	var variantResults []duckdb.VariantResult
	if cr.store != nil && !noVariantCache {
		start := time.Now()
		vc, err := cr.store.LoadVariantCache()
		if err != nil {
			logger.Warn("could not load variant cache", zap.Error(err))
		} else if vc.Len() > 0 {
			ann.SetVariantCache(vc)
			logger.Info("loaded variant cache",
				zap.Int("variants", vc.Len()),
				zap.Duration("elapsed", time.Since(start)))
		}
	}

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
		if err := runMAFOutput(logger, mafParser, ann, out, cgl); err != nil {
			return err
		}

		// Write new variant results to DuckDB
		if cr.store != nil && !noVariantCache && len(variantResults) > 0 {
			if err := cr.store.WriteVariantResults(variantResults); err != nil {
				logger.Warn("could not write variant results to cache", zap.Error(err))
			}
		}
		return nil
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

func runCompare(logger *zap.Logger, inputPath, assembly string, columns map[string]bool, showAll bool, noCache, clearCache bool) error {
	parser, err := maf.NewParser(inputPath)
	if err != nil {
		if os.IsNotExist(err) {
			return fmt.Errorf("%w (check that the file path is correct)", err)
		}
		return err
	}
	defer parser.Close()

	// Load transcript cache
	cr, err := loadCache(logger, assembly, noCache, clearCache)
	if err != nil {
		return err
	}
	if cr.store != nil {
		defer cr.store.Close()
	}

	ann := annotate.NewAnnotator(cr.cache)
	ann.SetLogger(logger)

	cmpWriter := output.NewCompareWriter(os.Stdout, columns, showAll)
	if err := cmpWriter.WriteHeader(); err != nil {
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

	if err := annotate.OrderedCollect(results, func(r annotate.WorkResult) error {
		mafAnn := r.Extra.(*maf.MAFAnnotation)
		if r.Err != nil {
			logger.Warn("failed to annotate variant",
				zap.String("chrom", r.Variant.Chrom),
				zap.Int64("pos", r.Variant.Pos),
				zap.Error(r.Err))
			return nil
		}
		return cmpWriter.WriteComparison(r.Variant, mafAnn, r.Anns)
	}); err != nil {
		return err
	}

	if parseErr != nil {
		return parseErr
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

	if err := annotate.OrderedCollect(results, func(r annotate.WorkResult) error {
		mafAnn := r.Extra.(*maf.MAFAnnotation)
		if r.Err != nil {
			logger.Warn("failed to annotate variant",
				zap.String("chrom", r.Variant.Chrom),
				zap.Int64("pos", r.Variant.Pos),
				zap.Error(r.Err))
			return mafWriter.WriteRow(mafAnn.RawFields, nil, r.Variant)
		}

		best := output.SelectBestAnnotation(mafAnn, r.Anns)
		if best != nil && cgl != nil {
			if ga, ok := cgl[best.GeneName]; ok {
				best.GeneType = ga.GeneType
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
