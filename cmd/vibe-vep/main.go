// Package main provides the vibe-vep command-line tool.
package main

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"text/tabwriter"
	"time"

	"go.uber.org/zap"
	"go.uber.org/zap/zapcore"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/datasource/alphamissense"
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
	rootCmd.AddCommand(newConfigCmd())
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
}

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
	)

	cmd := &cobra.Command{
		Use:   "maf <file>",
		Short: "Annotate variants in a MAF file",
		Long:  "Annotate variants in a MAF file with consequence predictions.",
		Example: `  vibe-vep annotate maf input.maf
  vibe-vep annotate maf -o output.maf input.maf
  vibe-vep annotate maf --save-results data_mutations.txt`,
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
			return runAnnotateMAF(logger, args[0],
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

func newAnnotateVCFCmd(verbose *bool) *cobra.Command {
	var (
		assembly      string
		outputFile    string
		canonicalOnly bool
		saveResults   bool
	)

	cmd := &cobra.Command{
		Use:   "vcf <file>",
		Short: "Annotate variants in a VCF file",
		Long:  "Annotate variants in a VCF file with consequence predictions.",
		Example: `  vibe-vep annotate vcf input.vcf
  vibe-vep annotate vcf -o output.vcf input.vcf
  cat input.vcf | vibe-vep annotate vcf -`,
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
			return runAnnotateVCF(logger, args[0],
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
	cache         *cache.Cache
	store         *duckdb.Store          // variant cache (DuckDB), nil if --no-cache
	alphaMissense *alphamissense.Store   // AlphaMissense scores, nil if not enabled
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

	cr := &cacheResult{cache: c, store: store}

	// --- AlphaMissense (optional) ---
	if viper.GetBool("annotations.alphamissense") {
		amStore, err := loadAlphaMissense(logger, cacheDir, assembly)
		if err != nil {
			logger.Warn("could not load AlphaMissense data", zap.Error(err))
		} else {
			cr.alphaMissense = amStore
		}
	}

	return cr, nil
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

// loadAlphaMissense opens the AlphaMissense DuckDB store, loading from TSV if needed.
func loadAlphaMissense(logger *zap.Logger, cacheDir, assembly string) (*alphamissense.Store, error) {
	dbPath := filepath.Join(cacheDir, "annotations.duckdb")
	amStore, err := alphamissense.Open(dbPath)
	if err != nil {
		return nil, fmt.Errorf("open AlphaMissense store: %w", err)
	}

	if !amStore.Loaded() {
		tsvPath := filepath.Join(cacheDir, AlphaMissenseFileName(assembly))
		if _, err := os.Stat(tsvPath); err != nil {
			amStore.Close()
			return nil, fmt.Errorf("AlphaMissense TSV not found at %s\nHint: run 'vibe-vep download' to fetch it", tsvPath)
		}

		logger.Info("loading AlphaMissense data (this may take a few minutes)...", zap.String("path", tsvPath))
		start := time.Now()
		if err := amStore.Load(tsvPath); err != nil {
			amStore.Close()
			return nil, fmt.Errorf("load AlphaMissense data: %w", err)
		}
		logger.Info("loaded AlphaMissense data", zap.Duration("elapsed", time.Since(start)))
	} else {
		logger.Info("AlphaMissense data already loaded")
	}

	return amStore, nil
}

// isMissense returns true if the consequence includes missense_variant.
func isMissense(consequence string) bool {
	for rest := consequence; rest != ""; {
		term := rest
		if i := strings.IndexByte(rest, ','); i >= 0 {
			term = rest[:i]
			rest = rest[i+1:]
		} else {
			rest = ""
		}
		if term == "missense_variant" {
			return true
		}
	}
	return false
}

// enrichAlphaMissense adds AlphaMissense scores to missense annotations.
func enrichAlphaMissense(amStore *alphamissense.Store, v *vcf.Variant, anns []*annotate.Annotation) {
	chrom := v.NormalizeChrom()
	// Add "chr" prefix if not present (AlphaMissense uses "chr" prefixed chromosomes)
	if len(chrom) > 0 && chrom[0] != 'c' {
		chrom = "chr" + chrom
	}
	for _, ann := range anns {
		if isMissense(ann.Consequence) {
			if r, ok := amStore.Lookup(chrom, v.Pos, v.Ref, v.Alt); ok {
				ann.AlphaMissenseScore = r.Score
				ann.AlphaMissenseClass = r.Class
			}
		}
	}
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
			if cr.alphaMissense != nil {
				cr.alphaMissense.Close()
				logger.Info("AlphaMissense data ready")
			}
			logger.Info("transcript cache ready",
				zap.Int("transcripts", cr.cache.TranscriptCount()))
			return nil
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")

	return cmd
}

func runAnnotateMAF(logger *zap.Logger, inputPath, assembly, outputFile string, canonicalOnly, saveResults, noCache, clearCache bool) error {
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
	if cr.alphaMissense != nil {
		defer cr.alphaMissense.Close()
	}

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

	// Load cancer gene list from config if available
	var cgl oncokb.CancerGeneList
	if cglPath := viper.GetString("oncokb.cancer-gene-list"); cglPath != "" {
		cgl, err = oncokb.LoadCancerGeneList(cglPath)
		if err != nil {
			logger.Warn("could not load cancer gene list", zap.String("path", cglPath), zap.Error(err))
		} else {
			logger.Info("loaded cancer gene list", zap.Int("genes", len(cgl)))
		}
	}

	var variantResults []duckdb.VariantResult
	var collectResults *[]duckdb.VariantResult
	if saveResults && cr.store != nil {
		collectResults = &variantResults
	}

	if err := runMAFOutput(logger, parser, ann, out, cgl, collectResults, cr.alphaMissense); err != nil {
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

func runAnnotateVCF(logger *zap.Logger, inputPath, assembly, outputFile string, canonicalOnly, saveResults, noCache, clearCache bool) error {
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
	if cr.alphaMissense != nil {
		defer cr.alphaMissense.Close()
	}

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
	if cr.alphaMissense != nil {
		writer.AddCSQField("AlphaMissense_score")
		writer.AddCSQField("AlphaMissense_class")
	}
	if err := writer.WriteHeader(); err != nil {
		return fmt.Errorf("writing header: %w", err)
	}

	// If AlphaMissense is enabled, use enriching annotation pipeline
	if cr.alphaMissense != nil {
		amStore := cr.alphaMissense
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
			enrichAlphaMissense(amStore, v, anns)
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
	if cr.alphaMissense != nil {
		defer cr.alphaMissense.Close()
	}

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

	// Annotate each variant and display
	w := tabwriter.NewWriter(os.Stdout, 0, 4, 2, ' ', 0)

	for _, v := range variants {
		fmt.Fprintf(os.Stdout, "Variant: %s:%d %s>%s\n\n", v.Chrom, v.Pos, v.Ref, v.Alt)

		anns, err := ann.Annotate(v)
		if err != nil {
			logger.Warn("annotation failed", zap.Error(err))
			continue
		}

		if cr.alphaMissense != nil {
			enrichAlphaMissense(cr.alphaMissense, v, anns)
		}

		if cr.alphaMissense != nil {
			fmt.Fprintln(w, "Gene\tTranscript\tCanon\tConsequence\tImpact\tHGVSc\tHGVSp\tAM_Score\tAM_Class")
			for _, a := range anns {
				canon := "no"
				if a.IsCanonical {
					canon = "YES"
				}
				amScore := ""
				if a.AlphaMissenseScore > 0 {
					amScore = fmt.Sprintf("%.4f", a.AlphaMissenseScore)
				}
				fmt.Fprintf(w, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
					a.GeneName, a.TranscriptID, canon,
					a.Consequence, a.Impact, a.HGVSc, a.HGVSp,
					amScore, a.AlphaMissenseClass)
			}
		} else {
			fmt.Fprintln(w, "Gene\tTranscript\tCanon\tConsequence\tImpact\tHGVSc\tHGVSp")
			for _, a := range anns {
				canon := "no"
				if a.IsCanonical {
					canon = "YES"
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
						canon = "YES"
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
	if cr.alphaMissense != nil {
		defer cr.alphaMissense.Close()
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
func runMAFOutput(logger *zap.Logger, parser *maf.Parser, ann *annotate.Annotator, out *os.File, cgl oncokb.CancerGeneList, newResults *[]duckdb.VariantResult, amStore *alphamissense.Store) error {
	mafWriter := output.NewMAFWriter(out, parser.Header(), parser.Columns())

	if cgl != nil {
		mafWriter.AddExtraColumn("Gene_Type")
	}
	if amStore != nil {
		mafWriter.AddExtraColumn("AlphaMissense_Score")
		mafWriter.AddExtraColumn("AlphaMissense_Class")
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

		// Collect all results for DuckDB persistence
		if newResults != nil && len(r.Anns) > 0 {
			chrom := r.Variant.NormalizeChrom()
			for _, a := range r.Anns {
				*newResults = append(*newResults, duckdb.VariantResult{
					Chrom: chrom, Pos: r.Variant.Pos, Ref: r.Variant.Ref, Alt: r.Variant.Alt, Ann: a,
				})
			}
		}

		best := output.SelectBestAnnotation(mafAnn, r.Anns)
		if best != nil && cgl != nil {
			if ga, ok := cgl[best.GeneName]; ok {
				best.GeneType = ga.GeneType
			}
		}
		if best != nil && amStore != nil {
			enrichAlphaMissense(amStore, r.Variant, []*annotate.Annotation{best})
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
