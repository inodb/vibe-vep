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
	"github.com/inodb/vibe-vep/internal/datasource/clinvar"
	"github.com/inodb/vibe-vep/internal/datasource/hotspots"
	"github.com/inodb/vibe-vep/internal/datasource/oncokb"
	"github.com/inodb/vibe-vep/internal/datasource/signal"
	"github.com/inodb/vibe-vep/internal/duckdb"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/output"
	pqexport "github.com/inodb/vibe-vep/internal/parquet"
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

// isColorTerminal returns true if stdout appears to be a color-capable terminal.
func isColorTerminal() bool {
	if os.Getenv("NO_COLOR") != "" {
		return false
	}
	fi, err := os.Stdout.Stat()
	if err != nil {
		return false
	}
	return fi.Mode()&os.ModeCharDevice != 0
}

func banner() string {
	const plain = "\n" +
		"  ██╗   ██╗██╗██████╗ ███████╗    ██╗   ██╗███████╗██████╗ \n" +
		"  ██║   ██║██║██╔══██╗██╔════╝    ██║   ██║██╔════╝██╔══██╗\n" +
		"  ██║   ██║██║██████╔╝█████╗      ██║   ██║█████╗  ██████╔╝\n" +
		"  ╚██╗ ██╔╝██║██╔══██╗██╔══╝     ╚██╗ ██╔╝██╔══╝  ██╔═══╝ \n" +
		"   ╚████╔╝ ██║██████╔╝███████╗    ╚████╔╝ ███████╗██║     \n" +
		"    ╚═══╝  ╚═╝╚═════╝ ╚══════╝    ╚═══╝  ╚══════╝╚═╝     \n" +
		"\n" +
		"  Variant Effect Predictor and Annotator for Oncology.\n" +
		"  Combines Ensembl VEP, Genome Nexus, and other annotation tools into one binary."

	if !isColorTerminal() {
		return plain
	}

	return "\n" +
		"\033[96m  ██╗   ██╗██╗██████╗ ███████╗    ██╗   ██╗███████╗██████╗ \n" +
		"\033[36m  ██║   ██║██║██╔══██╗██╔════╝    ██║   ██║██╔════╝██╔══██╗\n" +
		"\033[94m  ██║   ██║██║██████╔╝█████╗      ██║   ██║█████╗  ██████╔╝\n" +
		"\033[34m  ╚██╗ ██╔╝██║██╔══██╗██╔══╝     ╚██╗ ██╔╝██╔══╝  ██╔═══╝ \n" +
		"\033[95m   ╚████╔╝ ██║██████╔╝███████╗    ╚████╔╝ ███████╗██║     \n" +
		"\033[35m    ╚═══╝  ╚═╝╚═════╝ ╚══════╝    ╚═══╝  ╚══════╝╚═╝     \n" +
		"\033[0m\n" +
		"\033[2m  Variant Effect Predictor and Annotator for Oncology.\n" +
		"  Combines Ensembl VEP, Genome Nexus, and other annotation tools into one binary.\033[0m"
}

func newRootCmd() *cobra.Command {
	var (
		verbose    bool
		configFile string
	)

	rootCmd := &cobra.Command{
		Use:     "vibe-vep",
		Short:   "Variant Effect Predictor",
		Long:    banner(),
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
	rootCmd.AddCommand(newConvertCmd(&verbose))
	rootCmd.AddCommand(newDownloadCmd(&verbose))
	rootCmd.AddCommand(newExportCmd(&verbose))
	rootCmd.AddCommand(newPrepareCmd(&verbose))
	rootCmd.AddCommand(newVersionCmd(&verbose))

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
		pick          bool
		mostSevere    bool
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

func newExportCmd(verbose *bool) *cobra.Command {
	cmd := &cobra.Command{
		Use:   "export",
		Short: "Export annotations to various formats",
		Long:  "Export variant annotations to Parquet or other formats for downstream analysis.",
	}

	cmd.AddCommand(newExportParquetCmd(verbose))

	return cmd
}

func newExportParquetCmd(verbose *bool) *cobra.Command {
	var (
		assembly      string
		outputFile    string
		canonicalOnly bool
		pick          bool
		rowGroupSize  int
		fromCache     bool
	)

	cmd := &cobra.Command{
		Use:   "parquet [input.vcf|input.maf]",
		Short: "Export annotations as a sorted Parquet file",
		Long: `Export variant annotations to a sorted Parquet file for browser-based querying
via DuckDB-WASM. The output is sorted by (chrom, pos) for efficient row group pruning.

Use --from-cache to export directly from the DuckDB variant cache without loading
transcripts or annotation sources into memory.`,
		Example: `  vibe-vep export parquet --canonical --pick -o annotations.parquet input.maf
  vibe-vep export parquet -o output.parquet input.vcf
  vibe-vep export parquet --from-cache -o cached.parquet`,
		Args: cobra.RangeArgs(0, 1),
		PreRunE: func(cmd *cobra.Command, args []string) error {
			return viper.BindPFlags(cmd.Flags())
		},
		RunE: func(cmd *cobra.Command, args []string) error {
			logger, err := newLogger(*verbose)
			if err != nil {
				return fmt.Errorf("creating logger: %w", err)
			}
			defer logger.Sync()

			if viper.GetBool("from-cache") {
				return runExportParquetFromCache(logger,
					viper.GetString("assembly"),
					viper.GetString("output"),
					viper.GetInt("row-group-size"),
				)
			}

			if len(args) == 0 {
				return fmt.Errorf("input file required (or use --from-cache)")
			}
			return runExportParquet(logger, args[0],
				viper.GetString("assembly"),
				viper.GetString("output"),
				viper.GetBool("canonical"),
				viper.GetBool("pick"),
				viper.GetInt("row-group-size"),
				viper.GetBool("no-cache"),
				viper.GetBool("clear-cache"),
			)
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")
	cmd.Flags().StringVarP(&outputFile, "output", "o", "annotations.parquet", "Output Parquet file path")
	cmd.Flags().BoolVar(&canonicalOnly, "canonical", false, "Only report canonical transcript annotations")
	cmd.Flags().BoolVar(&pick, "pick", false, "One annotation per variant (best transcript)")
	cmd.Flags().IntVar(&rowGroupSize, "row-group-size", pqexport.DefaultRowGroupSize, "Rows per row group")
	cmd.Flags().BoolVar(&fromCache, "from-cache", false, "Export directly from DuckDB cache (no transcript loading)")
	addCacheFlags(cmd)

	return cmd
}

func runExportParquet(logger *zap.Logger, inputPath, assembly, outputFile string, canonicalOnly, pick bool, rowGroupSize int, noCache, clearCache bool) error {
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

	// Auto-detect input format by extension
	ext := strings.ToLower(filepath.Ext(inputPath))
	isVCF := ext == ".vcf" || ext == ".gz"

	var rows []pqexport.Row

	if isVCF {
		rows, err = exportVCFToRows(logger, inputPath, ann, cr.sources, pick)
	} else {
		rows, err = exportMAFToRows(logger, inputPath, ann, cr.sources, pick)
	}
	if err != nil {
		return err
	}

	if len(rows) == 0 {
		return fmt.Errorf("no variants to export")
	}

	// Sort by (chrom_numeric, pos, ref, alt, transcript_id)
	logger.Info("sorting rows", zap.Int("count", len(rows)))
	pqexport.SortRows(rows)

	// Write Parquet file
	f, err := os.Create(outputFile)
	if err != nil {
		return fmt.Errorf("creating output file: %w", err)
	}
	defer f.Close()

	logger.Info("writing Parquet file", zap.String("path", outputFile), zap.Int("rows", len(rows)), zap.Int("row_group_size", rowGroupSize))
	w := pqexport.NewWriter(f, rowGroupSize)
	if err := w.WriteRows(rows); err != nil {
		return fmt.Errorf("writing rows: %w", err)
	}
	if err := w.Close(); err != nil {
		return fmt.Errorf("closing parquet writer: %w", err)
	}

	logger.Info("export complete", zap.String("output", outputFile), zap.Int("total_rows", len(rows)))
	return nil
}

func runExportParquetFromCache(logger *zap.Logger, assembly, outputFile string, rowGroupSize int) error {
	cacheDir := DefaultGENCODEPath(assembly)
	dbPath := filepath.Join(cacheDir, "variant_cache.duckdb")

	store, err := duckdb.Open(dbPath)
	if err != nil {
		return fmt.Errorf("opening variant cache: %w", err)
	}
	defer store.Close()

	logger.Info("exporting from cache", zap.String("db", dbPath))

	results, err := store.ExportAllRows()
	if err != nil {
		return fmt.Errorf("reading cached results: %w", err)
	}
	if len(results) == 0 {
		return fmt.Errorf("no variants in cache")
	}

	rows := make([]pqexport.Row, len(results))
	for i, r := range results {
		rows[i] = pqexport.AnnotationToRow(r.Chrom, r.Pos, r.Ref, r.Alt, r.Ann)
	}

	// Already sorted by DuckDB ORDER BY, but apply numeric chrom sort
	logger.Info("sorting rows", zap.Int("count", len(rows)))
	pqexport.SortRows(rows)

	f, err := os.Create(outputFile)
	if err != nil {
		return fmt.Errorf("creating output file: %w", err)
	}
	defer f.Close()

	logger.Info("writing Parquet file", zap.String("path", outputFile), zap.Int("rows", len(rows)), zap.Int("row_group_size", rowGroupSize))
	w := pqexport.NewWriter(f, rowGroupSize)
	if err := w.WriteRows(rows); err != nil {
		return fmt.Errorf("writing rows: %w", err)
	}
	if err := w.Close(); err != nil {
		return fmt.Errorf("closing parquet writer: %w", err)
	}

	logger.Info("export complete", zap.String("output", outputFile), zap.Int("total_rows", len(rows)))
	return nil
}

func exportVCFToRows(logger *zap.Logger, inputPath string, ann *annotate.Annotator, sources []annotate.AnnotationSource, pick bool) ([]pqexport.Row, error) {
	parser, err := vcf.NewParser(inputPath)
	if err != nil {
		return nil, err
	}
	defer parser.Close()

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
			items <- annotate.WorkItem{Seq: seq, Variant: v}
			seq++
		}
	}()

	results := ann.ParallelAnnotate(items, 0)

	var rows []pqexport.Row
	progress := func(n int) {
		logger.Info("progress", zap.Int("variants_processed", n))
	}

	if err := annotate.OrderedCollectWithProgress(results, 2*time.Second, progress, func(r annotate.WorkResult) error {
		if r.Err != nil {
			logger.Warn("annotation failed", zap.Error(r.Err))
			return nil
		}

		anns := r.Anns
		for _, src := range sources {
			src.Annotate(r.Variant, anns)
		}

		if pick && len(anns) > 1 {
			anns = []*annotate.Annotation{output.PickBestAnnotation(anns)}
		}

		chrom := r.Variant.NormalizeChrom()
		for _, a := range anns {
			rows = append(rows, pqexport.AnnotationToRow(chrom, r.Variant.Pos, r.Variant.Ref, r.Variant.Alt, a))
		}
		return nil
	}); err != nil {
		return nil, err
	}

	if parseErr != nil {
		return nil, parseErr
	}

	return rows, nil
}

func exportMAFToRows(logger *zap.Logger, inputPath string, ann *annotate.Annotator, sources []annotate.AnnotationSource, pick bool) ([]pqexport.Row, error) {
	parser, err := maf.NewParser(inputPath)
	if err != nil {
		return nil, err
	}
	defer parser.Close()

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

	var rows []pqexport.Row
	progress := func(n int) {
		logger.Info("progress", zap.Int("variants_processed", n))
	}

	if err := annotate.OrderedCollectWithProgress(results, 2*time.Second, progress, func(r annotate.WorkResult) error {
		if r.Err != nil {
			logger.Warn("annotation failed", zap.Error(r.Err))
			return nil
		}

		anns := r.Anns

		// For MAF, select best annotation per variant (matching MAF behavior)
		if pick && len(anns) > 0 {
			best := output.PickBestAnnotation(anns)
			if best != nil {
				anns = []*annotate.Annotation{best}
			}
		}

		for _, src := range sources {
			src.Annotate(r.Variant, anns)
		}

		chrom := r.Variant.NormalizeChrom()
		for _, a := range anns {
			rows = append(rows, pqexport.AnnotationToRow(chrom, r.Variant.Pos, r.Variant.Ref, r.Variant.Alt, a))
		}
		return nil
	}); err != nil {
		return nil, err
	}

	if parseErr != nil {
		return nil, parseErr
	}

	return rows, nil
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
	cache   *cache.Cache
	store   *duckdb.Store // variant cache (DuckDB), nil if --no-cache
	sources []annotate.AnnotationSource
}

// closeSources closes any sources that implement io.Closer (e.g. AlphaMissense Store).
func (cr *cacheResult) closeSources() {
	for _, src := range cr.sources {
		if amSrc, ok := src.(*alphamissense.Source); ok {
			amSrc.Store().Close()
		}
	}
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

	// Build interval tree index for O(log n) transcript lookup
	{
		start := time.Now()
		c.BuildIndex()
		logger.Debug("built interval tree index", zap.Duration("elapsed", time.Since(start)))
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

	// --- Build annotation sources ---
	cr.sources = buildSources(logger, cacheDir, assembly)

	return cr, nil
}

// buildSources creates annotation sources from config.
func buildSources(logger *zap.Logger, cacheDir, assembly string) []annotate.AnnotationSource {
	var sources []annotate.AnnotationSource

	// OncoKB cancer gene list
	if cglPath := viper.GetString("oncokb.cancer-gene-list"); cglPath != "" {
		cgl, err := oncokb.LoadCancerGeneList(cglPath)
		if err != nil {
			logger.Warn("could not load cancer gene list", zap.String("path", cglPath), zap.Error(err))
		} else {
			logger.Info("loaded cancer gene list", zap.Int("genes", len(cgl)))
			sources = append(sources, oncokb.NewSource(cgl))
		}
	}

	// AlphaMissense
	if viper.GetBool("annotations.alphamissense") {
		amStore, err := loadAlphaMissense(logger, cacheDir, assembly)
		if err != nil {
			logger.Warn("could not load AlphaMissense data", zap.Error(err))
		} else {
			sources = append(sources, alphamissense.NewSource(amStore))
		}
	}

	// ClinVar
	if viper.GetBool("annotations.clinvar") {
		cvStore, err := loadClinVar(logger, cacheDir)
		if err != nil {
			logger.Warn("could not load ClinVar data", zap.Error(err))
		} else {
			sources = append(sources, clinvar.NewSource(cvStore))
		}
	}

	// Cancer Hotspots
	if hotspotsPath := viper.GetString("annotations.hotspots"); hotspotsPath != "" {
		store, err := hotspots.Load(hotspotsPath)
		if err != nil {
			logger.Warn("could not load hotspots data", zap.String("path", hotspotsPath), zap.Error(err))
		} else {
			logger.Info("loaded cancer hotspots", zap.Int("transcripts", store.TranscriptCount()), zap.Int("hotspots", store.HotspotCount()))
			sources = append(sources, hotspots.NewSource(store))
		}
	}

	// SIGNAL (GRCh37 only — no liftover support yet)
	if viper.GetBool("annotations.signal") {
		if assembly != "grch37" {
			logger.Warn("SIGNAL data uses GRCh37 coordinates, skipping for assembly", zap.String("assembly", assembly))
		} else {
			sigStore, err := loadSignal(logger, cacheDir)
			if err != nil {
				logger.Warn("could not load SIGNAL data", zap.Error(err))
			} else {
				sources = append(sources, signal.NewSource(sigStore))
			}
		}
	}

	return sources
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

	// Preload into memory for fast lookups (eliminates per-variant DuckDB overhead)
	logger.Info("preloading AlphaMissense data into memory...")
	start := time.Now()
	if err := amStore.PreloadToMemory(); err != nil {
		logger.Warn("could not preload AlphaMissense to memory, using DuckDB fallback", zap.Error(err))
	} else {
		logger.Info("preloaded AlphaMissense data",
			zap.Int64("variants", amStore.MemCacheSize()),
			zap.Duration("elapsed", time.Since(start)))
	}

	return amStore, nil
}

// loadClinVar loads ClinVar data from the VCF file.
// ClinVar has 4M+ entries — too large for gob caching on typical machines.
// Parsing the gzipped VCF takes ~25s, which is acceptable for an opt-in source.
func loadClinVar(logger *zap.Logger, cacheDir string) (*clinvar.Store, error) {
	vcfPath := filepath.Join(cacheDir, ClinVarFileName)
	if _, err := os.Stat(vcfPath); err != nil {
		return nil, fmt.Errorf("ClinVar VCF not found at %s\nHint: run 'vibe-vep download' to fetch it", vcfPath)
	}

	logger.Info("parsing ClinVar VCF...", zap.String("path", vcfPath))
	start := time.Now()
	store, err := clinvar.Load(vcfPath)
	if err != nil {
		return nil, fmt.Errorf("load ClinVar data: %w", err)
	}
	logger.Info("loaded ClinVar data",
		zap.Int("variants", store.Count()),
		zap.Duration("elapsed", time.Since(start)))

	return store, nil
}

// loadSignal loads SIGNAL data, using gob cache if available.
func loadSignal(logger *zap.Logger, cacheDir string) (*signal.Store, error) {
	tsvPath := filepath.Join(cacheDir, SignalFileName)
	gobPath := filepath.Join(cacheDir, "signal.gob")

	// Try gob cache first
	if signal.GobValid(gobPath, tsvPath) {
		logger.Info("loading SIGNAL from cache...")
		start := time.Now()
		store, err := signal.LoadGob(gobPath)
		if err == nil {
			logger.Info("loaded SIGNAL from cache",
				zap.Int("variants", store.Count()),
				zap.Duration("elapsed", time.Since(start)))
			return store, nil
		}
		logger.Warn("could not load SIGNAL cache, re-parsing TSV", zap.Error(err))
	}

	// Parse TSV
	if _, err := os.Stat(tsvPath); err != nil {
		return nil, fmt.Errorf("SIGNAL TSV not found at %s\nHint: run 'vibe-vep download' to fetch it", tsvPath)
	}

	logger.Info("parsing SIGNAL data...", zap.String("path", tsvPath))
	start := time.Now()
	store, err := signal.Load(tsvPath)
	if err != nil {
		return nil, fmt.Errorf("load SIGNAL data: %w", err)
	}
	logger.Info("loaded SIGNAL data",
		zap.Int("variants", store.Count()),
		zap.Duration("elapsed", time.Since(start)))

	// Save gob cache for fast reload
	if err := store.SaveGob(gobPath); err != nil {
		logger.Warn("could not save SIGNAL cache", zap.Error(err))
	}

	return store, nil
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
			cr.closeSources()
			logger.Info("transcript cache ready",
				zap.Int("transcripts", cr.cache.TranscriptCount()))
			return nil
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")

	return cmd
}

func newVersionCmd(verbose *bool) *cobra.Command {
	var mafColumns bool

	cmd := &cobra.Command{
		Use:   "version",
		Short: "Show version and data source information",
		Long:  "Show vibe-vep version, loaded data sources, and optionally the MAF column mapping.",
		PreRunE: func(cmd *cobra.Command, args []string) error {
			return viper.BindPFlags(cmd.Flags())
		},
		RunE: func(cmd *cobra.Command, args []string) error {
			fmt.Printf("%-15s%s\n", "vibe-vep", version)

			// Show GENCODE info
			assembly := viper.GetString("assembly")
			if assembly == "" {
				assembly = "GRCh38"
			}
			if _, _, _, found := FindGENCODEFiles(assembly); found {
				fmt.Printf("%-15s%s (%s)\n", "GENCODE", GencodeVersion, assembly)
			}

			// Show configured annotation sources (lightweight, no data loading).
			cacheDir := DefaultGENCODEPath(assembly)
			type sourceInfo struct {
				name, match, assembly, version, status string
				columns                                []annotate.ColumnDef
			}
			var infos []sourceInfo

			// fileModDate returns the file modification date as YYYY-MM-DD, or "" on error.
			fileModDate := func(path string) string {
				fi, err := os.Stat(path)
				if err != nil {
					return ""
				}
				return fi.ModTime().Format("2006-01-02")
			}

			// OncoKB
			if cglPath := viper.GetString("oncokb.cancer-gene-list"); cglPath != "" {
				status := "configured"
				ver := ""
				if fi, err := os.Stat(cglPath); err != nil {
					status = "file not found"
				} else {
					ver = fi.ModTime().Format("2006-01-02")
				}
				infos = append(infos, sourceInfo{"oncokb", string(annotate.MatchGene), "any", ver, status,
					[]annotate.ColumnDef{{Name: "gene_type", Description: "Gene classification (ONCOGENE/TSG)"}}})
			}

			// AlphaMissense
			if viper.GetBool("annotations.alphamissense") {
				dbPath := filepath.Join(cacheDir, "annotations.duckdb")
				status := "configured"
				ver := "2023"
				if _, err := os.Stat(dbPath); err != nil {
					status = "not prepared (run: vibe-vep prepare)"
				}
				infos = append(infos, sourceInfo{"alphamissense", string(annotate.MatchGenomic), "GRCh38", ver, status,
					[]annotate.ColumnDef{
						{Name: "score", Description: "Pathogenicity score (0-1)"},
						{Name: "class", Description: "likely_benign/ambiguous/likely_pathogenic"},
					}})
			}

			// ClinVar
			if viper.GetBool("annotations.clinvar") {
				cvPath := filepath.Join(cacheDir, ClinVarFileName)
				status := "configured"
				ver := fileModDate(cvPath)
				if ver == "" {
					status = "not downloaded (run: vibe-vep download)"
				}
				infos = append(infos, sourceInfo{"clinvar", string(annotate.MatchGenomic), "GRCh38", ver, status,
					[]annotate.ColumnDef{
						{Name: "clnsig", Description: "Clinical significance (e.g. Pathogenic, Benign)"},
						{Name: "clnrevstat", Description: "Review status"},
						{Name: "clndn", Description: "Disease name(s)"},
					}})
			}

			// Hotspots
			if hsPath := viper.GetString("annotations.hotspots"); hsPath != "" {
				status := "configured"
				ver := fileModDate(hsPath)
				if ver == "" {
					status = "file not found"
				}
				infos = append(infos, sourceInfo{"hotspots", string(annotate.MatchProteinPosition), "any", ver, status,
					[]annotate.ColumnDef{
						{Name: "hotspot", Description: "Y if position is a known cancer hotspot"},
						{Name: "type", Description: "Hotspot type: single residue, in-frame indel, 3d, splice"},
						{Name: "qvalue", Description: "Statistical significance (q-value)"},
					}})
			}

			// SIGNAL
			if viper.GetBool("annotations.signal") {
				if assembly != "grch37" {
					infos = append(infos, sourceInfo{"signal", string(annotate.MatchGenomic), "GRCh37", "", "skipped (GRCh37 only)", nil})
				} else {
					status := "configured"
					sigPath := filepath.Join(cacheDir, SignalFileName)
					ver := fileModDate(sigPath)
					if ver == "" {
						status = "not downloaded"
					}
					infos = append(infos, sourceInfo{"signal", string(annotate.MatchGenomic), "GRCh37", ver, status,
						[]annotate.ColumnDef{
							{Name: "mutation_status", Description: "Germline mutation status"},
							{Name: "count_carriers", Description: "Number of carriers in SIGNAL cohort"},
							{Name: "frequency", Description: "Overall allele frequency in SIGNAL cohort"},
						}})
				}
			}

			if len(infos) > 0 {
				fmt.Println()
				fmt.Println("Annotation Sources:")
				tw := tabwriter.NewWriter(os.Stdout, 0, 4, 2, ' ', 0)
				fmt.Fprintln(tw, "  NAME\tMATCH\tASSEMBLY\tVERSION\tSTATUS")
				for _, info := range infos {
					fmt.Fprintf(tw, "  %s\t%s\t%s\t%s\t%s\n", info.name, info.match, info.assembly, info.version, info.status)
				}
				tw.Flush()
			}

			if mafColumns {
				fmt.Println()
				fmt.Println("MAF Output Columns:")
				fmt.Println()
				tw := tabwriter.NewWriter(os.Stdout, 0, 4, 2, ' ', 0)
				fmt.Fprintln(tw, "COLUMN\tSOURCE\tDESCRIPTION")
				for _, col := range annotate.CoreColumns {
					fmt.Fprintf(tw, "vibe.%s\tvibe-vep\t%s\n", col.Name, col.Description)
				}
				for _, info := range infos {
					for _, col := range info.columns {
						fmt.Fprintf(tw, "vibe.%s.%s\t%s\t%s\n", info.name, col.Name, info.name, col.Description)
					}
				}
				tw.Flush()
			}

			return nil
		},
	}

	cmd.Flags().BoolVar(&mafColumns, "maf-columns", false, "Show MAF output column mapping")

	return cmd
}

func runAnnotateMAF(logger *zap.Logger, inputPath, assembly, outputFile string, canonicalOnly, saveResults, noCache, clearCache, mostSevere bool) error {
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

	if err := runMAFOutput(logger, parser, ann, out, cr.sources, collectResults, mostSevere); err != nil {
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
				for _, col := range src.Columns() {
					header += "\t" + src.Name() + "_" + col.Name
				}
			}
			fmt.Fprintln(w, header)
			for _, a := range anns {
				canon := "no"
				if a.IsCanonical {
					canon = "YES"
				}
				line := fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s",
					a.GeneName, a.TranscriptID, canon,
					a.Consequence, a.Impact, a.HGVSc, a.HGVSp)
				for _, src := range cr.sources {
					for _, col := range src.Columns() {
						line += "\t" + a.GetExtra(src.Name(), col.Name)
					}
				}
				fmt.Fprintln(w, line)
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
	defer cr.closeSources()

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
			return writer.WriteRow(r.Variant, nil)
		}

		for _, src := range cr.sources {
			src.Annotate(r.Variant, r.Anns)
		}

		best := output.PickBestAnnotation(r.Anns)
		return writer.WriteRow(r.Variant, best)
	}); err != nil {
		return err
	}

	if parseErr != nil {
		return parseErr
	}

	return writer.Flush()
}

// runMAFOutput runs MAF annotation mode, preserving all original columns.
func runMAFOutput(logger *zap.Logger, parser *maf.Parser, ann *annotate.Annotator, out *os.File, sources []annotate.AnnotationSource, newResults *[]duckdb.VariantResult, mostSevere bool) error {
	mafWriter := output.NewMAFWriter(out, parser.Header(), parser.Columns())
	mafWriter.SetSources(sources)

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
