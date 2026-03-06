// Package main provides the vibe-vep command-line tool.
package main

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"

	"go.uber.org/zap"
	"go.uber.org/zap/zapcore"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/datasource/hotspots"
	"github.com/inodb/vibe-vep/internal/datasource/oncokb"
	"github.com/inodb/vibe-vep/internal/duckdb"
	"github.com/inodb/vibe-vep/internal/genomicindex"
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
	rootCmd.AddCommand(newCompareCmd())
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

// cacheResult holds the loaded transcript cache and optional DuckDB variant store.
type cacheResult struct {
	cache   *cache.Cache
	store   *duckdb.Store // variant cache (DuckDB), nil if --no-cache
	sources []annotate.AnnotationSource
}

// closeSources closes any sources that implement io.Closer (e.g. GenomicSource).
func (cr *cacheResult) closeSources() {
	for _, src := range cr.sources {
		if gs, ok := src.(*genomicindex.GenomicSource); ok {
			gs.Store().Close()
		}
	}
}

// normalizeAssembly validates and normalizes the assembly name.
// Accepts GRCh37, GRCh38 (canonical) and common aliases hg19, hg38.
func normalizeAssembly(assembly string) (string, error) {
	switch strings.ToLower(assembly) {
	case "grch38", "hg38":
		return "GRCh38", nil
	case "grch37", "hg19":
		return "GRCh37", nil
	default:
		return "", fmt.Errorf("unsupported assembly %q (use GRCh37 or GRCh38)", assembly)
	}
}

// loadCache loads transcripts using gob transcript cache, and opens DuckDB for variant cache.
func loadCache(logger *zap.Logger, assembly string, noCache, clearCache bool) (*cacheResult, error) {
	var err error
	assembly, err = normalizeAssembly(assembly)
	if err != nil {
		return nil, err
	}
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
			logger.Warn("transcript cache load failed, falling back to GTF/FASTA (try --clear-cache to rebuild)",
				zap.Error(err))
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
		logger.Warn("could not open variant cache (try --clear-cache or delete "+dbPath+")",
			zap.Error(err))
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

	if len(cr.sources) > 0 {
		names := make([]string, len(cr.sources))
		for i, s := range cr.sources {
			names[i] = s.Name()
		}
		logger.Info("annotation sources loaded", zap.Strings("sources", names))
	}

	return cr, nil
}

// buildSources creates annotation sources from config.
func buildSources(logger *zap.Logger, cacheDir, assembly string) []annotate.AnnotationSource {
	var sources []annotate.AnnotationSource

	// OncoKB cancer gene list
	if cglPath := viper.GetString("oncokb.cancer-gene-list"); cglPath != "" {
		cgl, err := oncokb.LoadCancerGeneList(cglPath)
		if err != nil {
			logger.Warn("could not load cancer gene list (check oncokb.cancer-gene-list path in config)",
				zap.String("path", cglPath), zap.Error(err))
		} else {
			logger.Info("loaded cancer gene list", zap.Int("genes", len(cgl)))
			sources = append(sources, oncokb.NewSource(cgl))
		}
	}

	// Unified genomic index (AlphaMissense + ClinVar + SIGNAL)
	needGenomic := viper.GetBool("annotations.alphamissense") || viper.GetBool("annotations.clinvar") ||
		(viper.GetBool("annotations.signal") && assembly == "grch37")
	if needGenomic {
		gs, err := loadGenomicIndex(logger, cacheDir, assembly)
		if err != nil {
			logger.Warn("could not load genomic index (try: vibe-vep prepare --assembly "+assembly+")",
				zap.Error(err))
		} else {
			sources = append(sources, gs)
		}
	}

	// Cancer Hotspots
	if hotspotsPath := viper.GetString("annotations.hotspots"); hotspotsPath != "" {
		store, err := hotspots.Load(hotspotsPath)
		if err != nil {
			logger.Warn("could not load hotspots data (check annotations.hotspots path in config)",
				zap.String("path", hotspotsPath), zap.Error(err))
		} else {
			logger.Info("loaded cancer hotspots", zap.Int("transcripts", store.TranscriptCount()), zap.Int("hotspots", store.HotspotCount()))
			sources = append(sources, hotspots.NewSource(store))
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

// genomicIndexSources returns the BuildSources config for the given assembly and cache dir.
func genomicIndexSources(cacheDir, assembly string) genomicindex.BuildSources {
	bs := genomicindex.BuildSources{
		AlphaMissenseTSV: filepath.Join(cacheDir, AlphaMissenseFileName(assembly)),
		ClinVarVCF:       filepath.Join(cacheDir, ClinVarFileName),
		SignalTSV:        filepath.Join(cacheDir, SignalFileName),
	}
	return bs
}

// genomicIndexPath returns the path to the unified SQLite genomic index.
func genomicIndexPath(cacheDir string) string {
	return filepath.Join(cacheDir, "genomic_annotations.sqlite")
}

// loadGenomicIndex opens (or builds) the unified genomic annotation index.
func loadGenomicIndex(logger *zap.Logger, cacheDir, assembly string) (*genomicindex.GenomicSource, error) {
	dbPath := genomicIndexPath(cacheDir)
	bs := genomicIndexSources(cacheDir, assembly)

	if !genomicindex.Ready(dbPath, bs) {
		logger.Info("building genomic index (this may take several minutes)...")
		start := time.Now()
		if err := genomicindex.Build(dbPath, bs, func(msg string, args ...any) {
			logger.Info(fmt.Sprintf(msg, args...))
		}); err != nil {
			return nil, fmt.Errorf("build genomic index: %w\nHint: ensure source data files exist in %s (run: vibe-vep download --assembly %s)", err, cacheDir, assembly)
		}
		logger.Info("built genomic index", zap.Duration("elapsed", time.Since(start)))
	} else {
		logger.Info("genomic index up to date")
	}

	store, err := genomicindex.Open(dbPath)
	if err != nil {
		return nil, fmt.Errorf("open genomic index: %w", err)
	}

	return genomicindex.NewSource(store, "1.0"), nil
}
