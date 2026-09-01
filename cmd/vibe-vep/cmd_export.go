package main

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"time"

	"go.uber.org/zap"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/duckdb"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/output"
	pqexport "github.com/inodb/vibe-vep/internal/parquet"
	"github.com/inodb/vibe-vep/internal/vcf"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

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
	addAnnotateFlags(cmd)

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
	ann.SetDistance(configuredDistance())
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
		return fmt.Errorf("opening variant cache: %w\nHint: run an annotation with --save-results first to populate the cache", err)
	}
	defer store.Close()

	logger.Info("exporting from cache", zap.String("db", dbPath))

	results, err := store.ExportAllRows()
	if err != nil {
		return fmt.Errorf("reading cached results: %w", err)
	}
	if len(results) == 0 {
		return fmt.Errorf("no variants in cache\nHint: run an annotation with --save-results first, then export")
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
