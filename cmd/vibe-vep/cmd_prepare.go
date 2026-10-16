package main

import (
	"fmt"
	"os"
	"path/filepath"
	"time"

	"go.uber.org/zap"

	"github.com/inodb/vibe-vep/internal/datasource/mutationassessor"
	"github.com/inodb/vibe-vep/internal/genomicindex"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

func newPrepareCmd(verbose *bool) *cobra.Command {
	var assembly string

	cmd := &cobra.Command{
		Use:   "prepare",
		Short: "Build transcript cache and genomic index for fast startup",
		Long: `Load GENCODE GTF/FASTA and build the transcript cache and genomic annotation
index so subsequent annotate/compare runs start instantly.`,
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

			// Build the genomic annotation index.
			asm, err := normalizeAssembly(viper.GetString("assembly"))
			if err != nil {
				return err
			}
			cacheDir := DefaultGENCODEPath(asm)
			dbPath := genomicIndexPath(cacheDir)
			bs := genomicIndexSources(cacheDir, asm)

			if !genomicindex.Ready(dbPath, bs) {
				logger.Info("building genomic index...")
				start := time.Now()
				if err := genomicindex.Build(dbPath, bs, func(msg string, args ...any) {
					logger.Info(fmt.Sprintf(msg, args...))
				}); err != nil {
					return fmt.Errorf("build genomic index: %w", err)
				}
				logger.Info("built genomic index", zap.Duration("elapsed", time.Since(start)))
			} else {
				logger.Info("genomic index already up to date")
			}

			// Build the Mutation Assessor index, if its files were fetched.
			if err := prepareMutationAssessor(logger, cacheDir); err != nil {
				return err
			}

			return nil
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")

	return cmd
}

// prepareMutationAssessor builds the Mutation Assessor index from the
// downloaded files. Missing files are not an error: the source is optional,
// and saying so once is more use than failing the whole prepare.
func prepareMutationAssessor(logger *zap.Logger, cacheDir string) error {
	dbPath := filepath.Join(cacheDir, MutationAssessorDBName)
	files := MutationAssessorFiles(rawDirForCache(cacheDir))

	for _, f := range files {
		if _, err := os.Stat(f); err != nil {
			logger.Info("Mutation Assessor files not downloaded, skipping",
				zap.String("missing", f))
			return nil
		}
	}
	if mutationassessor.Ready(dbPath) {
		logger.Info("Mutation Assessor index already up to date")
		return nil
	}

	logger.Info("building Mutation Assessor index...")
	start := time.Now()
	if err := mutationassessor.Build(dbPath, files); err != nil {
		return fmt.Errorf("build Mutation Assessor index: %w", err)
	}
	logger.Info("built Mutation Assessor index", zap.Duration("elapsed", time.Since(start)))
	return nil
}
