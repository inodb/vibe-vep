package main

import (
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"time"

	"go.uber.org/zap"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/datasource/ensembltranscript"
	"github.com/inodb/vibe-vep/internal/datasource/genomeref"
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

			// Build the Ensembl transcript index, if its file was fetched.
			if err := prepareEnsemblTranscripts(logger, cacheDir); err != nil {
				return err
			}

			// Build the Mutation Assessor index, if its files were fetched.
			if err := prepareMutationAssessor(logger, cacheDir); err != nil {
				return err
			}

			// Build the reference genome index, if the genomic FASTA is here.
			if err := prepareGenomeReference(logger, cacheDir); err != nil {
				return err
			}

			return nil
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")

	return cmd
}

// prepareEnsemblTranscripts builds the Ensembl transcript index the /ensembl
// endpoints report from. A missing source file is not an error: those
// endpoints fall back to the GENCODE cache, which answers with less.
func prepareEnsemblTranscripts(logger *zap.Logger, cacheDir string) error {
	source := filepath.Join(rawDirForCache(cacheDir), EnsemblTranscriptsFileName)
	if _, err := os.Stat(source); err != nil {
		logger.Info("Ensembl transcript records not downloaded, skipping",
			zap.String("missing", source))
		return nil
	}

	dbPath := filepath.Join(cacheDir, EnsemblTranscriptDBName)
	logger.Info("building Ensembl transcript index...")
	start := time.Now()
	if err := ensembltranscript.Build(source, dbPath, func(msg string, args ...any) {
		logger.Info(fmt.Sprintf(msg, args...))
	}); err != nil {
		return fmt.Errorf("build Ensembl transcript index: %w", err)
	}
	logger.Info("built Ensembl transcript index", zap.Duration("elapsed", time.Since(start)))
	return nil
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

// prepareGenomeReference builds the reference index hgvsg is renormalised
// against. The genomic FASTA is optional, as it is for the intron flanks, so
// a missing one is not an error: hgvsg then names a change where the query
// put it rather than where HGVS requires.
func prepareGenomeReference(logger *zap.Logger, cacheDir string) error {
	fasta := findGenomeFASTA(rawDirForCache(cacheDir))
	if fasta == "" {
		logger.Info("genomic FASTA not downloaded, skipping the reference index")
		return nil
	}
	path := filepath.Join(cacheDir, genomeref.FileName)
	if store, err := genomeref.Open(path); err == nil {
		n := store.Chromosomes()
		store.Close()
		if n > 0 {
			logger.Info("reference index already built", zap.Int("chromosomes", n))
			return nil
		}
	}

	logger.Info("building reference index...", zap.String("from", fasta))
	start := time.Now()
	if err := genomeref.Build(fasta, path, wantedChromosome, cache.StreamGenomeFASTA); err != nil {
		return fmt.Errorf("build reference index: %w", err)
	}
	logger.Info("built reference index", zap.Duration("elapsed", time.Since(start)))
	return nil
}

// wantedChromosome accepts the assembled chromosomes and turns away the
// scaffolds and patches, which carry no annotation and would double the
// index for nothing.
func wantedChromosome(chrom string) bool {
	switch chrom {
	case "X", "Y", "M", "MT":
		return true
	}
	if n, err := strconv.Atoi(chrom); err == nil {
		return n >= 1 && n <= 22
	}
	return false
}
