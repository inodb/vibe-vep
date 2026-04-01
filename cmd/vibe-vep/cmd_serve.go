package main

import (
	"context"
	"fmt"
	"net"
	"net/http"
	"os"
	"os/signal"
	"path/filepath"
	"strings"
	"syscall"
	"time"

	"go.uber.org/zap"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/datasource/pfam"
	"github.com/inodb/vibe-vep/internal/datasource/ptm"
	"github.com/inodb/vibe-vep/internal/datasource/uniprot"
	"github.com/inodb/vibe-vep/internal/server"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

func newServeCmd(verbose *bool) *cobra.Command {
	var (
		assembly     string
		port         int
		host         string
		readTimeout  time.Duration
		writeTimeout time.Duration
	)

	cmd := &cobra.Command{
		Use:   "serve",
		Short: "Start HTTP annotation server",
		Long: `Start an HTTP server that provides variant annotation endpoints.

Loads transcript cache and annotation sources at startup, then serves
annotation requests over HTTP. Supports both Ensembl VEP REST API and
genome-nexus API compatible endpoints.

Endpoints:
  Ensembl VEP compatibility:
    GET  /ensembl/{assembly}/vep/human/region/{region}/{allele}
    POST /ensembl/{assembly}/vep/human/region
    GET  /ensembl/{assembly}/vep/human/hgvs/{notation}
    POST /ensembl/{assembly}/vep/human/hgvs

  Genome-nexus compatibility:
    GET  /genome-nexus/{assembly}/annotation/genomic/{genomicLocation}
    POST /genome-nexus/{assembly}/annotation/genomic
    GET  /genome-nexus/{assembly}/annotation/{variant}
    POST /genome-nexus/{assembly}/annotation

    Supported ?fields= enrichments (comma-separated):
      annotation_summary  Canonical transcript summary with variantClassification,
                          hgvspShort, variantType, aminoAcidRef/Alt
      clinvar             ClinVar clinical significance
      hotspots            Cancer mutation hotspots

    Always included when available:
      colocatedVariants   dbSNP RS identifiers
      hgnc_id             HGNC gene identifier (per transcript)
      protein_id          Ensembl protein ID (per transcript)
      sift/polyphen       SIFT and PolyPhen-2 scores (per transcript)

    Not yet implemented (genome-nexus fields):
      mutation_assessor, my_variant_info, oncokb, ptms, nucleotide_context
      entrezGeneId, refseq_transcript_ids (per transcript)

  Health/info:
    GET  /health
    GET  /info`,
		Example: `  # Start server (single assembly)
  vibe-vep serve --assembly GRCh38 --port 8080

  # Start server (both assemblies)
  vibe-vep serve --assembly GRCh38,GRCh37 --port 8080

  # Test endpoints
  curl "http://localhost:8080/ensembl/grch38/vep/human/region/7:140753336-140753336:1/T"
  curl "http://localhost:8080/genome-nexus/grch38/annotation/genomic/7,140753336,140753336,A,T"`,
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
			return runServe(logger, runServeConfig{
				assemblies:   viper.GetString("assembly"),
				host:         viper.GetString("host"),
				port:         viper.GetInt("port"),
				readTimeout:  viper.GetDuration("read-timeout"),
				writeTimeout: viper.GetDuration("write-timeout"),
				noCache:      viper.GetBool("no-cache"),
				clearCache:   viper.GetBool("clear-cache"),
			})
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Comma-separated assemblies to load (e.g. GRCh38,GRCh37)")
	cmd.Flags().IntVar(&port, "port", 8080, "Port to listen on")
	cmd.Flags().StringVar(&host, "host", "0.0.0.0", "Host to bind to")
	cmd.Flags().DurationVar(&readTimeout, "read-timeout", 30*time.Second, "HTTP read timeout")
	cmd.Flags().DurationVar(&writeTimeout, "write-timeout", 60*time.Second, "HTTP write timeout")
	addCacheFlags(cmd)

	return cmd
}

type runServeConfig struct {
	assemblies   string
	host         string
	port         int
	readTimeout  time.Duration
	writeTimeout time.Duration
	noCache      bool
	clearCache   bool
}

func runServe(logger *zap.Logger, cfg runServeConfig) error {
	srv := server.New(logger, version)

	// Load each assembly.
	assemblyNames := strings.Split(cfg.assemblies, ",")
	var cacheResults []*cacheResult

	for _, name := range assemblyNames {
		name = strings.TrimSpace(name)
		if name == "" {
			continue
		}

		normalized, err := normalizeAssembly(name)
		if err != nil {
			return err
		}

		logger.Info("loading assembly", zap.String("assembly", normalized))
		start := time.Now()

		cr, err := loadCache(logger, normalized, cfg.noCache, cfg.clearCache)
		if err != nil {
			return fmt.Errorf("loading assembly %s: %w", normalized, err)
		}
		cacheResults = append(cacheResults, cr)

		ann := annotate.NewAnnotator(cr.cache)
		ann.SetLogger(logger)

		srv.AddAssembly(normalized, cr.cache, ann, cr.sources)

		// Load PFAM domain data if available.
		if pfamStore := loadPfamStore(logger, normalized); pfamStore != nil {
			srv.SetPfamStore(normalized, pfamStore)
		}

		// Load PTM data if available.
		if ptmStore := loadPtmStore(logger, normalized); ptmStore != nil {
			srv.SetPtmStore(normalized, ptmStore)
		}

		// Load UniProt mapping if available.
		if uniprotStore := loadUniprotStore(logger, normalized); uniprotStore != nil {
			srv.SetUniprotStore(normalized, uniprotStore)
		}

		logger.Info("assembly loaded",
			zap.String("assembly", normalized),
			zap.Int("transcripts", cr.cache.TranscriptCount()),
			zap.Int("sources", len(cr.sources)),
			zap.Duration("elapsed", time.Since(start)))
	}

	// Ensure cleanup on shutdown.
	defer func() {
		for _, cr := range cacheResults {
			if cr.store != nil {
				cr.store.Close()
			}
			cr.closeSources()
		}
	}()

	addr := net.JoinHostPort(cfg.host, fmt.Sprintf("%d", cfg.port))

	httpServer := &http.Server{
		Addr:         addr,
		Handler:      srv.Handler(),
		ReadTimeout:  cfg.readTimeout,
		WriteTimeout: cfg.writeTimeout,
	}

	// Graceful shutdown on SIGINT/SIGTERM.
	errCh := make(chan error, 1)
	go func() {
		logger.Info("server starting", zap.String("addr", addr))
		errCh <- httpServer.ListenAndServe()
	}()

	sigCh := make(chan os.Signal, 1)
	signal.Notify(sigCh, syscall.SIGINT, syscall.SIGTERM)

	select {
	case err := <-errCh:
		if err != http.ErrServerClosed {
			return fmt.Errorf("server error: %w", err)
		}
	case sig := <-sigCh:
		logger.Info("shutting down", zap.String("signal", sig.String()))
		ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
		defer cancel()
		if err := httpServer.Shutdown(ctx); err != nil {
			return fmt.Errorf("shutdown error: %w", err)
		}
	}

	return nil
}

// loadPfamStore loads PFAM domain data from the raw download directory.
func loadPfamStore(logger *zap.Logger, assembly string) *pfam.Store {
	rawDir := RawDir(assembly)
	if rawDir == "" {
		return nil
	}

	pfamAPath := filepath.Join(rawDir, PfamAFileName)
	biomartPath := filepath.Join(rawDir, PfamBiomartFileName)

	// Check if both files exist.
	if _, err := os.Stat(pfamAPath); err != nil {
		logger.Debug("PFAM pfamA.txt not found, skipping", zap.String("path", pfamAPath))
		return nil
	}
	if _, err := os.Stat(biomartPath); err != nil {
		logger.Debug("PFAM ensembl_biomart_pfam.txt not found, skipping", zap.String("path", biomartPath))
		return nil
	}

	start := time.Now()
	store, err := pfam.Load(pfamAPath, biomartPath)
	if err != nil {
		logger.Warn("could not load PFAM data", zap.Error(err))
		return nil
	}

	logger.Info("loaded PFAM domain data",
		zap.Int("domains", store.DomainCount()),
		zap.Int("transcripts", store.TranscriptCount()),
		zap.Duration("elapsed", time.Since(start)))

	return store
}

// loadPtmStore loads PTM data from the raw download directory.
func loadPtmStore(logger *zap.Logger, assembly string) *ptm.Store {
	rawDir := RawDir(assembly)
	if rawDir == "" {
		return nil
	}

	ptmPath := filepath.Join(rawDir, PtmFileName)
	if _, err := os.Stat(ptmPath); err != nil {
		logger.Debug("PTM data not found, skipping", zap.String("path", ptmPath))
		return nil
	}

	start := time.Now()
	store, err := ptm.Load(ptmPath)
	if err != nil {
		logger.Warn("could not load PTM data", zap.Error(err))
		return nil
	}

	logger.Info("loaded PTM data",
		zap.Int("transcripts", store.TranscriptCount()),
		zap.Duration("elapsed", time.Since(start)))

	return store
}

// loadUniprotStore loads UniProt transcript mapping from the raw download directory.
func loadUniprotStore(logger *zap.Logger, assembly string) *uniprot.Store {
	rawDir := RawDir(assembly)
	if rawDir == "" {
		return nil
	}

	uniprotPath := filepath.Join(rawDir, UniprotMappingFileName)
	if _, err := os.Stat(uniprotPath); err != nil {
		logger.Debug("UniProt mapping not found, skipping", zap.String("path", uniprotPath))
		return nil
	}

	start := time.Now()
	store, err := uniprot.Load(uniprotPath)
	if err != nil {
		logger.Warn("could not load UniProt mapping", zap.Error(err))
		return nil
	}

	logger.Info("loaded UniProt transcript mapping",
		zap.Int("mappings", store.Count()),
		zap.Duration("elapsed", time.Since(start)))

	return store
}
