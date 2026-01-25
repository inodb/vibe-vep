package main

import (
	"flag"
	"fmt"
	"os"
	"path/filepath"

	"github.com/inodb/vibe-vep/internal/cache"
)

func runConvert(args []string) int {
	fs := flag.NewFlagSet("convert", flag.ExitOnError)

	var (
		inputPath  string
		outputPath string
		species    string
		assembly   string
		chrom      string
	)

	fs.StringVar(&inputPath, "input", "", "Input VEP cache directory or DuckDB file")
	fs.StringVar(&inputPath, "i", "", "Input VEP cache directory or DuckDB file (shorthand)")
	fs.StringVar(&outputPath, "output", "", "Output DuckDB file path")
	fs.StringVar(&outputPath, "o", "", "Output DuckDB file path (shorthand)")
	fs.StringVar(&species, "species", "homo_sapiens", "Species name")
	fs.StringVar(&assembly, "assembly", "GRCh38", "Genome assembly")
	fs.StringVar(&chrom, "chrom", "", "Only convert specific chromosome (optional)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Convert VEP Sereal cache to DuckDB format.

This command converts the VEP cache files (Sereal format) to a DuckDB database
for faster queries and S3 support.

Usage:
  vibe-vep convert [options]

Options:
`)
		fs.PrintDefaults()
		fmt.Fprintf(os.Stderr, `
Examples:
  # Convert full cache to DuckDB
  vibe-vep convert --input ~/.vep --output transcripts.duckdb

  # Convert specific chromosome
  vibe-vep convert --input ~/.vep --output chr12.duckdb --chrom 12

  # Convert with species/assembly
  vibe-vep convert -i ~/.vep -o cache.duckdb --species homo_sapiens --assembly GRCh38
`)
	}

	if err := fs.Parse(args); err != nil {
		return ExitUsage
	}

	// Validate required arguments
	if inputPath == "" {
		fmt.Fprintf(os.Stderr, "Error: --input is required\n\n")
		fs.Usage()
		return ExitUsage
	}
	if outputPath == "" {
		fmt.Fprintf(os.Stderr, "Error: --output is required\n\n")
		fs.Usage()
		return ExitUsage
	}

	// Ensure output has .duckdb extension
	if filepath.Ext(outputPath) != ".duckdb" && filepath.Ext(outputPath) != ".db" {
		outputPath = outputPath + ".duckdb"
	}

	// Remove existing output file if it exists
	if _, err := os.Stat(outputPath); err == nil {
		if err := os.Remove(outputPath); err != nil {
			fmt.Fprintf(os.Stderr, "Error removing existing file: %v\n", err)
			return ExitError
		}
	}

	fmt.Fprintf(os.Stderr, "Converting VEP cache to DuckDB...\n")
	fmt.Fprintf(os.Stderr, "  Input:  %s\n", inputPath)
	fmt.Fprintf(os.Stderr, "  Output: %s\n", outputPath)
	fmt.Fprintf(os.Stderr, "  Species: %s\n", species)
	fmt.Fprintf(os.Stderr, "  Assembly: %s\n", assembly)
	if chrom != "" {
		fmt.Fprintf(os.Stderr, "  Chromosome: %s\n", chrom)
	}

	// Load source cache
	srcCache := cache.New()
	srcLoader := cache.NewLoader(inputPath, species, assembly)

	// Read cache info if available
	info, err := srcLoader.ReadInfo()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Warning: Could not read cache info: %v\n", err)
	} else {
		fmt.Fprintf(os.Stderr, "  VEP Version: %s\n", info.Version)
		fmt.Fprintf(os.Stderr, "  Available chromosomes: %v\n", info.Chromosomes)
	}

	// Load transcripts
	if chrom != "" {
		fmt.Fprintf(os.Stderr, "\nLoading chromosome %s...\n", chrom)
		if err := srcLoader.Load(srcCache, chrom); err != nil {
			fmt.Fprintf(os.Stderr, "Error loading chromosome: %v\n", err)
			return ExitError
		}
	} else {
		fmt.Fprintf(os.Stderr, "\nLoading all chromosomes...\n")
		if err := srcLoader.LoadAll(srcCache); err != nil {
			fmt.Fprintf(os.Stderr, "Error loading cache: %v\n", err)
			return ExitError
		}
	}

	fmt.Fprintf(os.Stderr, "Loaded %d transcripts from %d chromosomes\n",
		srcCache.TranscriptCount(), len(srcCache.Chromosomes()))

	if srcCache.TranscriptCount() == 0 {
		fmt.Fprintf(os.Stderr, "Warning: No transcripts loaded\n")
		return ExitSuccess
	}

	// Create DuckDB database
	duckLoader, err := cache.NewDuckDBLoader(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating DuckDB: %v\n", err)
		return ExitError
	}
	defer duckLoader.Close()

	if err := duckLoader.CreateSchema(); err != nil {
		fmt.Fprintf(os.Stderr, "Error creating schema: %v\n", err)
		return ExitError
	}

	// Insert transcripts
	fmt.Fprintf(os.Stderr, "Writing transcripts to DuckDB...\n")
	var insertCount int
	for _, c := range srcCache.Chromosomes() {
		for _, t := range srcCache.FindTranscriptsByChrom(c) {
			if err := duckLoader.InsertTranscript(t); err != nil {
				fmt.Fprintf(os.Stderr, "Error inserting transcript %s: %v\n", t.ID, err)
				return ExitError
			}
			insertCount++
			if insertCount%1000 == 0 {
				fmt.Fprintf(os.Stderr, "  Inserted %d transcripts...\n", insertCount)
			}
		}
	}

	// Verify count
	finalCount, err := duckLoader.TranscriptCount()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error verifying count: %v\n", err)
		return ExitError
	}

	// Get file size
	stat, err := os.Stat(outputPath)
	var sizeStr string
	if err == nil {
		sizeMB := float64(stat.Size()) / (1024 * 1024)
		sizeStr = fmt.Sprintf("%.2f MB", sizeMB)
	} else {
		sizeStr = "unknown"
	}

	fmt.Fprintf(os.Stderr, "\nConversion complete!\n")
	fmt.Fprintf(os.Stderr, "  Transcripts: %d\n", finalCount)
	fmt.Fprintf(os.Stderr, "  Output size: %s\n", sizeStr)
	fmt.Fprintf(os.Stderr, "  Output file: %s\n", outputPath)

	return ExitSuccess
}
