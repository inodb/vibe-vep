// Package main provides the vibe-vep command-line tool.
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/output"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Exit codes
const (
	ExitSuccess = 0
	ExitError   = 1
	ExitUsage   = 2
)

// Version information (set at build time)
var (
	version = "dev"
	commit  = "none"
	date    = "unknown"
)

func main() {
	os.Exit(run())
}

func run() int {
	// Global flags
	var showVersion bool
	flag.BoolVar(&showVersion, "version", false, "Show version information")

	// Parse global flags first
	flag.Parse()

	if showVersion {
		fmt.Printf("vibe-vep version %s (%s) built %s\n", version, commit, date)
		return ExitSuccess
	}

	// Check for subcommand
	args := flag.Args()
	if len(args) < 1 {
		printUsage()
		return ExitUsage
	}

	switch args[0] {
	case "annotate":
		return runAnnotate(args[1:])
	case "convert":
		return runConvert(args[1:])
	case "help":
		printUsage()
		return ExitSuccess
	default:
		fmt.Fprintf(os.Stderr, "Error: unknown command %q\n\n", args[0])
		printUsage()
		return ExitUsage
	}
}

func printUsage() {
	fmt.Fprintf(os.Stderr, `vibe-vep - Variant Effect Predictor

Usage:
  vibe-vep [options] <command> [arguments]

Commands:
  annotate    Annotate variants in a VCF file
  convert     Convert VEP cache to DuckDB format
  help        Show this help message

Global Options:
  --version   Show version information

Examples:
  # Annotate a VCF file
  vibe-vep annotate input.vcf

  # Read from stdin
  vibe-vep annotate -

  # Use VCF output format
  vibe-vep annotate -f vcf input.vcf

  # Use VEP Sereal cache directory
  vibe-vep annotate --cache-dir ~/.vep input.vcf

  # Use DuckDB cache file
  vibe-vep annotate --cache transcripts.duckdb input.vcf

  # Convert VEP cache to DuckDB
  vibe-vep convert --input ~/.vep --output transcripts.duckdb

For more information on a command, use:
  vibe-vep <command> --help
`)
}

func runAnnotate(args []string) int {
	fs := flag.NewFlagSet("annotate", flag.ExitOnError)

	var (
		cachePath     string
		cacheDir      string
		species       string
		assembly      string
		outputFormat  string
		outputFile    string
		canonicalOnly bool
	)

	fs.StringVar(&cachePath, "cache", "", "Cache file (DuckDB) or directory (Sereal)")
	fs.StringVar(&cacheDir, "cache-dir", defaultCacheDir(), "VEP cache directory (Sereal format)")
	fs.StringVar(&species, "species", "homo_sapiens", "Species name")
	fs.StringVar(&assembly, "assembly", "GRCh38", "Genome assembly")
	fs.StringVar(&outputFormat, "f", "tab", "Output format: tab, vcf")
	fs.StringVar(&outputFormat, "output-format", "tab", "Output format: tab, vcf")
	fs.StringVar(&outputFile, "o", "", "Output file (default: stdout)")
	fs.StringVar(&outputFile, "output", "", "Output file (default: stdout)")
	fs.BoolVar(&canonicalOnly, "canonical", false, "Only report canonical transcript annotations")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Annotate variants in a VCF file with consequence predictions.

Usage:
  vibe-vep annotate [options] <vcf-file>

Arguments:
  <vcf-file>  Input VCF file (use '-' for stdin)

Options:
`)
		fs.PrintDefaults()
		fmt.Fprintf(os.Stderr, `
Examples:
  vibe-vep annotate input.vcf
  vibe-vep annotate --cache transcripts.duckdb input.vcf
  vibe-vep annotate -f vcf -o output.vcf input.vcf
  cat input.vcf | vibe-vep annotate -
`)
	}

	if err := fs.Parse(args); err != nil {
		return ExitUsage
	}

	if fs.NArg() < 1 {
		fmt.Fprintf(os.Stderr, "Error: VCF file argument required\n\n")
		fs.Usage()
		return ExitUsage
	}

	vcfPath := fs.Arg(0)

	// Create VCF parser
	parser, err := vcf.NewParser(vcfPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		if os.IsNotExist(err) {
			fmt.Fprintf(os.Stderr, "Hint: Check that the file path is correct\n")
		}
		return ExitError
	}
	defer parser.Close()

	// Determine cache path
	effectiveCachePath := cachePath
	if effectiveCachePath == "" {
		effectiveCachePath = cacheDir
	}

	// Load cache (DuckDB or Sereal based on path)
	var c *cache.Cache
	var duckLoader *cache.DuckDBLoader
	if cache.IsDuckDB(effectiveCachePath) {
		// Use DuckDB cache
		duckLoader, err = cache.NewDuckDBLoader(effectiveCachePath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening DuckDB cache: %v\n", err)
			return ExitError
		}
		defer duckLoader.Close()

		c = cache.New()
		if err := duckLoader.LoadAll(c); err != nil {
			fmt.Fprintf(os.Stderr, "Error loading DuckDB cache: %v\n", err)
			return ExitError
		}
	} else {
		// Use Sereal cache
		c, err = loadCache(effectiveCachePath, species, assembly)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error loading cache: %v\n", err)
			fmt.Fprintf(os.Stderr, "Hint: Download VEP cache with: vep_install -a cf -s %s -y %s -c %s\n",
				species, assembly, effectiveCachePath)
			return ExitError
		}
	}

	// Create annotator
	ann := annotate.NewAnnotator(c)
	ann.SetCanonicalOnly(canonicalOnly)
	ann.SetWarnings(os.Stderr)

	// Create output writer
	var out *os.File
	if outputFile == "" {
		out = os.Stdout
	} else {
		out, err = os.Create(outputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
			return ExitError
		}
		defer out.Close()
	}

	var writer annotate.AnnotationWriter
	switch outputFormat {
	case "tab":
		writer = output.NewTabWriter(out)
	case "vcf":
		fmt.Fprintf(os.Stderr, "Error: VCF output format not yet implemented\n")
		return ExitError
	default:
		fmt.Fprintf(os.Stderr, "Error: unknown output format %q\n", outputFormat)
		return ExitError
	}

	// Write header
	if err := writer.WriteHeader(); err != nil {
		fmt.Fprintf(os.Stderr, "Error writing header: %v\n", err)
		return ExitError
	}

	// Annotate all variants
	if err := ann.AnnotateAll(parser, writer); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		return ExitError
	}

	return ExitSuccess
}

// defaultCacheDir returns the default VEP cache directory.
func defaultCacheDir() string {
	home, err := os.UserHomeDir()
	if err != nil {
		return ".vep"
	}
	return home + "/.vep"
}

// loadCache loads VEP cache data for annotation.
func loadCache(cacheDir, species, assembly string) (*cache.Cache, error) {
	c := cache.New()
	loader := cache.NewLoader(cacheDir, species, assembly)

	// Try to load all chromosomes
	err := loader.LoadAll(c)
	if err != nil {
		// Check if it's just a missing cache (warn but continue)
		cachePath := fmt.Sprintf("%s/%s/%s_%s", cacheDir, species, species, assembly)
		if _, statErr := os.Stat(cachePath); os.IsNotExist(statErr) {
			fmt.Fprintf(os.Stderr, "Warning: VEP cache not found at %s\n", cachePath)
			fmt.Fprintf(os.Stderr, "Warning: Annotations will be limited without cache data\n")
			return c, nil
		}
		return nil, err
	}

	if c.TranscriptCount() == 0 {
		fmt.Fprintf(os.Stderr, "Warning: No transcripts loaded from cache\n")
	}

	return c, nil
}
