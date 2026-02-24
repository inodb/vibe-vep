// Package main provides the vibe-vep command-line tool.
package main

import (
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/maf"
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
	case "download":
		return runDownload(args[1:])
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
  annotate    Annotate variants in a VCF or MAF file
  download    Download GENCODE annotation files
  help        Show this help message

Global Options:
  --version   Show version information

Examples:
  # Download GENCODE annotations (one-time setup)
  vibe-vep download --assembly GRCh38

  # Annotate a VCF file (uses GENCODE cache automatically)
  vibe-vep annotate input.vcf

  # Annotate a MAF file
  vibe-vep annotate input.maf

  # Validate MAF annotations against VEP predictions
  vibe-vep annotate --validate data_mutations.txt

For more information on a command, use:
  vibe-vep <command> --help
`)
}

func runAnnotate(args []string) int {
	fs := flag.NewFlagSet("annotate", flag.ExitOnError)

	var (
		assembly      string
		outputFormat  string
		outputFile    string
		canonicalOnly bool
		inputFormat   string
		validate      bool
		validateAll   bool
	)

	fs.StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")
	fs.StringVar(&outputFormat, "f", "tab", "Output format: tab, vcf")
	fs.StringVar(&outputFormat, "output-format", "tab", "Output format: tab, vcf")
	fs.StringVar(&outputFile, "o", "", "Output file (default: stdout)")
	fs.StringVar(&outputFile, "output", "", "Output file (default: stdout)")
	fs.BoolVar(&canonicalOnly, "canonical", false, "Only report canonical transcript annotations")
	fs.StringVar(&inputFormat, "input-format", "", "Input format: vcf, maf (auto-detected if not specified)")
	fs.BoolVar(&validate, "validate", false, "Validate MAF annotations against VEP predictions (MAF input only)")
	fs.BoolVar(&validateAll, "validate-all", false, "Show all variants in validation output (default: mismatches only)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Annotate variants in a VCF or MAF file with consequence predictions.

Usage:
  vibe-vep annotate [options] <input-file>

Arguments:
  <input-file>  Input VCF or MAF file (use '-' for stdin)

Options:
`)
		fs.PrintDefaults()
		fmt.Fprintf(os.Stderr, `
Examples:
  vibe-vep annotate input.vcf
  vibe-vep annotate input.maf
  vibe-vep annotate -f vcf -o output.vcf input.vcf
  vibe-vep annotate --validate data_mutations.txt
  cat input.vcf | vibe-vep annotate -
`)
	}

	if err := fs.Parse(args); err != nil {
		return ExitUsage
	}

	if fs.NArg() < 1 {
		fmt.Fprintf(os.Stderr, "Error: input file argument required\n\n")
		fs.Usage()
		return ExitUsage
	}

	inputPath := fs.Arg(0)

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
		fmt.Fprintf(os.Stderr, "Error: unknown input format %q\n", detectedFormat)
		fmt.Fprintf(os.Stderr, "Hint: Use --input-format to specify vcf or maf\n")
		return ExitError
	}

	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		if os.IsNotExist(err) {
			fmt.Fprintf(os.Stderr, "Hint: Check that the file path is correct\n")
		}
		return ExitError
	}
	defer parser.Close()

	// Load GENCODE cache
	gtfPath, fastaPath, canonicalPath, found := FindGENCODEFiles(assembly)
	if !found {
		fmt.Fprintf(os.Stderr, "Error: No GENCODE cache found for %s\n", assembly)
		fmt.Fprintf(os.Stderr, "Hint: Download GENCODE annotations with: vibe-vep download --assembly %s\n", assembly)
		return ExitError
	}

	fmt.Fprintf(os.Stderr, "Using GENCODE cache for %s\n", assembly)
	fmt.Fprintf(os.Stderr, "  GTF: %s\n", gtfPath)
	if fastaPath != "" {
		fmt.Fprintf(os.Stderr, "  FASTA: %s\n", fastaPath)
	}

	c := cache.New()
	loader := cache.NewGENCODELoader(gtfPath, fastaPath)

	// Load canonical transcript overrides if available
	if canonicalPath != "" {
		fmt.Fprintf(os.Stderr, "  Canonical overrides: %s\n", canonicalPath)
		overrides, err := cache.LoadCanonicalOverrides(canonicalPath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Warning: could not load canonical overrides: %v\n", err)
		} else {
			loader.SetCanonicalOverrides(overrides)
			fmt.Fprintf(os.Stderr, "  Loaded %d canonical overrides\n", len(overrides))
		}
	}

	if err := loader.Load(c); err != nil {
		fmt.Fprintf(os.Stderr, "Error loading GENCODE cache: %v\n", err)
		return ExitError
	}
	fmt.Fprintf(os.Stderr, "Loaded %d transcripts\n", c.TranscriptCount())
	transcriptCache := c

	// Create annotator
	ann := annotate.NewAnnotator(transcriptCache)
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

	// Validation mode for MAF files
	if validate {
		if detectedFormat != "maf" {
			fmt.Fprintf(os.Stderr, "Error: --validate requires MAF input format\n")
			return ExitError
		}

		mafParser, ok := parser.(*maf.Parser)
		if !ok {
			fmt.Fprintf(os.Stderr, "Error: --validate requires MAF parser\n")
			return ExitError
		}

		return runValidation(mafParser, ann, out, validateAll)
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

// runValidation runs validation mode comparing MAF annotations to VEP predictions.
func runValidation(parser *maf.Parser, ann *annotate.Annotator, out *os.File, showAll bool) int {
	valWriter := output.NewValidationWriter(out, showAll)

	if err := valWriter.WriteHeader(); err != nil {
		fmt.Fprintf(os.Stderr, "Error writing header: %v\n", err)
		return ExitError
	}

	for {
		v, mafAnn, err := parser.NextWithAnnotation()
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error reading variant: %v\n", err)
			return ExitError
		}
		if v == nil {
			break
		}

		// Get VEP annotations
		vepAnns, err := ann.Annotate(v)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Warning: failed to annotate %s:%d: %v\n", v.Chrom, v.Pos, err)
			continue
		}

		if err := valWriter.WriteComparison(v, mafAnn, vepAnns); err != nil {
			fmt.Fprintf(os.Stderr, "Error writing comparison: %v\n", err)
			return ExitError
		}
	}

	if err := valWriter.Flush(); err != nil {
		fmt.Fprintf(os.Stderr, "Error flushing output: %v\n", err)
		return ExitError
	}

	valWriter.WriteSummary(os.Stderr)

	return ExitSuccess
}


// detectInputFormat detects the input file format based on extension or content.
func detectInputFormat(path string) string {
	// Check by extension
	lowerPath := strings.ToLower(path)

	// Handle gzipped files
	if strings.HasSuffix(lowerPath, ".gz") {
		lowerPath = lowerPath[:len(lowerPath)-3]
	}

	if strings.HasSuffix(lowerPath, ".vcf") {
		return "vcf"
	}
	if strings.HasSuffix(lowerPath, ".maf") {
		return "maf"
	}

	// Check for cBioPortal MAF filenames
	baseName := filepath.Base(lowerPath)
	if baseName == "data_mutations.txt" || baseName == "data_mutations_extended.txt" {
		return "maf"
	}

	// For .txt files or stdin, try to detect from content
	if path == "-" {
		// Default to VCF for stdin
		return "vcf"
	}

	// Try to peek at the file to detect format
	file, err := os.Open(path)
	if err != nil {
		return "vcf" // Default to VCF
	}
	defer file.Close()

	buf := make([]byte, 512)
	n, err := file.Read(buf)
	if err != nil || n == 0 {
		return "vcf"
	}

	content := string(buf[:n])

	// Check for VCF header
	if strings.HasPrefix(content, "##fileformat=VCF") || strings.HasPrefix(content, "#CHROM") {
		return "vcf"
	}

	// Check for MAF header columns
	if strings.Contains(content, "Hugo_Symbol") && strings.Contains(content, "Chromosome") {
		return "maf"
	}

	// Default to VCF
	return "vcf"
}
