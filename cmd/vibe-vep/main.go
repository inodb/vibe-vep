// Package main provides the vibe-vep command-line tool.
package main

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/output"
	"github.com/inodb/vibe-vep/internal/vcf"
	"github.com/spf13/cobra"
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
	rootCmd := &cobra.Command{
		Use:   "vibe-vep",
		Short: "Variant Effect Predictor",
		Long:  "vibe-vep - Variant Effect Predictor using GENCODE annotations",
		Version: fmt.Sprintf("%s (%s) built %s", version, commit, date),
	}

	rootCmd.AddCommand(newAnnotateCmd())
	rootCmd.AddCommand(newDownloadCmd())

	return rootCmd
}

func main() {
	if err := newRootCmd().Execute(); err != nil {
		os.Exit(ExitError)
	}
}

func newAnnotateCmd() *cobra.Command {
	var (
		assembly      string
		outputFormat  string
		outputFile    string
		canonicalOnly bool
		inputFormat   string
		validate      bool
		validateAll   bool
	)

	cmd := &cobra.Command{
		Use:   "annotate <input-file>",
		Short: "Annotate variants in a VCF or MAF file",
		Long:  "Annotate variants in a VCF or MAF file with consequence predictions.",
		Example: `  vibe-vep annotate input.vcf
  vibe-vep annotate input.maf
  vibe-vep annotate -f vcf -o output.vcf input.vcf
  vibe-vep annotate --validate data_mutations.txt
  cat input.vcf | vibe-vep annotate -`,
		Args: cobra.ExactArgs(1),
		RunE: func(cmd *cobra.Command, args []string) error {
			return runAnnotate(args[0], assembly, outputFormat, outputFile, canonicalOnly, inputFormat, validate, validateAll)
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")
	cmd.Flags().StringVarP(&outputFormat, "output-format", "f", "tab", "Output format: tab, vcf")
	cmd.Flags().StringVarP(&outputFile, "output", "o", "", "Output file (default: stdout)")
	cmd.Flags().BoolVar(&canonicalOnly, "canonical", false, "Only report canonical transcript annotations")
	cmd.Flags().StringVar(&inputFormat, "input-format", "", "Input format: vcf, maf (auto-detected if not specified)")
	cmd.Flags().BoolVar(&validate, "validate", false, "Validate MAF annotations against VEP predictions (MAF input only)")
	cmd.Flags().BoolVar(&validateAll, "validate-all", false, "Show all variants in validation output (default: mismatches only)")

	return cmd
}

func runAnnotate(inputPath, assembly, outputFormat, outputFile string, canonicalOnly bool, inputFormat string, validate, validateAll bool) error {
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
		return fmt.Errorf("unknown input format %q (use --input-format to specify vcf or maf)", detectedFormat)
	}

	if err != nil {
		if os.IsNotExist(err) {
			return fmt.Errorf("%w (check that the file path is correct)", err)
		}
		return err
	}
	defer parser.Close()

	// Load GENCODE cache
	gtfPath, fastaPath, canonicalPath, found := FindGENCODEFiles(assembly)
	if !found {
		return fmt.Errorf("no GENCODE cache found for %s\nHint: Download GENCODE annotations with: vibe-vep download --assembly %s", assembly, assembly)
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
		return fmt.Errorf("loading GENCODE cache: %w", err)
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
			return fmt.Errorf("creating output file: %w", err)
		}
		defer out.Close()
	}

	// Validation mode for MAF files
	if validate {
		if detectedFormat != "maf" {
			return fmt.Errorf("--validate requires MAF input format")
		}

		mafParser, ok := parser.(*maf.Parser)
		if !ok {
			return fmt.Errorf("--validate requires MAF parser")
		}

		return runValidation(mafParser, ann, out, validateAll)
	}

	var writer annotate.AnnotationWriter
	switch outputFormat {
	case "tab":
		writer = output.NewTabWriter(out)
	case "vcf":
		return fmt.Errorf("VCF output format not yet implemented")
	default:
		return fmt.Errorf("unknown output format %q", outputFormat)
	}

	// Write header
	if err := writer.WriteHeader(); err != nil {
		return fmt.Errorf("writing header: %w", err)
	}

	// Annotate all variants
	if err := ann.AnnotateAll(parser, writer); err != nil {
		return err
	}

	return nil
}

// runValidation runs validation mode comparing MAF annotations to VEP predictions.
func runValidation(parser *maf.Parser, ann *annotate.Annotator, out *os.File, showAll bool) error {
	valWriter := output.NewValidationWriter(out, showAll)

	if err := valWriter.WriteHeader(); err != nil {
		return fmt.Errorf("writing header: %w", err)
	}

	for {
		v, mafAnn, err := parser.NextWithAnnotation()
		if err != nil {
			return fmt.Errorf("reading variant: %w", err)
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
			return fmt.Errorf("writing comparison: %w", err)
		}
	}

	if err := valWriter.Flush(); err != nil {
		return fmt.Errorf("flushing output: %w", err)
	}

	valWriter.WriteSummary(os.Stderr)

	return nil
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
