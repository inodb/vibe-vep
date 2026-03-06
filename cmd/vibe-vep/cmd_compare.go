package main

import (
	"fmt"
	"os"
	"strings"

	"github.com/inodb/vibe-vep/internal/output"
	"github.com/spf13/cobra"
)

func newCompareCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "compare",
		Short: "Compare variant annotations",
		Long: `Compare variant annotations between two files.

Subcommands:
  maf   Compare two MAF files side by side
  vcf   Compare two VCF files side by side`,
		Example: `  vibe-vep compare maf file1.maf file2.maf
  vibe-vep compare maf --categorize --columns Consequence,HGVSp_Short,HGVSc file1.maf file2.maf
  vibe-vep compare vcf file1.vcf file2.vcf`,
	}

	cmd.AddCommand(newCompareMAFCmd())
	cmd.AddCommand(newCompareVCFCmd())

	return cmd
}

// newCompareMAFCmd creates the "compare maf" subcommand for two-file MAF diff.
func newCompareMAFCmd() *cobra.Command {
	var (
		columns    string
		colMap     string
		all        bool
		maxDiffs   int
		categorize bool
	)

	cmd := &cobra.Command{
		Use:   "maf <left.maf> <right.maf>",
		Short: "Compare two MAF files side by side",
		Long: `Compare two MAF files column by column, showing differences and concordance rates.

Variants are matched by genomic position (Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2).
By default, all columns present in both files are compared.

With --categorize, known annotation columns (Consequence, HGVSp, HGVSc) get semantic
categorization (match, position_shift, transcript_model_change, etc.) instead of
simple string equality.`,
		Example: `  vibe-vep compare maf file1.maf file2.maf
  vibe-vep compare maf --columns Hugo_Symbol,Variant_Classification,HGVSp_Short file1.maf file2.maf
  vibe-vep compare maf --categorize --columns Consequence,HGVSp_Short,HGVSc file1.maf file2.maf
  vibe-vep compare maf --map 'Protein_Change=HGVSp_Short' file1.maf file2.maf
  vibe-vep compare maf --all file1.maf file2.maf
  vibe-vep compare maf --max-diffs 100 file1.maf file2.maf`,
		Args: cobra.ExactArgs(2),
		RunE: func(cmd *cobra.Command, args []string) error {
			return runCompareMAF(args[0], args[1], columns, colMap, all, maxDiffs, categorize)
		},
	}

	cmd.Flags().StringVar(&columns, "columns", "", "Columns to compare (comma-separated; default: all shared columns)")
	cmd.Flags().StringVar(&colMap, "map", "", "Map columns with different names: 'left_name=right_name,...'")
	cmd.Flags().BoolVar(&all, "all", false, "Show all rows, not just differences")
	cmd.Flags().IntVar(&maxDiffs, "max-diffs", 0, "Limit output to first N diffs (0 = unlimited)")
	cmd.Flags().BoolVar(&categorize, "categorize", false, "Use semantic categorization for annotation columns")

	return cmd
}

// newCompareVCFCmd creates the "compare vcf" subcommand for two-file VCF diff.
func newCompareVCFCmd() *cobra.Command {
	var (
		columns  string
		colMap   string
		all      bool
		maxDiffs int
	)

	cmd := &cobra.Command{
		Use:   "vcf <left.vcf> <right.vcf>",
		Short: "Compare two VCF files side by side",
		Long: `Compare two VCF files by INFO fields, showing differences and concordance rates.

Variants are matched by genomic position (CHROM, POS, REF, ALT).
By default, all shared INFO field keys are compared.`,
		Example: `  vibe-vep compare vcf file1.vcf file2.vcf
  vibe-vep compare vcf --columns CSQ,AF,DP file1.vcf file2.vcf
  vibe-vep compare vcf --all file1.vcf file2.vcf`,
		Args: cobra.ExactArgs(2),
		RunE: func(cmd *cobra.Command, args []string) error {
			return runCompareVCF(args[0], args[1], columns, colMap, all, maxDiffs)
		},
	}

	cmd.Flags().StringVar(&columns, "columns", "", "INFO keys to compare (comma-separated; default: all shared keys)")
	cmd.Flags().StringVar(&colMap, "map", "", "Map columns with different names: 'left_name=right_name,...'")
	cmd.Flags().BoolVar(&all, "all", false, "Show all rows, not just differences")
	cmd.Flags().IntVar(&maxDiffs, "max-diffs", 0, "Limit output to first N diffs (0 = unlimited)")

	return cmd
}

func runCompareMAF(leftPath, rightPath, columnsStr, mapStr string, showAll bool, maxDiffs int, categorize bool) error {
	// Parse --columns
	var columns []string
	if columnsStr != "" {
		for _, c := range strings.Split(columnsStr, ",") {
			c = strings.TrimSpace(c)
			if c != "" {
				columns = append(columns, c)
			}
		}
	}

	// Parse --map
	colMap, err := output.ParseColumnMap(mapStr)
	if err != nil {
		return fmt.Errorf("parsing --map: %w", err)
	}

	leftHeader, leftVariants, leftKeys, err := output.ReadMAFFile(leftPath)
	if err != nil {
		return fmt.Errorf("reading left file: %w", err)
	}

	rightHeader, rightVariants, rightKeys, err := output.ReadMAFFile(rightPath)
	if err != nil {
		return fmt.Errorf("reading right file: %w", err)
	}

	var cat *output.Categorizer
	if categorize {
		cat = &output.Categorizer{}
	}

	return output.CompareFiles(
		leftHeader, rightHeader,
		leftVariants, rightVariants,
		leftKeys, rightKeys,
		leftPath, rightPath,
		columns, colMap, showAll, maxDiffs,
		cat,
		os.Stdout, os.Stderr,
	)
}

func runCompareVCF(leftPath, rightPath, columnsStr, mapStr string, showAll bool, maxDiffs int) error {
	// Parse --columns
	var columns []string
	if columnsStr != "" {
		for _, c := range strings.Split(columnsStr, ",") {
			c = strings.TrimSpace(c)
			if c != "" {
				columns = append(columns, c)
			}
		}
	}

	// Parse --map
	colMap, err := output.ParseColumnMap(mapStr)
	if err != nil {
		return fmt.Errorf("parsing --map: %w", err)
	}

	leftHeader, leftVariants, leftKeys, err := output.ReadVCFFile(leftPath)
	if err != nil {
		return fmt.Errorf("reading left file: %w", err)
	}

	rightHeader, rightVariants, rightKeys, err := output.ReadVCFFile(rightPath)
	if err != nil {
		return fmt.Errorf("reading right file: %w", err)
	}

	return output.CompareFiles(
		leftHeader, rightHeader,
		leftVariants, rightVariants,
		leftKeys, rightKeys,
		leftPath, rightPath,
		columns, colMap, showAll, maxDiffs,
		nil,
		os.Stdout, os.Stderr,
	)
}
