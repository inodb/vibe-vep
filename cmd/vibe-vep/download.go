package main

import (
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"strings"
	"time"

	"go.uber.org/zap"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

// GENCODE FTP base URL and per-assembly versions.
// VEP v111 (Ensembl 111) uses GENCODE v45 for GRCh38 and v19 for GRCh37.
const (
	gencodeBase = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"

	// GencodeVersionGRCh38 is the GENCODE release for GRCh38 (matches VEP v111).
	GencodeVersionGRCh38 = "v45"

	// GencodeVersionGRCh37 is the native GENCODE release for GRCh37 (Ensembl 75).
	GencodeVersionGRCh37 = "v19"
)

// GencodeVersionForAssembly returns the GENCODE version for the given assembly.
func GencodeVersionForAssembly(assembly string) string {
	if strings.EqualFold(assembly, "GRCh37") {
		return GencodeVersionGRCh37
	}
	return GencodeVersionGRCh38
}

// getGENCODEURLs returns the GTF and FASTA URLs for the given assembly.
// GRCh37 uses native GENCODE v19 (not liftover), GRCh38 uses GENCODE v45.
func getGENCODEURLs(assembly string) (gtfURL, fastaURL string) {
	ver := GencodeVersionForAssembly(assembly)
	release := strings.TrimPrefix(ver, "v")
	base := fmt.Sprintf("%s/release_%s", gencodeBase, release)

	gtfURL = fmt.Sprintf("%s/gencode.%s.annotation.gtf.gz", base, ver)
	fastaURL = fmt.Sprintf("%s/gencode.%s.pc_transcripts.fa.gz", base, ver)
	return
}

func newDownloadCmd(verbose *bool) *cobra.Command {
	var (
		assembly  string
		outputDir string
		gtfOnly   bool
	)

	cmd := &cobra.Command{
		Use:   "download",
		Short: "Download GENCODE annotation files",
		Long:  "Download GENCODE annotation files for variant annotation.",
		Example: `  # Download GRCh38 annotations (default)
  vibe-vep download

  # Download GRCh37 annotations
  vibe-vep download --assembly GRCh37

  # Download to a custom directory
  vibe-vep download --output /data/gencode`,
		Args: cobra.NoArgs,
		RunE: func(cmd *cobra.Command, args []string) error {
			logger, err := newLogger(*verbose)
			if err != nil {
				return fmt.Errorf("creating logger: %w", err)
			}
			defer logger.Sync()
			return runDownload(logger, assembly, outputDir, gtfOnly)
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")
	cmd.Flags().StringVar(&outputDir, "output", "", "Output directory (default: ~/.vibe-vep/)")
	cmd.Flags().BoolVar(&gtfOnly, "gtf-only", false, "Only download GTF annotations (skip FASTA sequences)")

	return cmd
}

func runDownload(logger *zap.Logger, assembly, outputDir string, gtfOnly bool) error {
	var err error
	assembly, err = normalizeAssembly(assembly)
	if err != nil {
		return err
	}

	// Determine output directory
	if outputDir == "" {
		home, err := os.UserHomeDir()
		if err != nil {
			return fmt.Errorf("cannot determine home directory: %w", err)
		}
		outputDir = filepath.Join(home, ".vibe-vep")
	}

	// Create assembly-specific subdirectory
	assemblyLower := strings.ToLower(assembly)
	destDir := filepath.Join(outputDir, assemblyLower)

	if err := os.MkdirAll(destDir, 0755); err != nil {
		return fmt.Errorf("cannot create directory %s: %w", destDir, err)
	}

	gtfURL, fastaURL := getGENCODEURLs(assembly)

	fmt.Printf("Downloading GENCODE %s annotations for %s...\n", GencodeVersionForAssembly(assembly), assembly)
	fmt.Printf("Destination: %s\n\n", destDir)

	// Download GTF
	gtfFile := filepath.Join(destDir, filepath.Base(gtfURL))
	if err := downloadFile(gtfURL, gtfFile); err != nil {
		return fmt.Errorf("downloading GTF: %w", err)
	}

	// Download FASTA (unless gtf-only)
	if !gtfOnly {
		fastaFile := filepath.Join(destDir, filepath.Base(fastaURL))
		if err := downloadFile(fastaURL, fastaFile); err != nil {
			return fmt.Errorf("downloading FASTA: %w", err)
		}
	}

	// Download Genome Nexus canonical transcript overrides
	canonicalURL := cache.CanonicalFileURL(assembly)
	canonicalFile := filepath.Join(destDir, cache.CanonicalFileName())
	if err := downloadFile(canonicalURL, canonicalFile); err != nil {
		logger.Warn("could not download canonical transcript overrides", zap.Error(err))
		// Non-fatal: tool still works without overrides
	}

	// Download AlphaMissense data if enabled in config
	if viper.GetBool("annotations.alphamissense") {
		amURL := getAlphaMissenseURL(assembly)
		amFile := filepath.Join(destDir, filepath.Base(amURL))
		fmt.Printf("\nAlphaMissense annotation enabled in config, downloading...\n")
		if err := downloadFile(amURL, amFile); err != nil {
			logger.Warn("could not download AlphaMissense data", zap.Error(err))
			// Non-fatal: tool still works without AlphaMissense
		}
	}

	// Download ClinVar data if enabled in config
	if viper.GetBool("annotations.clinvar") {
		cvURL := getClinVarURL(assembly)
		cvFile := filepath.Join(destDir, ClinVarFileName)
		fmt.Printf("\nClinVar annotation enabled in config, downloading...\n")
		if err := downloadFile(cvURL, cvFile); err != nil {
			logger.Warn("could not download ClinVar data", zap.Error(err))
		}
	}

	// Download gnomAD data if enabled in config
	if viper.GetBool("annotations.gnomad") {
		gnomadURL := getGnomadURL(assembly)
		gnomadFile := filepath.Join(destDir, GnomadFileName(assembly))
		fmt.Printf("\ngnomAD annotation enabled in config, downloading...\n")
		if err := downloadFile(gnomadURL, gnomadFile); err != nil {
			logger.Warn("could not download gnomAD data", zap.Error(err))
		}
	}

	// Download dbSNP data if enabled in config
	if viper.GetBool("annotations.dbsnp") {
		dbsnpURL := getDbSnpURL(assembly)
		dbsnpFile := filepath.Join(destDir, DbSnpFileName)
		fmt.Printf("\ndbSNP annotation enabled in config, downloading...\n")
		if err := downloadFile(dbsnpURL, dbsnpFile); err != nil {
			logger.Warn("could not download dbSNP data", zap.Error(err))
		}
	}

	fmt.Printf("\nDownload complete!\n")
	fmt.Printf("To annotate variants, run:\n")
	fmt.Printf("  vibe-vep annotate input.vcf\n")

	return nil
}

// downloadFile downloads a file from URL to the destination path with progress.
func downloadFile(url, destPath string) error {
	// Check if file already exists
	if info, err := os.Stat(destPath); err == nil {
		fmt.Printf("  %s already exists (%s), skipping\n", filepath.Base(destPath), formatSize(info.Size()))
		return nil
	}

	fmt.Printf("  Downloading %s...\n", filepath.Base(destPath))

	// Create HTTP client with timeout
	client := &http.Client{
		Timeout: 30 * time.Minute, // Long timeout for large files
	}

	resp, err := client.Get(url)
	if err != nil {
		return fmt.Errorf("HTTP request failed: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("HTTP error: %s", resp.Status)
	}

	// Create destination file
	tmpPath := destPath + ".tmp"
	f, err := os.Create(tmpPath)
	if err != nil {
		return fmt.Errorf("create file: %w", err)
	}

	// Copy with progress
	var downloaded int64
	contentLength := resp.ContentLength

	// Create a progress writer
	pw := &progressWriter{
		total:      contentLength,
		downloaded: &downloaded,
		lastPrint:  time.Now(),
	}

	_, err = io.Copy(f, io.TeeReader(resp.Body, pw))
	f.Close()

	if err != nil {
		os.Remove(tmpPath)
		return fmt.Errorf("download failed: %w", err)
	}

	// Rename temp file to final destination
	if err := os.Rename(tmpPath, destPath); err != nil {
		os.Remove(tmpPath)
		return fmt.Errorf("rename file: %w", err)
	}

	fmt.Printf("    Done: %s\n", formatSize(downloaded))
	return nil
}

// progressWriter tracks download progress.
type progressWriter struct {
	total      int64
	downloaded *int64
	lastPrint  time.Time
}

func (pw *progressWriter) Write(p []byte) (int, error) {
	n := len(p)
	*pw.downloaded += int64(n)

	// Print progress every second
	if time.Since(pw.lastPrint) > time.Second {
		if pw.total > 0 {
			pct := float64(*pw.downloaded) / float64(pw.total) * 100
			fmt.Printf("\r    Progress: %s / %s (%.1f%%)  ",
				formatSize(*pw.downloaded), formatSize(pw.total), pct)
		} else {
			fmt.Printf("\r    Progress: %s  ", formatSize(*pw.downloaded))
		}
		pw.lastPrint = time.Now()
	}

	return n, nil
}

// formatSize formats bytes as human-readable size.
func formatSize(bytes int64) string {
	const unit = 1024
	if bytes < unit {
		return fmt.Sprintf("%d B", bytes)
	}
	div, exp := int64(unit), 0
	for n := bytes / unit; n >= unit; n /= unit {
		div *= unit
		exp++
	}
	return fmt.Sprintf("%.1f %cB", float64(bytes)/float64(div), "KMGTPE"[exp])
}

// Annotation source data file names.
const (
	ClinVarFileName = "clinvar.vcf.gz"
	SignalFileName  = "signaldb_all_variants_frequencies.txt"
	DbNSFPDirName   = "dbnsfp"
	DbSnpFileName   = "dbsnp.vcf.gz"
)

// GnomadFileName returns the expected filename for the given assembly.
func GnomadFileName(assembly string) string {
	return filepath.Base(getGnomadURL(assembly))
}

// getClinVarURL returns the download URL for ClinVar data.
func getClinVarURL(assembly string) string {
	switch strings.ToUpper(assembly) {
	case "GRCH37":
		return "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"
	default:
		return "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
	}
}

// getAlphaMissenseURL returns the download URL for AlphaMissense data.
func getAlphaMissenseURL(assembly string) string {
	switch strings.ToUpper(assembly) {
	case "GRCH37":
		return "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz"
	default:
		return "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
	}
}

// AlphaMissenseFileName returns the expected filename for the assembly.
func AlphaMissenseFileName(assembly string) string {
	return filepath.Base(getAlphaMissenseURL(assembly))
}

// getDbSnpURL returns the download URL for dbSNP VCF.
func getDbSnpURL(assembly string) string {
	switch strings.ToUpper(assembly) {
	case "GRCH37":
		return "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz"
	default:
		return "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
	}
}

// getGnomadURL returns the download URL for gnomAD sites VCF.
// GRCh38 uses gnomAD v4.1 (joint genomes+exomes), GRCh37 uses gnomAD v2.1.1 (genomes).
func getGnomadURL(assembly string) string {
	switch strings.ToUpper(assembly) {
	case "GRCH37":
		return "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz"
	default:
		return "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.vcf.bgz"
	}
}

// DefaultGENCODEPath returns the default path for GENCODE cache files.
func DefaultGENCODEPath(assembly string) string {
	home, err := os.UserHomeDir()
	if err != nil {
		return ""
	}
	return filepath.Join(home, ".vibe-vep", strings.ToLower(assembly))
}

// FindGENCODEFiles looks for GENCODE files in the default location.
// Returns gtfPath, fastaPath, canonicalPath, and whether files were found.
func FindGENCODEFiles(assembly string) (gtfPath, fastaPath, canonicalPath string, found bool) {
	dir := DefaultGENCODEPath(assembly)
	if dir == "" {
		return "", "", "", false
	}

	// Look for GTF file (both assemblies use gencode.v*.annotation.gtf.gz now)
	matches, err := filepath.Glob(filepath.Join(dir, "gencode.v*.annotation.gtf.gz"))
	if err != nil || len(matches) == 0 {
		return "", "", "", false
	}
	gtfPath = matches[0]

	// Look for FASTA file

	matches, err = filepath.Glob(filepath.Join(dir, "gencode.v*.pc_transcripts.fa.gz"))
	if err == nil && len(matches) > 0 {
		fastaPath = matches[0]
	}

	// Look for canonical transcript overrides
	cPath := filepath.Join(dir, cache.CanonicalFileName())
	if _, err := os.Stat(cPath); err == nil {
		canonicalPath = cPath
	}

	return gtfPath, fastaPath, canonicalPath, true
}
