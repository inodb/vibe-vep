package main

import (
	"crypto/sha256"
	"encoding/hex"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"time"

	"go.uber.org/zap"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

// GENCODE ↔ Ensembl version mapping.
//
// Each GENCODE release corresponds to a specific Ensembl release since GENCODE
// gene models are produced by the Ensembl annotation team. This mapping is used
// to download matching Ensembl data (variation predictions, etc.).
//
// Reference: https://www.gencodegenes.org/human/releases.html
var gencodeEnsemblMap = map[string]int{
	"v45": 111, // Jul 2023
	"v46": 112, // May 2024
	"v47": 113, // Aug 2024
	"v19": 75,  // GRCh37 native (Dec 2013)
}

// GENCODE FTP base URL and per-assembly versions.
const (
	gencodeBase = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"

	// GencodeVersionGRCh38 is the GENCODE release for GRCh38.
	GencodeVersionGRCh38 = "v45"

	// GencodeVersionGRCh37 is the native GENCODE release for GRCh37.
	GencodeVersionGRCh37 = "v19"
)

// EnsemblReleaseForAssembly returns the Ensembl release number matching
// the GENCODE version used for the given assembly.
func EnsemblReleaseForAssembly(assembly string) int {
	ver := GencodeVersionForAssembly(assembly)
	if rel, ok := gencodeEnsemblMap[ver]; ok {
		return rel
	}
	return 111 // fallback
}

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
		Short: "Download GENCODE and annotation source data files",
		Long: `Download GENCODE annotation files and optional annotation source data.

Core files (always downloaded):
  GENCODE GTF + FASTA transcripts, canonical transcript overrides

Optional annotation sources (enabled via config):
  annotations.alphamissense  AlphaMissense pathogenicity scores (~643 MB)
  annotations.clinvar        ClinVar clinical significance (~182 MB)
  annotations.gnomad         gnomAD allele frequencies (~17 GB)
  annotations.sift           SIFT predictions via Ensembl (~4.1 GB, shared with polyphen)
  annotations.polyphen       PolyPhen-2 predictions via Ensembl (~4.1 GB, shared with sift)
  annotations.dbsnp          dbSNP RS identifiers (~17 GB)`,
		Example: `  # Download GRCh38 annotations (default)
  vibe-vep download

  # Download GRCh37 annotations
  vibe-vep download --assembly GRCh37

  # Enable and download SIFT + PolyPhen-2
  vibe-vep config set annotations.sift true
  vibe-vep config set annotations.polyphen true
  vibe-vep download`,
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
		outputDir = DataDir()
		if outputDir == "" {
			return fmt.Errorf("cannot determine data directory (set VIBE_VEP_DATA_DIR or HOME)")
		}
	}

	// Raw downloads go into {assembly}/raw/
	assemblyLower := strings.ToLower(assembly)
	assemblyDir := filepath.Join(outputDir, assemblyLower)
	rawDir := filepath.Join(assemblyDir, "raw")

	if err := os.MkdirAll(rawDir, 0755); err != nil {
		return fmt.Errorf("cannot create directory %s: %w", rawDir, err)
	}

	gtfURL, fastaURL := getGENCODEURLs(assembly)

	fmt.Printf("Downloading GENCODE %s annotations for %s...\n", GencodeVersionForAssembly(assembly), assembly)
	fmt.Printf("Destination: %s\n\n", rawDir)

	// checksums collects sha256 → filename for the checksum manifest.
	checksums := map[string]string{}
	addChecksum := func(path, sum string) {
		if sum != "" {
			checksums[filepath.Base(path)] = sum
		}
	}

	// Download GTF
	gtfFile := filepath.Join(rawDir, filepath.Base(gtfURL))
	sum, err := downloadFile(gtfURL, gtfFile)
	if err != nil {
		return fmt.Errorf("downloading GTF: %w", err)
	}
	addChecksum(gtfFile, sum)

	// Download FASTA (unless gtf-only)
	if !gtfOnly {
		fastaFile := filepath.Join(rawDir, filepath.Base(fastaURL))
		sum, err := downloadFile(fastaURL, fastaFile)
		if err != nil {
			return fmt.Errorf("downloading FASTA: %w", err)
		}
		addChecksum(fastaFile, sum)
	}

	// Download Genome Nexus canonical transcript overrides
	canonicalURL := cache.CanonicalFileURL(assembly)
	canonicalFile := filepath.Join(rawDir, cache.CanonicalFileName())
	if sum, err := downloadFile(canonicalURL, canonicalFile); err != nil {
		logger.Warn("could not download canonical transcript overrides", zap.Error(err))
	} else {
		addChecksum(canonicalFile, sum)
	}

	// Download AlphaMissense data if enabled in config
	if viper.GetBool("annotations.alphamissense") {
		amURL := getAlphaMissenseURL(assembly)
		amFile := filepath.Join(rawDir, filepath.Base(amURL))
		fmt.Printf("\nAlphaMissense annotation enabled in config, downloading...\n")
		if sum, err := downloadFile(amURL, amFile); err != nil {
			logger.Warn("could not download AlphaMissense data", zap.Error(err))
		} else {
			addChecksum(amFile, sum)
		}
	}

	// Download ClinVar data if enabled in config
	if viper.GetBool("annotations.clinvar") {
		cvURL := getClinVarURL(assembly)
		cvFile := filepath.Join(rawDir, ClinVarFileName)
		fmt.Printf("\nClinVar annotation enabled in config, downloading...\n")
		if sum, err := downloadFile(cvURL, cvFile); err != nil {
			logger.Warn("could not download ClinVar data", zap.Error(err))
		} else {
			addChecksum(cvFile, sum)
		}
	}

	// Download gnomAD data if enabled in config
	if viper.GetBool("annotations.gnomad") {
		gnomadURL := getGnomadURL(assembly)
		gnomadFile := filepath.Join(rawDir, GnomadFileName(assembly))
		fmt.Printf("\ngnomAD annotation enabled in config, downloading...\n")
		if sum, err := downloadFile(gnomadURL, gnomadFile); err != nil {
			logger.Warn("could not download gnomAD data", zap.Error(err))
		} else {
			addChecksum(gnomadFile, sum)
		}
	}

	// Download dbSNP data if enabled in config
	if viper.GetBool("annotations.dbsnp") {
		dbsnpURL := getDbSnpURL(assembly)
		dbsnpFile := filepath.Join(rawDir, DbSnpFileName)
		fmt.Printf("\ndbSNP annotation enabled in config, downloading...\n")
		if sum, err := downloadFile(dbsnpURL, dbsnpFile); err != nil {
			logger.Warn("could not download dbSNP data", zap.Error(err))
		} else {
			addChecksum(dbsnpFile, sum)
		}
	}

	// Download Ensembl SIFT/PolyPhen prediction data if enabled in config
	if viper.GetBool("annotations.sift") || viper.GetBool("annotations.polyphen") {
		baseURL := getEnsemblVariationBaseURL(assembly)
		fmt.Printf("\nSIFT/PolyPhen annotation enabled in config, downloading Ensembl prediction data...\n")
		md5File := filepath.Join(rawDir, EnsemblTranslationMD5Name)
		if sum, err := downloadFile(baseURL+"/"+EnsemblTranslationMD5Name, md5File); err != nil {
			logger.Warn("could not download Ensembl translation MD5 mapping", zap.Error(err))
		} else {
			addChecksum(md5File, sum)
		}
		predFile := filepath.Join(rawDir, EnsemblPredictionsName)
		if sum, err := downloadFile(baseURL+"/"+EnsemblPredictionsName, predFile); err != nil {
			logger.Warn("could not download Ensembl protein function predictions", zap.Error(err))
		} else {
			addChecksum(predFile, sum)
		}
	}

	// Write checksum manifest to assembly dir (not raw/, so it survives clean).
	if len(checksums) > 0 {
		if err := writeChecksums(filepath.Join(assemblyDir, "checksums.sha256"), checksums); err != nil {
			logger.Warn("could not write checksums", zap.Error(err))
		}
	}

	fmt.Printf("\nDownload complete!\n")
	fmt.Printf("To build indexes, run:\n")
	fmt.Printf("  vibe-vep prepare --assembly %s\n", assembly)

	return nil
}

// writeChecksums writes a sha256sum-compatible manifest file.
func writeChecksums(path string, checksums map[string]string) error {
	// Sort by filename for stable output.
	names := make([]string, 0, len(checksums))
	for name := range checksums {
		names = append(names, name)
	}
	sort.Strings(names)

	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	for _, name := range names {
		fmt.Fprintf(f, "%s  %s\n", checksums[name], name)
	}
	return nil
}

// downloadFile downloads a file from URL to the destination path with progress,
// computing a SHA-256 checksum along the way. Returns the hex-encoded checksum.
func downloadFile(url, destPath string) (string, error) {
	// Check if file already exists
	if info, err := os.Stat(destPath); err == nil {
		fmt.Printf("  %s already exists (%s), skipping\n", filepath.Base(destPath), formatSize(info.Size()))
		return checksumFile(destPath)
	}

	fmt.Printf("  Downloading %s...\n", filepath.Base(destPath))

	// Create HTTP client with timeout
	client := &http.Client{
		Timeout: 30 * time.Minute, // Long timeout for large files
	}

	resp, err := client.Get(url)
	if err != nil {
		return "", fmt.Errorf("HTTP request failed: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return "", fmt.Errorf("HTTP error: %s", resp.Status)
	}

	// Create destination file
	tmpPath := destPath + ".tmp"
	f, err := os.Create(tmpPath)
	if err != nil {
		return "", fmt.Errorf("create file: %w", err)
	}

	// Copy with progress + SHA-256
	var downloaded int64
	contentLength := resp.ContentLength

	pw := &progressWriter{
		total:      contentLength,
		downloaded: &downloaded,
		lastPrint:  time.Now(),
	}

	h := sha256.New()
	_, err = io.Copy(f, io.TeeReader(resp.Body, io.MultiWriter(pw, h)))
	f.Close()

	if err != nil {
		os.Remove(tmpPath)
		return "", fmt.Errorf("download failed: %w", err)
	}

	// Rename temp file to final destination
	if err := os.Rename(tmpPath, destPath); err != nil {
		os.Remove(tmpPath)
		return "", fmt.Errorf("rename file: %w", err)
	}

	checksum := hex.EncodeToString(h.Sum(nil))
	fmt.Printf("    Done: %s (sha256:%s)\n", formatSize(downloaded), checksum[:12])
	return checksum, nil
}

// checksumFile computes the SHA-256 checksum of an existing file.
func checksumFile(path string) (string, error) {
	f, err := os.Open(path)
	if err != nil {
		return "", err
	}
	defer f.Close()
	h := sha256.New()
	if _, err := io.Copy(h, f); err != nil {
		return "", err
	}
	return hex.EncodeToString(h.Sum(nil)), nil
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
	ClinVarFileName      = "clinvar.vcf.gz"
	SignalFileName       = "signaldb_all_variants_frequencies.txt"
	DbSnpFileName              = "dbsnp.vcf.gz"
	EnsemblPredDBName          = "ensembl_sift_polyphen.sqlite"
	EnsemblTranslationMD5Name  = "translation_md5.txt.gz"
	EnsemblPredictionsName     = "protein_function_predictions.txt.gz"
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

// getEnsemblVariationBaseURL returns the base URL for the Ensembl variation MySQL dump,
// using the Ensembl release that matches our GENCODE version.
func getEnsemblVariationBaseURL(assembly string) string {
	rel := EnsemblReleaseForAssembly(assembly)
	switch strings.ToUpper(assembly) {
	case "GRCH37":
		return fmt.Sprintf("https://ftp.ensembl.org/pub/release-%d/mysql/homo_sapiens_variation_%d_37", rel, rel)
	default:
		return fmt.Sprintf("https://ftp.ensembl.org/pub/release-%d/mysql/homo_sapiens_variation_%d_38", rel, rel)
	}
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

// DataDir returns the base data directory.
// Uses VIBE_VEP_DATA_DIR if set, otherwise defaults to ~/.vibe-vep.
func DataDir() string {
	if dir := os.Getenv("VIBE_VEP_DATA_DIR"); dir != "" {
		return dir
	}
	home, err := os.UserHomeDir()
	if err != nil {
		return ""
	}
	return filepath.Join(home, ".vibe-vep")
}

// DefaultGENCODEPath returns the default path for GENCODE cache files.
func DefaultGENCODEPath(assembly string) string {
	dir := DataDir()
	if dir == "" {
		return ""
	}
	return filepath.Join(dir, strings.ToLower(assembly))
}

// RawDir returns the path to the raw download subdirectory for an assembly.
func RawDir(assembly string) string {
	dir := DefaultGENCODEPath(assembly)
	if dir == "" {
		return ""
	}
	return filepath.Join(dir, "raw")
}

// FindGENCODEFiles looks for GENCODE files in the default location.
// Checks raw/ subdirectory first, falls back to flat layout for backward compatibility.
// Returns gtfPath, fastaPath, canonicalPath, and whether files were found.
func FindGENCODEFiles(assembly string) (gtfPath, fastaPath, canonicalPath string, found bool) {
	dir := DefaultGENCODEPath(assembly)
	if dir == "" {
		return "", "", "", false
	}

	// Try raw/ subdirectory first, fall back to flat layout.
	searchDirs := []string{filepath.Join(dir, "raw"), dir}

	for _, d := range searchDirs {
		matches, err := filepath.Glob(filepath.Join(d, "gencode.v*.annotation.gtf.gz"))
		if err != nil || len(matches) == 0 {
			continue
		}
		gtfPath = matches[0]

		matches, err = filepath.Glob(filepath.Join(d, "gencode.v*.pc_transcripts.fa.gz"))
		if err == nil && len(matches) > 0 {
			fastaPath = matches[0]
		}

		cPath := filepath.Join(d, cache.CanonicalFileName())
		if _, err := os.Stat(cPath); err == nil {
			canonicalPath = cPath
		}

		return gtfPath, fastaPath, canonicalPath, true
	}

	return "", "", "", false
}
