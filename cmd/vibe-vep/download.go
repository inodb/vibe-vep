package main

import (
	"flag"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/inodb/vibe-vep/internal/cache"
)

// GENCODE FTP URLs
const (
	gencodeBaseURL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46"
	gencodeVersion = "v46"
)

// getGENCODEURLs returns the GTF and FASTA URLs for the given assembly.
func getGENCODEURLs(assembly string) (gtfURL, fastaURL string) {
	switch strings.ToUpper(assembly) {
	case "GRCH37":
		gtfURL = fmt.Sprintf("%s/GRCh37_mapping/gencode.%slift37.annotation.gtf.gz", gencodeBaseURL, gencodeVersion)
		fastaURL = fmt.Sprintf("%s/GRCh37_mapping/gencode.%slift37.pc_transcripts.fa.gz", gencodeBaseURL, gencodeVersion)
	case "GRCH38":
		gtfURL = fmt.Sprintf("%s/gencode.%s.annotation.gtf.gz", gencodeBaseURL, gencodeVersion)
		fastaURL = fmt.Sprintf("%s/gencode.%s.pc_transcripts.fa.gz", gencodeBaseURL, gencodeVersion)
	default:
		// Default to GRCh38
		gtfURL = fmt.Sprintf("%s/gencode.%s.annotation.gtf.gz", gencodeBaseURL, gencodeVersion)
		fastaURL = fmt.Sprintf("%s/gencode.%s.pc_transcripts.fa.gz", gencodeBaseURL, gencodeVersion)
	}
	return
}

func runDownload(args []string) int {
	fs := flag.NewFlagSet("download", flag.ExitOnError)

	var (
		assembly  string
		outputDir string
		gtfOnly   bool
	)

	fs.StringVar(&assembly, "assembly", "GRCh38", "Genome assembly: GRCh37 or GRCh38")
	fs.StringVar(&outputDir, "output", "", "Output directory (default: ~/.vibe-vep/)")
	fs.BoolVar(&gtfOnly, "gtf-only", false, "Only download GTF annotations (skip FASTA sequences)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Download GENCODE annotation files for variant annotation.

Usage:
  vibe-vep download [options]

Options:
`)
		fs.PrintDefaults()
		fmt.Fprintf(os.Stderr, `
Examples:
  # Download GRCh38 annotations (default)
  vibe-vep download

  # Download GRCh37 annotations
  vibe-vep download --assembly GRCh37

  # Download to a custom directory
  vibe-vep download --output /data/gencode

Files downloaded:
  - gencode.v46.annotation.gtf.gz (~50MB for GRCh38)
  - gencode.v46.pc_transcripts.fa.gz (~70MB for GRCh38)

After downloading, vibe-vep will automatically detect and use these files.
`)
	}

	if err := fs.Parse(args); err != nil {
		return ExitUsage
	}

	// Determine output directory
	if outputDir == "" {
		home, err := os.UserHomeDir()
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error: cannot determine home directory: %v\n", err)
			return ExitError
		}
		outputDir = filepath.Join(home, ".vibe-vep")
	}

	// Create assembly-specific subdirectory
	assemblyLower := strings.ToLower(assembly)
	destDir := filepath.Join(outputDir, assemblyLower)

	if err := os.MkdirAll(destDir, 0755); err != nil {
		fmt.Fprintf(os.Stderr, "Error: cannot create directory %s: %v\n", destDir, err)
		return ExitError
	}

	gtfURL, fastaURL := getGENCODEURLs(assembly)

	fmt.Printf("Downloading GENCODE %s annotations for %s...\n", gencodeVersion, assembly)
	fmt.Printf("Destination: %s\n\n", destDir)

	// Download GTF
	gtfFile := filepath.Join(destDir, filepath.Base(gtfURL))
	if err := downloadFile(gtfURL, gtfFile); err != nil {
		fmt.Fprintf(os.Stderr, "Error downloading GTF: %v\n", err)
		return ExitError
	}

	// Download FASTA (unless gtf-only)
	if !gtfOnly {
		fastaFile := filepath.Join(destDir, filepath.Base(fastaURL))
		if err := downloadFile(fastaURL, fastaFile); err != nil {
			fmt.Fprintf(os.Stderr, "Error downloading FASTA: %v\n", err)
			return ExitError
		}
	}

	// Download Genome Nexus canonical transcript overrides
	canonicalURL := cache.CanonicalFileURL(assembly)
	canonicalFile := filepath.Join(destDir, cache.CanonicalFileName())
	if err := downloadFile(canonicalURL, canonicalFile); err != nil {
		fmt.Fprintf(os.Stderr, "Warning: could not download canonical transcript overrides: %v\n", err)
		// Non-fatal: tool still works without overrides
	}

	fmt.Printf("\nDownload complete!\n")
	fmt.Printf("To annotate variants, run:\n")
	fmt.Printf("  vibe-vep annotate input.vcf\n")

	return ExitSuccess
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

	// Look for GTF file
	assemblyLower := strings.ToLower(assembly)
	var gtfPattern string
	if assemblyLower == "grch37" {
		gtfPattern = "gencode.v*lift37.annotation.gtf.gz"
	} else {
		gtfPattern = "gencode.v*.annotation.gtf.gz"
	}

	matches, err := filepath.Glob(filepath.Join(dir, gtfPattern))
	if err != nil || len(matches) == 0 {
		return "", "", "", false
	}
	gtfPath = matches[0]

	// Look for FASTA file
	var fastaPattern string
	if assemblyLower == "grch37" {
		fastaPattern = "gencode.v*lift37.pc_transcripts.fa.gz"
	} else {
		fastaPattern = "gencode.v*.pc_transcripts.fa.gz"
	}

	matches, err = filepath.Glob(filepath.Join(dir, fastaPattern))
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
