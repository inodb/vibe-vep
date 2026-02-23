// Package cache provides VEP cache loading functionality.
package cache

import (
	"bufio"
	"fmt"
	"io"
	"net/http"
	"os"
	"strings"
	"time"
)

// CanonicalOverrides maps gene symbol -> canonical transcript ID.
type CanonicalOverrides map[string]string

// Genome Nexus canonical transcript file URLs.
const (
	canonicalFileGRCh38 = "https://raw.githubusercontent.com/genome-nexus/genome-nexus-importer/master/data/grch38_ensembl95/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt"
	canonicalFileGRCh37 = "https://raw.githubusercontent.com/genome-nexus/genome-nexus-importer/master/data/grch37_ensembl92/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt"
	canonicalFileName   = "ensembl_biomart_canonical_transcripts_per_hgnc.txt"
)

// CanonicalFileURL returns the URL for the canonical transcript file for the given assembly.
func CanonicalFileURL(assembly string) string {
	if strings.EqualFold(assembly, "GRCh37") {
		return canonicalFileGRCh37
	}
	return canonicalFileGRCh38
}

// CanonicalFileName returns the filename for the canonical transcript file.
func CanonicalFileName() string {
	return canonicalFileName
}

// LoadCanonicalOverrides loads canonical transcript overrides from a Genome Nexus TSV file.
// The file has columns: hgnc_symbol (col 0) and genome_nexus_canonical_transcript (col 4).
func LoadCanonicalOverrides(path string) (CanonicalOverrides, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open canonical overrides file: %w", err)
	}
	defer f.Close()

	return parseCanonicalOverrides(f)
}

// parseCanonicalOverrides parses the TSV content.
func parseCanonicalOverrides(reader io.Reader) (CanonicalOverrides, error) {
	overrides := make(CanonicalOverrides)
	scanner := bufio.NewScanner(reader)

	// Skip header line
	if !scanner.Scan() {
		return overrides, nil
	}

	for scanner.Scan() {
		line := scanner.Text()
		if line == "" {
			continue
		}

		fields := strings.Split(line, "\t")
		if len(fields) < 5 {
			continue
		}

		hgnc := fields[0]
		transcript := fields[4]

		if hgnc == "" || transcript == "" || transcript == "nan" {
			continue
		}

		// Strip version from transcript ID
		overrides[hgnc] = stripVersion(transcript)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scan canonical overrides: %w", err)
	}

	return overrides, nil
}

// DownloadCanonicalOverrides downloads the canonical transcript file to the given path.
func DownloadCanonicalOverrides(assembly, destPath string) error {
	url := CanonicalFileURL(assembly)

	client := &http.Client{Timeout: 5 * time.Minute}
	resp, err := client.Get(url)
	if err != nil {
		return fmt.Errorf("download canonical overrides: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("download canonical overrides: HTTP %s", resp.Status)
	}

	f, err := os.Create(destPath + ".tmp")
	if err != nil {
		return fmt.Errorf("create file: %w", err)
	}

	if _, err := io.Copy(f, resp.Body); err != nil {
		f.Close()
		os.Remove(destPath + ".tmp")
		return fmt.Errorf("write canonical overrides: %w", err)
	}
	f.Close()

	if err := os.Rename(destPath+".tmp", destPath); err != nil {
		os.Remove(destPath + ".tmp")
		return fmt.Errorf("rename canonical overrides: %w", err)
	}

	return nil
}
