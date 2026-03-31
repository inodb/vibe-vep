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

// LoadBiomartCanonicals loads both MSK and Ensembl canonical transcripts
// and Entrez gene IDs from a Genome Nexus biomart TSV file.
func LoadBiomartCanonicals(path string) (mskcc, ensembl CanonicalOverrides, entrez GeneEntrezMap, err error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, nil, fmt.Errorf("open biomart canonical file: %w", err)
	}
	defer f.Close()

	return parseBiomartCanonicals(f)
}

// GeneEntrezMap maps gene symbol -> Entrez gene ID (as string).
type GeneEntrezMap map[string]string

// Biomart TSV column indices (0-indexed).
const (
	biomartGeneCol    = 0
	biomartEnsemblCol = 2  // ensembl_canonical_transcript
	biomartMSKCol     = 11 // mskcc_canonical_transcript
	biomartEntrezCol  = 27 // entrez_gene_id
)

// parseBiomartCanonicals parses a Genome Nexus biomart TSV, extracting
// the gene symbol (col 0), Ensembl canonical transcript (col 2),
// MSKCC canonical transcript (col 11), and Entrez gene ID (col 27).
func parseBiomartCanonicals(reader io.Reader) (mskcc, ensembl CanonicalOverrides, entrez GeneEntrezMap, err error) {
	mskcc = make(CanonicalOverrides)
	ensembl = make(CanonicalOverrides)
	entrez = make(GeneEntrezMap)
	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 1024*1024), 1024*1024) // 1MB line buffer for wide biomart files

	// Skip header line
	if !scanner.Scan() {
		return mskcc, ensembl, entrez, nil
	}

	for scanner.Scan() {
		line := scanner.Text()
		if line == "" {
			continue
		}

		fields := strings.Split(line, "\t")
		if len(fields) <= biomartMSKCol {
			continue
		}

		hgnc := fields[biomartGeneCol]
		if hgnc == "" {
			continue
		}

		// Ensembl canonical (col 2)
		if ensemblTx := fields[biomartEnsemblCol]; ensemblTx != "" && ensemblTx != "nan" {
			ensembl[hgnc] = stripVersion(ensemblTx)
		}

		// MSK canonical (col 11)
		if mskTx := fields[biomartMSKCol]; mskTx != "" && mskTx != "nan" {
			mskcc[hgnc] = stripVersion(mskTx)
		}

		// Entrez gene ID (col 27)
		if len(fields) > biomartEntrezCol {
			if eid := fields[biomartEntrezCol]; eid != "" && eid != "nan" {
				entrez[hgnc] = eid
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, nil, nil, fmt.Errorf("scan biomart canonicals: %w", err)
	}

	return mskcc, ensembl, entrez, nil
}

// LoadMSKCCOverrides loads MSKCC isoform overrides from the genome-nexus-importer format.
// Format: gene_name, refseq_id, enst_id, note (tab-separated, header row).
func LoadMSKCCOverrides(path string) (CanonicalOverrides, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open MSKCC overrides file: %w", err)
	}
	defer f.Close()

	return ParseMSKCCOverrides(f)
}

// ParseMSKCCOverrides parses MSKCC isoform override TSV content.
// Expected columns: gene_name (0), refseq_id (1), enst_id (2), note (3).
func ParseMSKCCOverrides(reader io.Reader) (CanonicalOverrides, error) {
	overrides := make(CanonicalOverrides)
	scanner := bufio.NewScanner(reader)

	// Skip header
	if !scanner.Scan() {
		return overrides, nil
	}

	for scanner.Scan() {
		line := scanner.Text()
		if line == "" {
			continue
		}

		fields := strings.Split(line, "\t")
		if len(fields) < 3 {
			continue
		}

		gene := fields[0]
		enst := fields[2]

		if gene == "" || enst == "" {
			continue
		}

		overrides[gene] = stripVersion(enst)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scan MSKCC overrides: %w", err)
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
