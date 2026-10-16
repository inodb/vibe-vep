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
	biomartGeneCol        = 0
	biomartEnsemblCol     = 2  // ensembl_canonical_transcript
	biomartGenomeNexusCol = 5  // genome_nexus_canonical_transcript
	biomartUniprotCol     = 8  // uniprot_canonical_transcript
	biomartMSKCol         = 11 // mskcc_canonical_transcript
	biomartEntrezCol      = 31 // entrez_gene_id

	// biomartPrevSymbolsCol lists, comma separated, the symbols a gene used
	// to be called.
	biomartPrevSymbolsCol = 23 // previous_symbols
)

// parseBiomartCanonicals parses a Genome Nexus biomart TSV, extracting
// the gene symbol (col 0), Ensembl canonical transcript (col 2),
// MSKCC canonical transcript (col 11), and Entrez gene ID (col 31).
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

		// Entrez gene ID (col 31)
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

// LoadPreviousGeneSymbols maps a retired gene symbol to the approved one, so
// that a GENCODE release naming a gene WHSC1 can still be reported as NSD2.
//
// A retired symbol that some other gene has since taken as its own approved
// name is left alone -- ERF is now a gene in its own right even though it was
// once ETF1 -- and so is one that two genes both claim, since there would be
// no way to say which was meant.
func LoadPreviousGeneSymbols(path string) (map[string]string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open biomart canonicals: %w", err)
	}
	defer f.Close()
	return parsePreviousGeneSymbols(f)
}

func parsePreviousGeneSymbols(reader io.Reader) (map[string]string, error) {
	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)

	// Header.
	if !scanner.Scan() {
		return map[string]string{}, nil
	}

	type row struct {
		official string
		previous []string
	}
	var rows []row
	approved := make(map[string]bool)

	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) <= biomartPrevSymbolsCol {
			continue
		}
		official := strings.TrimSpace(fields[biomartGeneCol])
		if official == "" {
			continue
		}
		approved[official] = true

		var previous []string
		for _, p := range strings.Split(fields[biomartPrevSymbolsCol], ",") {
			if p = strings.TrimSpace(p); p != "" {
				previous = append(previous, p)
			}
		}
		if len(previous) > 0 {
			rows = append(rows, row{official: official, previous: previous})
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scan biomart canonicals: %w", err)
	}

	candidates := make(map[string]map[string]bool)
	for _, r := range rows {
		for _, prev := range r.previous {
			if approved[prev] {
				continue
			}
			if candidates[prev] == nil {
				candidates[prev] = make(map[string]bool)
			}
			candidates[prev][r.official] = true
		}
	}

	symbols := make(map[string]string, len(candidates))
	for prev, official := range candidates {
		if len(official) != 1 {
			continue
		}
		for name := range official {
			symbols[prev] = name
		}
	}
	return symbols, nil
}

// canonicalSourceCols names the isoform override tables the biomart file
// carries, by the name the API asks for them by.
var canonicalSourceCols = map[string]int{
	IsoformSourceEnsembl:     biomartEnsemblCol,
	IsoformSourceGenomeNexus: biomartGenomeNexusCol,
	IsoformSourceUniprot:     biomartUniprotCol,
	IsoformSourceMSKCC:       biomartMSKCol,
}

// The isoform override sources the API recognises.
const (
	IsoformSourceEnsembl     = "ensembl"
	IsoformSourceGenomeNexus = "genome_nexus"
	IsoformSourceUniprot     = "uniprot"
	IsoformSourceMSKCC       = "mskcc"
)

// NormalizeIsoformSource resolves what an API caller asked for to a table
// that exists. An unnamed source, and anything not recognised, means MSKCC.
func NormalizeIsoformSource(source string) string {
	normalized := strings.ToLower(strings.TrimSpace(source))
	switch normalized {
	case IsoformSourceEnsembl, IsoformSourceGenomeNexus, IsoformSourceUniprot:
		return normalized
	default:
		return IsoformSourceMSKCC
	}
}

// LoadCanonicalOverridesBySource reads every isoform override table the
// biomart file carries, keyed by source name, so that a request naming one of
// them gets that table rather than a silent fallback to another.
func LoadCanonicalOverridesBySource(path string) (map[string]CanonicalOverrides, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open biomart canonicals: %w", err)
	}
	defer f.Close()
	return parseCanonicalOverridesBySource(f)
}

func parseCanonicalOverridesBySource(reader io.Reader) (map[string]CanonicalOverrides, error) {
	bySource := make(map[string]CanonicalOverrides, len(canonicalSourceCols))
	widest := 0
	for source, col := range canonicalSourceCols {
		bySource[source] = make(CanonicalOverrides)
		if col > widest {
			widest = col
		}
	}

	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)

	// Header.
	if !scanner.Scan() {
		return bySource, nil
	}

	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) <= widest {
			continue
		}
		gene := strings.TrimSpace(fields[biomartGeneCol])
		if gene == "" {
			continue
		}
		for source, col := range canonicalSourceCols {
			txID := strings.TrimSpace(fields[col])
			if txID == "" || txID == "nan" {
				continue
			}
			bySource[source][gene] = stripVersion(txID)
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scan biomart canonicals: %w", err)
	}
	return bySource, nil
}

// GeneMetadata is what the biomart HGNC table says about a gene beyond its
// canonical transcript: the identifiers and alternative names the /ensembl
// gene endpoints and VEP's hgnc_id field report.
type GeneMetadata struct {
	// HGNCID is the numeric part of the HGNC identifier. The table writes it
	// as "HGNC:1097"; VEP reports it, and so genome-nexus does, as "1097".
	HGNCID string

	// Synonyms and PreviousSymbols are the other names a gene answers to,
	// in the order the table lists them.
	Synonyms        []string
	PreviousSymbols []string
}

const (
	biomartHGNCIDCol   = 14 // hgnc_id
	biomartSynonymsCol = 21 // synonyms
)

// LoadGeneMetadata reads the biomart HGNC table, keyed by approved symbol.
func LoadGeneMetadata(path string) (map[string]GeneMetadata, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open biomart canonicals: %w", err)
	}
	defer f.Close()
	return parseGeneMetadata(f)
}

// splitNameList reads one of the table's comma-separated name lists.
func splitNameList(s string) []string {
	var out []string
	for _, name := range strings.Split(s, ",") {
		if name = strings.TrimSpace(name); name != "" {
			out = append(out, name)
		}
	}
	return out
}

func parseGeneMetadata(reader io.Reader) (map[string]GeneMetadata, error) {
	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)

	// Header.
	if !scanner.Scan() {
		return map[string]GeneMetadata{}, nil
	}

	meta := make(map[string]GeneMetadata)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) <= biomartPrevSymbolsCol {
			continue
		}
		symbol := strings.TrimSpace(fields[biomartGeneCol])
		if symbol == "" {
			continue
		}
		meta[symbol] = GeneMetadata{
			HGNCID:          strings.TrimPrefix(strings.TrimSpace(fields[biomartHGNCIDCol]), "HGNC:"),
			Synonyms:        splitNameList(fields[biomartSynonymsCol]),
			PreviousSymbols: splitNameList(fields[biomartPrevSymbolsCol]),
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scan biomart canonicals: %w", err)
	}
	return meta, nil
}
