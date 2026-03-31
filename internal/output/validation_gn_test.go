package output_test

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
	"time"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/datasource/ensemblpred"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/output"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// TestValidateGN fetches annotations from the genome-nexus GRCh38 production API,
// re-annotates the same variants with vibe-vep, and compares key fields.
//
// Three-step process (each step caches its output in testdata/gn/):
//  1. Extract unique genomic locations from datahub_gdc MAF files
//  2. Fetch genome-nexus annotations via the REST API
//  3. Re-annotate with vibe-vep and compare field-by-field
//
// Skipped with -short. Run with:
//
//	go test ./internal/output/ -run TestValidateGN -v -count=1 -timeout 30m
func TestValidateGN(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping genome-nexus validation in short mode")
	}

	gnDir := findOrCreateGNDir(t)
	variantsFile := filepath.Join(gnDir, "test_variants.txt")
	annotationsFile := filepath.Join(gnDir, "gn_annotations.jsonl")
	reportFile := filepath.Join(gnDir, "comparison_report.md")
	jsonFile := filepath.Join(gnDir, "comparison_report.json")

	// Step 1: Extract variants from MAF files (skip if already done).
	if !fileExists(variantsFile) {
		extractGNVariants(t, variantsFile, 1000)
	} else {
		t.Logf("using cached variants: %s", variantsFile)
	}

	// Step 2: Fetch GN annotations (skip if already done).
	if !fileExists(annotationsFile) {
		fetchGNAnnotations(t, variantsFile, annotationsFile)
	} else {
		t.Logf("using cached GN annotations: %s", annotationsFile)
	}

	// Step 3: Re-annotate with vibe-vep and compare.
	compareGNAnnotations(t, annotationsFile, reportFile, jsonFile)
}

// extractGNVariants extracts unique genomic locations from datahub_gdc MAF files.
func extractGNVariants(t *testing.T, outputPath string, limit int) {
	t.Helper()

	studyDir := findStudyDir(t, "datahub_gdc")
	mafFiles, err := filepath.Glob(filepath.Join(studyDir, "*_data_mutations.txt"))
	if err != nil || len(mafFiles) == 0 {
		t.Skipf("no datahub_gdc MAF files found in %s", studyDir)
	}

	f, err := os.Create(outputPath)
	if err != nil {
		t.Fatalf("create output: %v", err)
	}
	defer f.Close()
	bw := bufio.NewWriter(f)

	seen := make(map[string]bool)
	count := 0

	for _, path := range mafFiles {
		if count >= limit {
			break
		}
		parser, err := maf.NewParser(path)
		if err != nil {
			t.Logf("WARN: skip %s: %v", path, err)
			continue
		}
		cols := parser.Columns()

		for {
			_, ann, err := parser.NextWithAnnotation()
			if err != nil {
				break
			}
			if ann == nil {
				break
			}
			fields := ann.RawFields
			chrom := strings.TrimPrefix(safeField(fields, cols.Chromosome), "chr")
			start := safeField(fields, cols.StartPosition)
			end := safeField(fields, cols.EndPosition)
			ref := safeField(fields, cols.ReferenceAllele)
			alt := safeField(fields, cols.TumorSeqAllele2)

			if chrom == "" || start == "" || alt == "" {
				continue
			}
			if end == "" {
				end = start
			}
			key := chrom + "," + start + "," + end + "," + ref + "," + alt
			if seen[key] {
				continue
			}
			seen[key] = true
			fmt.Fprintln(bw, key)
			count++
			if count >= limit {
				break
			}
		}
		parser.Close()
	}

	if err := bw.Flush(); err != nil {
		t.Fatalf("flush: %v", err)
	}
	t.Logf("extracted %d unique variants to %s", count, outputPath)
}

// fetchGNAnnotations fetches annotations from the genome-nexus production API.
func fetchGNAnnotations(t *testing.T, inputPath, outputPath string) {
	t.Helper()

	variants := readLines(t, inputPath)
	if len(variants) == 0 {
		t.Fatalf("no variants in %s", inputPath)
	}
	t.Logf("fetching GN annotations for %d variants", len(variants))

	f, err := os.Create(outputPath)
	if err != nil {
		t.Fatalf("create output: %v", err)
	}
	defer f.Close()
	bw := bufio.NewWriter(f)

	client := &http.Client{Timeout: 60 * time.Second}
	const (
		apiURL    = "https://grch38.genomenexus.org"
		batchSize = 200
		delay     = 500 * time.Millisecond
	)
	endpoint := apiURL + "/annotation/genomic"

	fetched := 0
	errors := 0

	for i := 0; i < len(variants); i += batchSize {
		end := i + batchSize
		if end > len(variants) {
			end = len(variants)
		}
		batch := variants[i:end]

		locations := make([]gnLocation, len(batch))
		for j, v := range batch {
			locations[j] = parseVariantLine(v)
		}

		body, _ := json.Marshal(locations)
		req, err := http.NewRequest("POST", endpoint, strings.NewReader(string(body)))
		if err != nil {
			t.Fatalf("create request: %v", err)
		}
		req.Header.Set("Content-Type", "application/json")
		req.Header.Set("Accept", "application/json")

		resp, err := client.Do(req)
		if err != nil {
			t.Logf("WARN: batch %d-%d failed: %v", i, end, err)
			errors += len(batch)
			time.Sleep(delay)
			continue
		}

		if resp.StatusCode != http.StatusOK {
			respBody, _ := io.ReadAll(resp.Body)
			resp.Body.Close()
			t.Logf("WARN: batch %d-%d: status %d: %s", i, end, resp.StatusCode, string(respBody))
			// Fall back to individual GET requests.
			for _, v := range batch {
				url := apiURL + "/annotation/genomic/" + v
				r, err := client.Get(url)
				if err != nil {
					errors++
					continue
				}
				b, _ := io.ReadAll(r.Body)
				r.Body.Close()
				if r.StatusCode == http.StatusOK {
					fmt.Fprintln(bw, string(b))
					fetched++
				} else {
					errors++
				}
			}
			time.Sleep(delay)
			continue
		}

		var annotations []json.RawMessage
		if err := json.NewDecoder(resp.Body).Decode(&annotations); err != nil {
			resp.Body.Close()
			t.Logf("WARN: batch %d-%d decode failed: %v", i, end, err)
			errors += len(batch)
			continue
		}
		resp.Body.Close()

		for _, ann := range annotations {
			fmt.Fprintln(bw, string(ann))
			fetched++
		}

		t.Logf("fetched %d/%d variants (%d errors)", fetched, len(variants), errors)
		if i+batchSize < len(variants) {
			time.Sleep(delay)
		}
	}

	if err := bw.Flush(); err != nil {
		t.Fatalf("flush: %v", err)
	}
	t.Logf("wrote %d annotations to %s (%d errors)", fetched, outputPath, errors)
}

// compareGNAnnotations re-annotates variants with vibe-vep and compares against GN.
func compareGNAnnotations(t *testing.T, gnFile, reportPath, jsonPath string) {
	t.Helper()

	c, _, _ := loadGENCODECache(t, "GRCh38")
	ann := annotate.NewAnnotator(c)

	// Load ensembl SIFT/PolyPhen source if available.
	sources := loadEnsemblPredSource(t)

	f, err := os.Open(gnFile)
	if err != nil {
		t.Fatalf("open GN file: %v", err)
	}
	defer f.Close()

	report := output.NewGNComparisonReport()
	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	processed := 0
	for scanner.Scan() {
		line := scanner.Bytes()
		if len(line) == 0 {
			continue
		}

		gnAnn, err := output.ParseGNAnnotation(line)
		if err != nil {
			report.AddComparison(output.GNVariantComparison{
				Error: fmt.Sprintf("parse error: %v", err),
			})
			continue
		}

		v := gnAnnotationToVariant(gnAnn)
		if v == nil {
			report.AddComparison(output.GNVariantComparison{
				VariantID: gnAnn.Variant,
				Error:     "could not determine variant coordinates",
			})
			continue
		}

		vepAnns, err := ann.Annotate(v)
		if err != nil {
			report.AddComparison(output.GNVariantComparison{
				VariantID: gnAnn.Variant,
				Error:     fmt.Sprintf("annotation error: %v", err),
			})
			continue
		}

		for _, src := range sources {
			src.Annotate(v, vepAnns)
		}

		cmp := output.CompareGNToVEP(gnAnn, vepAnns)
		report.AddComparison(cmp)
		processed++
	}

	if err := scanner.Err(); err != nil {
		t.Fatalf("read GN file: %v", err)
	}

	t.Logf("compared %d variants: %d full matches (%.1f%%), %d partial, %d errors",
		report.TotalVariants, report.FullMatches,
		pct(report.FullMatches, report.TotalVariants),
		report.PartialMatches, report.Errors)

	// Log per-field match rates.
	for _, field := range output.CompareGNFields {
		if stats, ok := report.FieldStats[field]; ok {
			t.Logf("  %-20s %5.1f%% match (%d/%d, %d both_empty)",
				field, stats.MatchRate(), stats.Matches, stats.Total-stats.BothEmpty, stats.BothEmpty)
		}
	}

	// Write markdown report.
	md := output.FormatGNReport(report)
	if err := os.WriteFile(reportPath, []byte(md), 0644); err != nil {
		t.Fatalf("write report: %v", err)
	}
	t.Logf("wrote report: %s", reportPath)

	// Write JSON report (without individual comparisons for size).
	reportCopy := *report
	reportCopy.Comparisons = nil
	jsonData, err := json.MarshalIndent(&reportCopy, "", "  ")
	if err != nil {
		t.Fatalf("marshal JSON: %v", err)
	}
	if err := os.WriteFile(jsonPath, jsonData, 0644); err != nil {
		t.Fatalf("write JSON: %v", err)
	}
	t.Logf("wrote JSON: %s", jsonPath)
}

// --- helpers ---

type gnLocation struct {
	Chromosome      string `json:"chromosome"`
	Start           int64  `json:"start"`
	End             int64  `json:"end"`
	ReferenceAllele string `json:"referenceAllele"`
	VariantAllele   string `json:"variantAllele"`
}

func parseVariantLine(line string) gnLocation {
	parts := strings.SplitN(line, ",", 5)
	gl := gnLocation{}
	if len(parts) >= 1 {
		gl.Chromosome = parts[0]
	}
	if len(parts) >= 2 {
		gl.Start, _ = strconv.ParseInt(parts[1], 10, 64)
	}
	if len(parts) >= 3 {
		gl.End, _ = strconv.ParseInt(parts[2], 10, 64)
	}
	if len(parts) >= 4 {
		gl.ReferenceAllele = parts[3]
	}
	if len(parts) >= 5 {
		gl.VariantAllele = parts[4]
	}
	return gl
}

// gnAnnotationToVariant converts a GN annotation to a vcf.Variant for re-annotation.
func gnAnnotationToVariant(gn *output.GNAnnotation) *vcf.Variant {
	if gn.OriginalVariantQuery != "" {
		parts := strings.SplitN(gn.OriginalVariantQuery, ",", 5)
		if len(parts) == 5 {
			start, err1 := strconv.ParseInt(parts[1], 10, 64)
			if err1 == nil {
				ref := parts[3]
				alt := parts[4]
				if ref == "-" {
					ref = ""
				}
				if alt == "-" {
					alt = ""
				}
				return &vcf.Variant{
					Chrom: strings.TrimPrefix(parts[0], "chr"),
					Pos:   start,
					Ref:   ref,
					Alt:   alt,
				}
			}
		}
	}
	if gn.SeqRegionName != "" && gn.Start > 0 && gn.AlleleString != "" {
		parts := strings.SplitN(gn.AlleleString, "/", 2)
		if len(parts) == 2 {
			ref := parts[0]
			alt := parts[1]
			if ref == "-" {
				ref = ""
			}
			if alt == "-" {
				alt = ""
			}
			return &vcf.Variant{
				Chrom: gn.SeqRegionName,
				Pos:   gn.Start,
				Ref:   ref,
				Alt:   alt,
			}
		}
	}
	return nil
}

func safeField(fields []string, i int) string {
	if i >= 0 && i < len(fields) {
		return fields[i]
	}
	return ""
}

func readLines(t *testing.T, path string) []string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	var lines []string
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line != "" && !strings.HasPrefix(line, "#") {
			lines = append(lines, line)
		}
	}
	if err := scanner.Err(); err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	return lines
}

func findOrCreateGNDir(t *testing.T) string {
	t.Helper()
	for _, base := range []string{"testdata", filepath.Join("..", "..", "testdata")} {
		dir := filepath.Join(base, "gn")
		if _, err := os.Stat(base); err == nil {
			if err := os.MkdirAll(dir, 0755); err != nil {
				t.Fatalf("create gn dir: %v", err)
			}
			return dir
		}
	}
	t.Fatalf("testdata directory not found")
	return ""
}

func fileExists(path string) bool {
	_, err := os.Stat(path)
	return err == nil
}

func pct(num, denom int) float64 {
	if denom == 0 {
		return 0
	}
	return float64(num) / float64(denom) * 100.0
}

// loadEnsemblPredSource tries to load the Ensembl SIFT/PolyPhen prediction source.
// Returns nil slice if not available (non-fatal).
func loadEnsemblPredSource(t *testing.T) []annotate.AnnotationSource {
	t.Helper()
	home, err := os.UserHomeDir()
	if err != nil {
		return nil
	}
	dbPath := filepath.Join(home, ".vibe-vep", "grch38", "ensembl_predictions.sqlite")
	store, err := ensemblpred.Open(dbPath)
	if err != nil {
		t.Logf("ensembl predictions not available (skipping SIFT/PolyPhen comparison): %v", err)
		return nil
	}
	t.Logf("loaded Ensembl SIFT/PolyPhen predictions")
	return []annotate.AnnotationSource{ensemblpred.NewSource(store)}
}
