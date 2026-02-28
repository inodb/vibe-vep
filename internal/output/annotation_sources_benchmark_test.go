package output_test

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strings"
	"testing"
	"time"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/datasource/alphamissense"
	"github.com/inodb/vibe-vep/internal/maf"
)

// TestAnnotationSourcesBenchmark runs AlphaMissense coverage and performance benchmarks
// against all TCGA MAF files and generates testdata/tcga/annotation_sources_report.md.
//
// Skipped with -short. Run with:
//
//	CGO_ENABLED=1 go test ./internal/output/ -run TestAnnotationSourcesBenchmark -v -count=1 -timeout 30m
func TestAnnotationSourcesBenchmark(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping annotation sources benchmark in short mode")
	}

	// Find TCGA MAF files.
	tcgaDir := findTCGADir(t)
	mafFiles, err := filepath.Glob(filepath.Join(tcgaDir, "*_data_mutations.txt"))
	if err != nil {
		t.Fatalf("glob MAF files: %v", err)
	}
	if len(mafFiles) == 0 {
		t.Skip("no TCGA MAF files found in", tcgaDir)
	}
	sort.Strings(mafFiles)

	// Load GENCODE cache.
	gtfPath, fastaPath, canonicalPath := findGENCODEFiles(t)
	c := cache.New()
	loader := cache.NewGENCODELoader(gtfPath, fastaPath)
	if canonicalPath != "" {
		overrides, err := cache.LoadCanonicalOverrides(canonicalPath)
		if err != nil {
			t.Logf("warning: could not load canonical overrides: %v", err)
		} else {
			loader.SetCanonicalOverrides(overrides)
		}
	}
	if err := loader.Load(c); err != nil {
		t.Fatalf("load GENCODE cache: %v", err)
	}
	t.Logf("loaded %d transcripts", c.TranscriptCount())

	// Open AlphaMissense database.
	home, err := os.UserHomeDir()
	if err != nil {
		t.Fatalf("get home dir: %v", err)
	}
	dbPath := filepath.Join(home, ".vibe-vep", "grch38", "annotations.duckdb")
	if _, err := os.Stat(dbPath); os.IsNotExist(err) {
		t.Skipf("AlphaMissense database not found at %s — run: vibe-vep prepare", dbPath)
	}
	amStore, err := alphamissense.Open(dbPath)
	if err != nil {
		t.Fatalf("open AlphaMissense: %v", err)
	}
	defer amStore.Close()
	if !amStore.Loaded() {
		t.Skip("AlphaMissense database is empty")
	}
	amCount, err := amStore.Count()
	if err != nil {
		t.Fatalf("count AlphaMissense rows: %v", err)
	}
	t.Logf("AlphaMissense variants in database: %d", amCount)

	var results []amStudyResult

	for _, mafFile := range mafFiles {
		name := studyName(mafFile)
		t.Run(name, func(t *testing.T) {
			r := benchmarkAMStudy(t, mafFile, c, amStore)
			r.name = name
			results = append(results, r)
		})
	}

	// Write report.
	reportPath := filepath.Join(tcgaDir, "annotation_sources_report.md")
	writeAnnotationSourcesReport(t, reportPath, results, c.TranscriptCount(), amCount)
}

// benchmarkAMStudy annotates a MAF file and benchmarks AlphaMissense lookups.
// Uses a two-pass approach: first annotate to find missense variants, then batch
// lookup all missense variants in a single DuckDB query.
func benchmarkAMStudy(t *testing.T, mafFile string, c *cache.Cache, amStore *alphamissense.Store) amStudyResult {
	t.Helper()

	parser, err := maf.NewParser(mafFile)
	if err != nil {
		t.Fatalf("open MAF: %v", err)
	}
	defer parser.Close()

	ann := annotate.NewAnnotator(c)

	// Pass 1: annotate all variants and collect missense lookup keys.
	type missenseKey struct {
		key alphamissense.LookupKey
	}
	var (
		totalVariants int
		missenseKeys  []alphamissense.LookupKey
		baseTime      time.Duration
	)
	// Deduplicate keys for efficient batch lookup.
	seen := make(map[alphamissense.LookupKey]bool)

	for {
		v, _, err := parser.NextWithAnnotation()
		if err != nil {
			t.Fatalf("read variant: %v", err)
		}
		if v == nil {
			break
		}
		totalVariants++

		annStart := time.Now()
		vepAnns, err := ann.Annotate(v)
		baseTime += time.Since(annStart)
		if err != nil {
			continue
		}

		for _, a := range vepAnns {
			if !isMissenseConsequence(a.Consequence) {
				continue
			}
			chrom := v.NormalizeChrom()
			if len(chrom) > 0 && chrom[0] != 'c' {
				chrom = "chr" + chrom
			}
			k := alphamissense.LookupKey{Chrom: chrom, Pos: v.Pos, Ref: v.Ref, Alt: v.Alt}
			if !seen[k] {
				seen[k] = true
				missenseKeys = append(missenseKeys, k)
			}
			break // one lookup per variant
		}
	}

	missenseCount := len(missenseKeys)

	// Pass 2: batch lookup all missense variants.
	lookupStart := time.Now()
	amResults, err := amStore.BatchLookup(missenseKeys)
	amLookupTime := time.Since(lookupStart)
	if err != nil {
		t.Fatalf("batch lookup: %v", err)
	}

	amHits := len(amResults)
	classCounts := make(map[string]int)
	for _, r := range amResults {
		classCounts[r.Class]++
	}

	coverage := 0.0
	if missenseCount > 0 {
		coverage = float64(amHits) / float64(missenseCount) * 100
	}
	t.Logf("%d variants, %d missense, %d AM hits (%.1f%%), batch lookup %s",
		totalVariants, missenseCount, amHits, coverage, amLookupTime.Round(time.Millisecond))

	return amStudyResult{
		variants:     totalVariants,
		missense:     missenseCount,
		amHits:       amHits,
		classCounts:  classCounts,
		baseTime:     baseTime,
		amLookupTime: amLookupTime,
	}
}

type amStudyResult struct {
	name         string
	variants     int
	missense     int
	amHits       int
	classCounts  map[string]int
	baseTime     time.Duration
	amLookupTime time.Duration
}

// isMissenseConsequence returns true if the consequence includes missense_variant.
func isMissenseConsequence(consequence string) bool {
	for rest := consequence; rest != ""; {
		term := rest
		if i := strings.IndexByte(rest, ','); i >= 0 {
			term = rest[:i]
			rest = rest[i+1:]
		} else {
			rest = ""
		}
		if term == "missense_variant" {
			return true
		}
	}
	return false
}

func writeAnnotationSourcesReport(t *testing.T, path string, results []amStudyResult, transcriptCount int, amCount int64) {
	t.Helper()

	var sb strings.Builder
	sb.WriteString("# Annotation Sources Report\n\n")
	sb.WriteString(fmt.Sprintf("Generated: %s  \n", time.Now().UTC().Format("2006-01-02 15:04 UTC")))
	sb.WriteString(fmt.Sprintf("GENCODE transcripts: %d  \n", transcriptCount))
	sb.WriteString(fmt.Sprintf("AlphaMissense variants in database: %d  \n", amCount))
	sb.WriteString(fmt.Sprintf("Workers: %d (GOMAXPROCS)\n\n", runtime.NumCPU()))

	// --- AlphaMissense Coverage ---
	sb.WriteString("## AlphaMissense Coverage\n\n")
	sb.WriteString("| Study | Variants | Missense | AM Hits | Coverage | likely_benign | ambiguous | likely_pathogenic |\n")
	sb.WriteString("|-------|----------|----------|---------|----------|---------------|-----------|-------------------|\n")

	var totVariants, totMissense, totHits int
	totClass := make(map[string]int)
	for _, r := range results {
		coverage := 0.0
		if r.missense > 0 {
			coverage = float64(r.amHits) / float64(r.missense) * 100
		}
		sb.WriteString(fmt.Sprintf("| %s | %d | %d | %d | %.1f%% | %d | %d | %d |\n",
			r.name, r.variants, r.missense, r.amHits, coverage,
			r.classCounts["likely_benign"],
			r.classCounts["ambiguous"],
			r.classCounts["likely_pathogenic"]))
		totVariants += r.variants
		totMissense += r.missense
		totHits += r.amHits
		for cls, n := range r.classCounts {
			totClass[cls] += n
		}
	}
	totCoverage := 0.0
	if totMissense > 0 {
		totCoverage = float64(totHits) / float64(totMissense) * 100
	}
	sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%d** | **%d** | **%.1f%%** | **%d** | **%d** | **%d** |\n\n",
		totVariants, totMissense, totHits, totCoverage,
		totClass["likely_benign"],
		totClass["ambiguous"],
		totClass["likely_pathogenic"]))

	// --- AlphaMissense Performance ---
	sb.WriteString("## AlphaMissense Performance\n\n")
	sb.WriteString("| Study | Variants | Base Time | AM Lookup Time | AM Overhead | Lookups/sec |\n")
	sb.WriteString("|-------|----------|-----------|----------------|-------------|-------------|\n")

	var totBaseTime, totAMLookupTime time.Duration
	for _, r := range results {
		overhead := 0.0
		if r.baseTime > 0 {
			overhead = float64(r.amLookupTime) / float64(r.baseTime) * 100
		}
		lookupsPerSec := 0.0
		if r.amLookupTime > 0 {
			lookupsPerSec = float64(r.missense) / r.amLookupTime.Seconds()
		}
		sb.WriteString(fmt.Sprintf("| %s | %d | %s | %s | %.1f%% | %.0f |\n",
			r.name, r.variants,
			r.baseTime.Round(time.Millisecond),
			r.amLookupTime.Round(time.Millisecond),
			overhead, lookupsPerSec))
		totBaseTime += r.baseTime
		totAMLookupTime += r.amLookupTime
	}
	totOverhead := 0.0
	if totBaseTime > 0 {
		totOverhead = float64(totAMLookupTime) / float64(totBaseTime) * 100
	}
	totLookupsPerSec := 0.0
	if totAMLookupTime > 0 {
		totLookupsPerSec = float64(totMissense) / totAMLookupTime.Seconds()
	}
	sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%s** | **%s** | **%.1f%%** | **%.0f** |\n",
		totVariants,
		totBaseTime.Round(time.Millisecond),
		totAMLookupTime.Round(time.Millisecond),
		totOverhead, totLookupsPerSec))

	if err := os.WriteFile(path, []byte(sb.String()), 0644); err != nil {
		t.Fatalf("write report: %v", err)
	}
	t.Logf("report written to %s", path)
}
