package output_test

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
	"time"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/output"
)

// TestValidationBenchmark runs comparison against all TCGA MAF files and generates
// a markdown report at testdata/tcga/validation_report.md.
//
// Skipped with -short (large files, not for CI). Run with:
//
//	go test ./internal/output/ -run TestValidationBenchmark -v -count=1
func TestValidationBenchmark(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping validation benchmark in short mode")
	}

	// Find TCGA MAF files
	tcgaDir := findTCGADir(t)
	mafFiles, err := filepath.Glob(filepath.Join(tcgaDir, "*_data_mutations.txt"))
	if err != nil {
		t.Fatalf("glob MAF files: %v", err)
	}
	if len(mafFiles) == 0 {
		t.Skip("no TCGA MAF files found in", tcgaDir)
	}
	sort.Strings(mafFiles)

	// Load GENCODE cache once
	gtfPath, fastaPath, canonicalPath := findGENCODEFiles(t)
	c := cache.New()
	loader := cache.NewGENCODELoader(gtfPath, fastaPath)
	if canonicalPath != "" {
		overrides, err := cache.LoadCanonicalOverridesWithSource(canonicalPath, "mskcc")
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

	allCols := map[string]bool{"consequence": true, "hgvsp": true, "hgvsc": true}
	var results []studyResult

	for _, mafFile := range mafFiles {
		name := studyName(mafFile)
		t.Run(name, func(t *testing.T) {
			parser, err := maf.NewParser(mafFile)
			if err != nil {
				t.Fatalf("open MAF: %v", err)
			}
			defer parser.Close()

			ann := annotate.NewAnnotator(c)
			cmpWriter := output.NewCompareWriter(os.Stderr, allCols, false)

			start := time.Now()
			for {
				v, mafAnn, err := parser.NextWithAnnotation()
				if err != nil {
					t.Fatalf("read variant: %v", err)
				}
				if v == nil {
					break
				}
				vepAnns, err := ann.Annotate(v)
				if err != nil {
					continue
				}
				if err := cmpWriter.WriteComparison(v, mafAnn, vepAnns); err != nil {
					t.Fatalf("write comparison: %v", err)
				}
			}
			elapsed := time.Since(start)

			counts := cmpWriter.Counts()
			total := cmpWriter.Total()

			// Extract consequence counts
			conseqMatch := counts["consequence"][output.CatMatch] + counts["consequence"][output.CatUpstreamReclass]
			conseqMismatch := counts["consequence"][output.CatMismatch]

			t.Logf("%s: %d variants, %.1f%% consequence match, %s",
				name, total, float64(conseqMatch)/float64(total)*100, elapsed)

			results = append(results, studyResult{
				name:           name,
				variants:       total,
				categoryCounts: counts,
				duration:       elapsed,
				conseqMatch:    conseqMatch,
				conseqMismatch: conseqMismatch,
			})
		})
	}

	// Write markdown report
	reportPath := filepath.Join(tcgaDir, "validation_report.md")
	writeReport(t, reportPath, results, c.TranscriptCount())
}

// studyName extracts the study name from a MAF file path.
func studyName(path string) string {
	base := filepath.Base(path)
	return strings.TrimSuffix(base, "_data_mutations.txt")
}

type studyResult struct {
	name           string
	variants       int
	categoryCounts map[string]map[output.Category]int
	duration       time.Duration
	conseqMatch    int
	conseqMismatch int
}

// writeReport generates a markdown validation report.
func writeReport(t *testing.T, path string, results []studyResult, transcriptCount int) {
	t.Helper()

	var sb strings.Builder
	sb.WriteString("# TCGA Validation Report\n\n")
	sb.WriteString(fmt.Sprintf("Generated: %s  \n", time.Now().UTC().Format("2006-01-02 15:04 UTC")))
	sb.WriteString(fmt.Sprintf("GENCODE transcripts loaded: %d\n\n", transcriptCount))

	// Consequence match table
	sb.WriteString("## Consequence Match\n\n")
	sb.WriteString("| Study | Variants | Match | Mismatch | Match Rate |\n")
	sb.WriteString("|-------|----------|-------|----------|------------|\n")

	var totVariants, totMatch, totMismatch int
	for _, r := range results {
		rate := float64(r.conseqMatch) / float64(r.variants) * 100
		sb.WriteString(fmt.Sprintf("| %s | %d | %d | %d | %.1f%% |\n",
			r.name, r.variants, r.conseqMatch, r.conseqMismatch, rate))
		totVariants += r.variants
		totMatch += r.conseqMatch
		totMismatch += r.conseqMismatch
	}
	totRate := float64(totMatch) / float64(totVariants) * 100
	sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%d** | **%d** | **%.1f%%** |\n\n",
		totVariants, totMatch, totMismatch, totRate))

	// Per-column category breakdown
	colOrder := []string{"consequence", "hgvsp", "hgvsc"}
	colNames := map[string]string{
		"consequence": "Consequence",
		"hgvsp":       "HGVSp",
		"hgvsc":       "HGVSc",
	}

	for _, col := range colOrder {
		sb.WriteString(fmt.Sprintf("## %s Category Breakdown\n\n", colNames[col]))

		// Collect all categories across studies
		allCats := make(map[output.Category]bool)
		for _, r := range results {
			for cat := range r.categoryCounts[col] {
				allCats[cat] = true
			}
		}
		var catList []output.Category
		for cat := range allCats {
			catList = append(catList, cat)
		}
		sort.Slice(catList, func(i, j int) bool { return catList[i] < catList[j] })

		// Header
		sb.WriteString("| Study |")
		for _, cat := range catList {
			sb.WriteString(fmt.Sprintf(" %s |", cat))
		}
		sb.WriteString("\n|-------|")
		for range catList {
			sb.WriteString("------|")
		}
		sb.WriteString("\n")

		// Rows
		totals := make(map[output.Category]int)
		for _, r := range results {
			sb.WriteString(fmt.Sprintf("| %s |", r.name))
			for _, cat := range catList {
				count := r.categoryCounts[col][cat]
				totals[cat] += count
				sb.WriteString(fmt.Sprintf(" %d |", count))
			}
			sb.WriteString("\n")
		}

		// Total row
		sb.WriteString("| **Total** |")
		for _, cat := range catList {
			sb.WriteString(fmt.Sprintf(" **%d** |", totals[cat]))
		}
		sb.WriteString("\n\n")
	}

	// Performance table
	sb.WriteString("## Performance\n\n")
	sb.WriteString("| Study | Variants | Time | Variants/sec |\n")
	sb.WriteString("|-------|----------|------|-------------|\n")

	var totDuration time.Duration
	for _, r := range results {
		vps := float64(r.variants) / r.duration.Seconds()
		sb.WriteString(fmt.Sprintf("| %s | %d | %s | %.0f |\n",
			r.name, r.variants, r.duration.Round(time.Millisecond), vps))
		totDuration += r.duration
	}
	totVPS := float64(totVariants) / totDuration.Seconds()
	sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%s** | **%.0f** |\n",
		totVariants, totDuration.Round(time.Millisecond), totVPS))

	if err := os.WriteFile(path, []byte(sb.String()), 0644); err != nil {
		t.Fatalf("write report: %v", err)
	}
	t.Logf("report written to %s", path)
}

// findTCGADir locates the testdata/tcga directory.
func findTCGADir(t *testing.T) string {
	t.Helper()
	for _, p := range []string{
		filepath.Join("testdata", "tcga"),
		filepath.Join("..", "..", "testdata", "tcga"),
	} {
		if _, err := os.Stat(p); err == nil {
			return p
		}
	}
	t.Fatal("testdata/tcga directory not found")
	return ""
}

// findGENCODEFiles locates GENCODE cache files for GRCh38.
func findGENCODEFiles(t *testing.T) (gtfPath, fastaPath, canonicalPath string) {
	t.Helper()

	home, err := os.UserHomeDir()
	if err != nil {
		t.Fatalf("get home dir: %v", err)
	}
	dir := filepath.Join(home, ".vibe-vep", "grch38")

	matches, err := filepath.Glob(filepath.Join(dir, "gencode.v*.annotation.gtf.gz"))
	if err != nil || len(matches) == 0 {
		t.Fatalf("GENCODE GTF not found in %s — run: vibe-vep download", dir)
	}
	gtfPath = matches[0]

	matches, err = filepath.Glob(filepath.Join(dir, "gencode.v*.pc_transcripts.fa.gz"))
	if err == nil && len(matches) > 0 {
		fastaPath = matches[0]
	}

	cPath := filepath.Join(dir, cache.CanonicalFileName())
	if _, err := os.Stat(cPath); err == nil {
		canonicalPath = cPath
	}

	return
}
