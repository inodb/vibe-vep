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
	"github.com/inodb/vibe-vep/internal/datasource/oncokb"
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

	// Load cancer gene list
	cgl := loadCancerGeneList(t)

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

			// Track per-gene mismatches for cancer genes
			geneMismatches := make(map[string]map[string]int) // gene → column → mismatch count
			geneTotal := make(map[string]int)                 // gene → total count

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

				// Track cancer gene mismatches
				if cgl != nil && cgl.IsCancerGene(mafAnn.HugoSymbol) {
					gene := mafAnn.HugoSymbol
					geneTotal[gene]++
					counts := cmpWriter.LastCategories()
					for col, cat := range counts {
						if cat == output.CatMismatch {
							if geneMismatches[gene] == nil {
								geneMismatches[gene] = make(map[string]int)
							}
							geneMismatches[gene][col]++
						}
					}
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
				geneMismatches: geneMismatches,
				geneTotal:      geneTotal,
			})
		})
	}

	// Write markdown report
	reportPath := filepath.Join(tcgaDir, "validation_report.md")
	writeReport(t, reportPath, results, c.TranscriptCount(), cgl)
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
	geneMismatches map[string]map[string]int // gene → column → count
	geneTotal      map[string]int            // gene → total variants
}

// writeReport generates a markdown validation report.
func writeReport(t *testing.T, path string, results []studyResult, transcriptCount int, cgl oncokb.CancerGeneList) {
	t.Helper()

	var sb strings.Builder
	sb.WriteString("# TCGA Validation Report\n\n")
	sb.WriteString(fmt.Sprintf("Generated: %s  \n", time.Now().UTC().Format("2006-01-02 15:04 UTC")))
	sb.WriteString(fmt.Sprintf("GENCODE transcripts loaded: %d\n\n", transcriptCount))

	// Match rates table (consequence + HGVSp + HGVSc)
	sb.WriteString("## Match Rates\n\n")
	sb.WriteString("| Study | Variants | Conseq Match | Conseq Mismatch | Conseq Rate | HGVSp Match | HGVSp Mismatch | HGVSp Rate | HGVSc Match | HGVSc Mismatch | HGVSc Rate |\n")
	sb.WriteString("|-------|----------|-------------|-----------------|-------------|-------------|----------------|------------|-------------|----------------|------------|\n")

	var totVariants, totConseqMatch, totConseqMismatch int
	var totHGVSpMatch, totHGVSpMismatch, totHGVScMatch, totHGVScMismatch int
	for _, r := range results {
		conseqRate := float64(r.conseqMatch) / float64(r.variants) * 100
		hgvspMatch := r.categoryCounts["hgvsp"][output.CatMatch]
		hgvspMismatch := r.categoryCounts["hgvsp"][output.CatMismatch]
		hgvspRate := float64(hgvspMatch) / float64(r.variants) * 100
		hgvscMatch := r.categoryCounts["hgvsc"][output.CatMatch]
		hgvscMismatch := r.categoryCounts["hgvsc"][output.CatMismatch]
		hgvscRate := float64(hgvscMatch) / float64(r.variants) * 100
		sb.WriteString(fmt.Sprintf("| %s | %d | %d | %d | %.1f%% | %d | %d | %.1f%% | %d | %d | %.1f%% |\n",
			r.name, r.variants,
			r.conseqMatch, r.conseqMismatch, conseqRate,
			hgvspMatch, hgvspMismatch, hgvspRate,
			hgvscMatch, hgvscMismatch, hgvscRate))
		totVariants += r.variants
		totConseqMatch += r.conseqMatch
		totConseqMismatch += r.conseqMismatch
		totHGVSpMatch += hgvspMatch
		totHGVSpMismatch += hgvspMismatch
		totHGVScMatch += hgvscMatch
		totHGVScMismatch += hgvscMismatch
	}
	totConseqRate := float64(totConseqMatch) / float64(totVariants) * 100
	totHGVSpRate := float64(totHGVSpMatch) / float64(totVariants) * 100
	totHGVScRate := float64(totHGVScMatch) / float64(totVariants) * 100
	sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%d** | **%d** | **%.1f%%** | **%d** | **%d** | **%.1f%%** | **%d** | **%d** | **%.1f%%** |\n\n",
		totVariants,
		totConseqMatch, totConseqMismatch, totConseqRate,
		totHGVSpMatch, totHGVSpMismatch, totHGVSpRate,
		totHGVScMatch, totHGVScMismatch, totHGVScRate))

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

	// Cancer gene mismatches
	if cgl != nil {
		sb.WriteString("## Cancer Gene Mismatches\n\n")

		// Aggregate across studies
		aggGeneMismatches := make(map[string]map[string]int)
		aggGeneTotal := make(map[string]int)
		for _, r := range results {
			for gene, counts := range r.geneMismatches {
				if aggGeneMismatches[gene] == nil {
					aggGeneMismatches[gene] = make(map[string]int)
				}
				for col, n := range counts {
					aggGeneMismatches[gene][col] += n
				}
			}
			for gene, n := range r.geneTotal {
				aggGeneTotal[gene] += n
			}
		}

		if len(aggGeneMismatches) == 0 {
			sb.WriteString("No mismatches in cancer genes.\n\n")
		} else {
			sb.WriteString("| Gene | Variants | Conseq Mismatches | HGVSp Mismatches | HGVSc Mismatches |\n")
			sb.WriteString("|------|----------|-------------------|------------------|------------------|\n")

			// Sort genes by total mismatches descending
			type geneEntry struct {
				gene  string
				total int
			}
			var genes []geneEntry
			for gene, counts := range aggGeneMismatches {
				total := counts["consequence"] + counts["hgvsp"] + counts["hgvsc"]
				genes = append(genes, geneEntry{gene, total})
			}
			sort.Slice(genes, func(i, j int) bool {
				if genes[i].total != genes[j].total {
					return genes[i].total > genes[j].total
				}
				return genes[i].gene < genes[j].gene
			})

			for _, ge := range genes {
				counts := aggGeneMismatches[ge.gene]
				sb.WriteString(fmt.Sprintf("| %s | %d | %d | %d | %d |\n",
					ge.gene, aggGeneTotal[ge.gene],
					counts["consequence"], counts["hgvsp"], counts["hgvsc"]))
			}
			sb.WriteString("\n")
		}
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

// loadCancerGeneList attempts to load the OncoKB cancer gene list from the repo root.
func loadCancerGeneList(t *testing.T) oncokb.CancerGeneList {
	t.Helper()
	for _, rel := range []string{
		filepath.Join("cancerGeneList.tsv"),
		filepath.Join("..", "..", "cancerGeneList.tsv"),
	} {
		if _, err := os.Stat(rel); err == nil {
			cgl, err := oncokb.LoadCancerGeneList(rel)
			if err != nil {
				t.Logf("warning: could not load cancer gene list: %v", err)
				return nil
			}
			t.Logf("loaded %d cancer genes", len(cgl))
			return cgl
		}
	}
	t.Log("cancerGeneList.tsv not found, skipping cancer gene tracking")
	return nil
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
