package output_test

import (
	"encoding/json"
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
	"github.com/inodb/vibe-vep/internal/datasource/oncokb"
	"github.com/inodb/vibe-vep/internal/duckdb"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/output"
)

// TestValidationBenchmark runs comparison against all TCGA MAF files (GRCh38)
// and generates a markdown report at testdata/tcga/validation_report.md.
//
// Skipped with -short (large files, not for CI). Run with:
//
//	go test ./internal/output/ -run TestValidationBenchmark$ -v -count=1
func TestValidationBenchmark(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping validation benchmark in short mode")
	}

	studyDir := findStudyDir(t, "tcga")
	mafFiles, err := filepath.Glob(filepath.Join(studyDir, "*_data_mutations.txt"))
	if err != nil {
		t.Fatalf("glob MAF files: %v", err)
	}
	if len(mafFiles) == 0 {
		t.Skip("no TCGA MAF files found in", studyDir)
	}
	sort.Strings(mafFiles)

	c, cacheDuration, loadSource := loadGENCODECache(t, "GRCh38")
	cgl := loadCancerGeneList(t)
	results := runValidationStudies(t, mafFiles, c, cgl)

	reportPath := filepath.Join(studyDir, "validation_report.md")
	meta := reportMeta{
		Assembly:        "GRCh38",
		TranscriptCount: c.TranscriptCount(),
		CacheDuration:   cacheDuration,
		LoadSource:      loadSource,
	}
	writeReport(t, reportPath, results, meta.TranscriptCount, meta.CacheDuration, meta.LoadSource, cgl, meta.Assembly)
	writeReportJSON(t, findDocsDataPath(t, "grch38_validation.json"), results, meta, cgl)
}

// TestValidationBenchmarkGRCh37 runs comparison against GRCh37 MAF files
// and generates a markdown report at testdata/grch37/validation_report.md.
//
// Skipped with -short or when GRCh37 GENCODE cache is not available. Run with:
//
//	go test ./internal/output/ -run TestValidationBenchmarkGRCh37 -v -count=1 -timeout 30m
func TestValidationBenchmarkGRCh37(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping validation benchmark in short mode")
	}

	studyDir := findStudyDir(t, "grch37")
	mafFiles, err := filepath.Glob(filepath.Join(studyDir, "*_data_mutations.txt"))
	if err != nil {
		t.Fatalf("glob MAF files: %v", err)
	}
	if len(mafFiles) == 0 {
		t.Skipf("no GRCh37 MAF files found in %s — run: scripts/download_grch37.sh", studyDir)
	}
	sort.Strings(mafFiles)

	c, cacheDuration, loadSource := loadGENCODECache(t, "GRCh37")
	cgl := loadCancerGeneList(t)
	results := runValidationStudies(t, mafFiles, c, cgl)

	reportPath := filepath.Join(studyDir, "validation_report.md")
	meta := reportMeta{
		Assembly:        "GRCh37",
		TranscriptCount: c.TranscriptCount(),
		CacheDuration:   cacheDuration,
		LoadSource:      loadSource,
	}
	writeReport(t, reportPath, results, meta.TranscriptCount, meta.CacheDuration, meta.LoadSource, cgl, meta.Assembly)
	writeReportJSON(t, findDocsDataPath(t, "grch37_validation.json"), results, meta, cgl)
}

// TestValidationBenchmarkGDC runs comparison against all 33 GDC TCGA MAF files (GRCh38)
// and generates a markdown report at testdata/datahub_gdc/validation_report.md.
//
// Skipped with -short. Run with:
//
//	go test ./internal/output/ -run TestValidationBenchmarkGDC -v -count=1 -timeout 30m
func TestValidationBenchmarkGDC(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping validation benchmark in short mode")
	}

	studyDir := findStudyDir(t, "datahub_gdc")
	mafFiles, err := filepath.Glob(filepath.Join(studyDir, "*_data_mutations.txt"))
	if err != nil {
		t.Fatalf("glob MAF files: %v", err)
	}
	if len(mafFiles) == 0 {
		t.Skipf("no GDC MAF files found in %s — run: scripts/download_datahub_gdc.sh", studyDir)
	}
	sort.Strings(mafFiles)

	c, cacheDuration, loadSource := loadGENCODECache(t, "GRCh38")
	cgl := loadCancerGeneList(t)
	results := runValidationStudies(t, mafFiles, c, cgl)

	reportPath := filepath.Join(studyDir, "validation_report.md")
	meta := reportMeta{
		Assembly:        "GRCh38",
		TranscriptCount: c.TranscriptCount(),
		CacheDuration:   cacheDuration,
		LoadSource:      loadSource,
	}
	writeReport(t, reportPath, results, meta.TranscriptCount, meta.CacheDuration, meta.LoadSource, cgl, meta.Assembly)
	writeReportJSON(t, filepath.Join(studyDir, "validation_report.json"), results, meta, cgl)
}

// TestValidationBenchmarkAll runs comparison against all datahub MAF files
// (GDC + public GRCh38) and generates a report at testdata/datahub_all/GRCh38/validation_report.md.
//
// Skipped with -short. Run with:
//
//	go test ./internal/output/ -run TestValidationBenchmarkAll -v -count=1 -timeout 120m
func TestValidationBenchmarkAll(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping validation benchmark in short mode")
	}

	studyDir := findStudyDir(t, filepath.Join("datahub_all", "GRCh38"))
	mafFiles, err := filepath.Glob(filepath.Join(studyDir, "*_data_mutations.txt"))
	if err != nil {
		t.Fatalf("glob MAF files: %v", err)
	}
	if len(mafFiles) == 0 {
		t.Skipf("no datahub_all MAF files found in %s — run: scripts/download_datahub_all.sh", studyDir)
	}
	sort.Strings(mafFiles)

	c, cacheDuration, loadSource := loadGENCODECache(t, "GRCh38")
	cgl := loadCancerGeneList(t)
	results := runValidationStudies(t, mafFiles, c, cgl)

	reportPath := filepath.Join(studyDir, "validation_report.md")
	meta := reportMeta{
		Assembly:        "GRCh38",
		TranscriptCount: c.TranscriptCount(),
		CacheDuration:   cacheDuration,
		LoadSource:      loadSource,
	}
	writeReport(t, reportPath, results, meta.TranscriptCount, meta.CacheDuration, meta.LoadSource, cgl, meta.Assembly)
	writeReportJSON(t, filepath.Join(studyDir, "validation_report.json"), results, meta, cgl)
}

// TestValidationBenchmarkAllGRCh37 runs comparison against all GRCh37 public studies
// from datahub and generates a report at testdata/datahub_all/GRCh37/validation_report.md.
//
// Skipped with -short. Run with:
//
//	go test ./internal/output/ -run TestValidationBenchmarkAllGRCh37 -v -count=1 -timeout 120m
func TestValidationBenchmarkAllGRCh37(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping validation benchmark in short mode")
	}

	studyDir := findStudyDir(t, filepath.Join("datahub_all", "GRCh37"))
	mafFiles, err := filepath.Glob(filepath.Join(studyDir, "*_data_mutations.txt"))
	if err != nil {
		t.Fatalf("glob MAF files: %v", err)
	}
	if len(mafFiles) == 0 {
		t.Skipf("no GRCh37 datahub_all MAF files found in %s — run: scripts/download_datahub_all.sh", studyDir)
	}
	sort.Strings(mafFiles)

	c, cacheDuration, loadSource := loadGENCODECache(t, "GRCh37")
	cgl := loadCancerGeneList(t)
	results := runValidationStudies(t, mafFiles, c, cgl)

	reportPath := filepath.Join(studyDir, "validation_report.md")
	meta := reportMeta{
		Assembly:        "GRCh37",
		TranscriptCount: c.TranscriptCount(),
		CacheDuration:   cacheDuration,
		LoadSource:      loadSource,
	}
	writeReport(t, reportPath, results, meta.TranscriptCount, meta.CacheDuration, meta.LoadSource, cgl, meta.Assembly)
	writeReportJSON(t, filepath.Join(studyDir, "validation_report.json"), results, meta, cgl)
}

// mismatchRecord captures a single mismatched variant for analysis.
type mismatchRecord struct {
	Study        string `json:"study"`
	Chrom        string `json:"chrom"`
	Pos          int64  `json:"pos"`
	Ref          string `json:"ref"`
	Alt          string `json:"alt"`
	TranscriptID string `json:"transcript_id"`
	Column       string `json:"column"`
	MAFValue     string `json:"maf_value"`
	VEPValue     string `json:"vep_value"`
	Category     string `json:"category"`
}

// TestMismatchCollection collects all mismatched variants from the largest
// available dataset and writes them to a JSON file for analysis.
//
//	go test ./internal/output/ -run TestMismatchCollection -v -count=1 -timeout 120m
func TestMismatchCollection(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping mismatch collection in short mode")
	}

	// Use best available dataset
	studyDir, datasetName := findBestDataset(t)
	mafFiles, err := filepath.Glob(filepath.Join(studyDir, "*_data_mutations.txt"))
	if err != nil {
		t.Fatalf("glob MAF files: %v", err)
	}
	if len(mafFiles) == 0 {
		t.Skipf("no MAF files found in %s", studyDir)
	}
	sort.Strings(mafFiles)
	t.Logf("using dataset: %s (%d files)", datasetName, len(mafFiles))

	c, _, _ := loadGENCODECache(t, "GRCh38")

	categorizer := &output.Categorizer{}
	comparedColumns := []string{"Consequence", "HGVSp_Short", "HGVSc"}
	var mismatches []mismatchRecord

	for _, mafFile := range mafFiles {
		name := studyName(mafFile)
		t.Run(name, func(t *testing.T) {
			parser, err := maf.NewParser(mafFile)
			if err != nil {
				t.Fatalf("open MAF: %v", err)
			}
			defer parser.Close()

			ann := annotate.NewAnnotator(c)
			fileMismatches := 0

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

				bestAnn := output.SelectBestAnnotation(mafAnn, vepAnns)

				left := map[string]string{
					"Consequence": mafAnn.Consequence,
					"HGVSp_Short": mafAnn.HGVSpShort,
					"HGVSc":       mafAnn.HGVSc,
				}
				right := map[string]string{
					"Consequence": "",
					"HGVSp_Short": "",
					"HGVSc":       "",
				}
				if bestAnn != nil {
					right["Consequence"] = bestAnn.Consequence
					right["HGVSp_Short"] = output.HGVSpToShort(bestAnn.HGVSp)
					right["HGVSc"] = bestAnn.HGVSc
				}

				cats := categorizer.CategorizeRow(comparedColumns, left, right,
					comparedColumns, comparedColumns)

				for _, col := range comparedColumns {
					if cats[col] == output.CatMismatch {
						mismatches = append(mismatches, mismatchRecord{
							Study:        name,
							Chrom:        v.Chrom,
							Pos:          v.Pos,
							Ref:          v.Ref,
							Alt:          v.Alt,
							TranscriptID: mafAnn.TranscriptID,
							Column:       col,
							MAFValue:     left[col],
							VEPValue:     right[col],
							Category:     string(cats[col]),
						})
						fileMismatches++
					}
				}
			}
			t.Logf("%s: %d mismatches collected", name, fileMismatches)
		})
	}

	// Group by consequence type for pattern recognition
	type consequenceGroup struct {
		Count    int              `json:"count"`
		Examples []mismatchRecord `json:"examples"`
	}
	groups := make(map[string]*consequenceGroup)
	for _, m := range mismatches {
		if m.Column != "Consequence" {
			continue
		}
		key := m.MAFValue + " -> " + m.VEPValue
		g, ok := groups[key]
		if !ok {
			g = &consequenceGroup{}
			groups[key] = g
		}
		g.Count++
		if len(g.Examples) < 3 {
			g.Examples = append(g.Examples, m)
		}
	}

	report := struct {
		Dataset    string                        `json:"dataset"`
		Total      int                           `json:"total_mismatches"`
		Mismatches []mismatchRecord              `json:"mismatches"`
		Groups     map[string]*consequenceGroup  `json:"consequence_groups"`
	}{
		Dataset:    datasetName,
		Total:      len(mismatches),
		Mismatches: mismatches,
		Groups:     groups,
	}

	data, err := json.MarshalIndent(report, "", "  ")
	if err != nil {
		t.Fatalf("marshal mismatches: %v", err)
	}
	outPath := filepath.Join(studyDir, "mismatches.json")
	if err := os.WriteFile(outPath, data, 0644); err != nil {
		t.Fatalf("write mismatches: %v", err)
	}
	t.Logf("wrote %d mismatches to %s", len(mismatches), outPath)
	t.Logf("consequence mismatch groups: %d", len(groups))
}

// TestMismatchCollectionGRCh37 collects mismatches from GRCh37 public studies.
//
//	go test ./internal/output/ -run TestMismatchCollectionGRCh37 -v -count=1 -timeout 120m
func TestMismatchCollectionGRCh37(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping mismatch collection in short mode")
	}

	studyDir := findStudyDir(t, filepath.Join("datahub_all", "GRCh37"))
	mafFiles, err := filepath.Glob(filepath.Join(studyDir, "*_data_mutations.txt"))
	if err != nil {
		t.Fatalf("glob MAF files: %v", err)
	}
	if len(mafFiles) == 0 {
		t.Skipf("no GRCh37 MAF files found in %s — run: scripts/download_datahub_all.sh", studyDir)
	}
	sort.Strings(mafFiles)
	t.Logf("using GRCh37 dataset (%d files)", len(mafFiles))

	c, _, _ := loadGENCODECache(t, "GRCh37")

	categorizer := &output.Categorizer{}
	comparedColumns := []string{"Consequence", "HGVSp_Short", "HGVSc"}
	var mismatches []mismatchRecord

	for _, mafFile := range mafFiles {
		name := studyName(mafFile)
		t.Run(name, func(t *testing.T) {
			parser, err := maf.NewParser(mafFile)
			if err != nil {
				t.Fatalf("open MAF: %v", err)
			}
			defer parser.Close()

			ann := annotate.NewAnnotator(c)
			fileMismatches := 0

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

				bestAnn := output.SelectBestAnnotation(mafAnn, vepAnns)

				left := map[string]string{
					"Consequence": mafAnn.Consequence,
					"HGVSp_Short": mafAnn.HGVSpShort,
					"HGVSc":       mafAnn.HGVSc,
				}
				right := map[string]string{
					"Consequence": "",
					"HGVSp_Short": "",
					"HGVSc":       "",
				}
				if bestAnn != nil {
					right["Consequence"] = bestAnn.Consequence
					right["HGVSp_Short"] = output.HGVSpToShort(bestAnn.HGVSp)
					right["HGVSc"] = bestAnn.HGVSc
				}

				cats := categorizer.CategorizeRow(comparedColumns, left, right,
					comparedColumns, comparedColumns)

				for _, col := range comparedColumns {
					if cats[col] == output.CatMismatch {
						mismatches = append(mismatches, mismatchRecord{
							Study:        name,
							Chrom:        v.Chrom,
							Pos:          v.Pos,
							Ref:          v.Ref,
							Alt:          v.Alt,
							TranscriptID: mafAnn.TranscriptID,
							Column:       col,
							MAFValue:     left[col],
							VEPValue:     right[col],
							Category:     string(cats[col]),
						})
						fileMismatches++
					}
				}
			}
			t.Logf("%s: %d mismatches collected", name, fileMismatches)
		})
	}

	// Group by consequence type
	type consequenceGroup struct {
		Count    int              `json:"count"`
		Examples []mismatchRecord `json:"examples"`
	}
	groups := make(map[string]*consequenceGroup)
	for _, m := range mismatches {
		if m.Column != "Consequence" {
			continue
		}
		key := m.MAFValue + " -> " + m.VEPValue
		g, ok := groups[key]
		if !ok {
			g = &consequenceGroup{}
			groups[key] = g
		}
		g.Count++
		if len(g.Examples) < 3 {
			g.Examples = append(g.Examples, m)
		}
	}

	report := struct {
		Dataset    string                       `json:"dataset"`
		Assembly   string                       `json:"assembly"`
		Total      int                          `json:"total_mismatches"`
		Mismatches []mismatchRecord             `json:"mismatches"`
		Groups     map[string]*consequenceGroup `json:"consequence_groups"`
	}{
		Dataset:    "datahub_all_grch37",
		Assembly:   "GRCh37",
		Total:      len(mismatches),
		Mismatches: mismatches,
		Groups:     groups,
	}

	data, err := json.MarshalIndent(report, "", "  ")
	if err != nil {
		t.Fatalf("marshal mismatches: %v", err)
	}
	outPath := filepath.Join(studyDir, "mismatches.json")
	if err := os.WriteFile(outPath, data, 0644); err != nil {
		t.Fatalf("write mismatches: %v", err)
	}
	t.Logf("wrote %d mismatches to %s", len(mismatches), outPath)
	t.Logf("consequence mismatch groups: %d", len(groups))
}

// findBestDataset returns the study directory and name for the largest
// available dataset, preferring datahub_all > datahub_gdc > tcga.
func findBestDataset(t *testing.T) (dir, name string) {
	t.Helper()
	candidates := []struct {
		subdir string
		name   string
	}{
		{filepath.Join("datahub_all", "GRCh38"), "datahub_all"},
		{"datahub_gdc", "datahub_gdc"},
		{"tcga", "datahub_gdc_tcga"},
	}
	for _, c := range candidates {
		for _, p := range []string{
			filepath.Join("testdata", c.subdir),
			filepath.Join("..", "..", "testdata", c.subdir),
		} {
			if files, _ := filepath.Glob(filepath.Join(p, "*_data_mutations.txt")); len(files) > 0 {
				return p, c.name
			}
		}
	}
	t.Fatal("no MAF datasets found — run one of the download scripts")
	return "", ""
}

// loadGENCODECache loads GENCODE transcripts for the given assembly,
// preferring gob cache over raw GTF/FASTA. Skips the test if files are not found.
func loadGENCODECache(t *testing.T, assembly string) (c *cache.Cache, duration time.Duration, source string) {
	t.Helper()

	gtfPath, fastaPath, canonicalPath := findGENCODEFiles(t, assembly)
	cacheDir := filepath.Dir(gtfPath)
	c = cache.New()

	cacheStart := time.Now()

	tc := duckdb.NewTranscriptCache(cacheDir)
	gtfFP, err1 := duckdb.StatFile(gtfPath)
	fastaFP, err2 := duckdb.StatFile(fastaPath)
	canonicalFP := duckdb.FileFingerprint{}
	if canonicalPath != "" {
		canonicalFP, _ = duckdb.StatFile(canonicalPath)
	}

	if err1 == nil && err2 == nil && tc.Valid(gtfFP, fastaFP, canonicalFP) {
		if err := tc.Load(c); err != nil {
			t.Fatalf("load transcript cache: %v", err)
		}
		source = "gob cache"
	} else {
		loader := cache.NewGENCODELoader(gtfPath, fastaPath)
		if canonicalPath != "" {
			mskOverrides, ensOverrides, _, err := cache.LoadBiomartCanonicals(canonicalPath)
			if err != nil {
				t.Logf("warning: could not load biomart canonicals: %v", err)
			} else {
				loader.SetCanonicalOverrides(mskOverrides, ensOverrides)
			}
		}
		if err := loader.Load(c); err != nil {
			t.Fatalf("load GENCODE: %v", err)
		}
		source = "GTF/FASTA"
	}

	c.BuildIndex()
	duration = time.Since(cacheStart)
	t.Logf("loaded %d transcripts from %s in %s", c.TranscriptCount(), source, duration.Round(time.Millisecond))
	return c, duration, source
}

// runValidationStudies iterates MAF files, annotates variants, categorizes
// results, and runs parallel benchmarks. Returns per-study results.
func runValidationStudies(t *testing.T, mafFiles []string, c *cache.Cache, cgl oncokb.CancerGeneList) []studyResult {
	t.Helper()

	categorizer := &output.Categorizer{}
	comparedColumns := []string{"Consequence", "HGVSp_Short", "HGVSc"}
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

			catCounts := make(map[string]map[output.Category]int)
			for _, col := range comparedColumns {
				catCounts[col] = make(map[output.Category]int)
			}
			geneMismatches := make(map[string]map[string]int)
			geneTotal := make(map[string]int)
			total := 0

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
				total++

				bestAnn := output.SelectBestAnnotation(mafAnn, vepAnns)

				left := map[string]string{
					"Consequence": mafAnn.Consequence,
					"HGVSp_Short": mafAnn.HGVSpShort,
					"HGVSc":       mafAnn.HGVSc,
				}
				right := map[string]string{
					"Consequence": "",
					"HGVSp_Short": "",
					"HGVSc":       "",
				}
				if bestAnn != nil {
					right["Consequence"] = bestAnn.Consequence
					right["HGVSp_Short"] = output.HGVSpToShort(bestAnn.HGVSp)
					right["HGVSc"] = bestAnn.HGVSc
				}

				cats := categorizer.CategorizeRow(comparedColumns, left, right,
					comparedColumns, comparedColumns)

				for _, col := range comparedColumns {
					catCounts[col][cats[col]]++
				}

				if cgl != nil && cgl.IsCancerGene(mafAnn.HugoSymbol) {
					gene := mafAnn.HugoSymbol
					geneTotal[gene]++
					for col, cat := range cats {
						if cat == output.CatMismatch {
							if geneMismatches[gene] == nil {
								geneMismatches[gene] = make(map[string]int)
							}
							geneMismatches[gene][col]++
						}
					}
				}
			}
			seqDuration := time.Since(start)

			conseqMatch := catCounts["Consequence"][output.CatMatch] + catCounts["Consequence"][output.CatUpstreamReclass]
			conseqMismatch := catCounts["Consequence"][output.CatMismatch]

			parDuration := benchmarkParallel(t, mafFile, ann, runtime.NumCPU())

			seqVPS := float64(total) / seqDuration.Seconds()
			parVPS := float64(total) / parDuration.Seconds()
			speedup := seqDuration.Seconds() / parDuration.Seconds()

			t.Logf("%s: %d variants, %.1f%% consequence match", name, total, float64(conseqMatch)/float64(total)*100)
			t.Logf("  sequential: %s (%.0f variants/sec)", seqDuration.Round(time.Millisecond), seqVPS)
			t.Logf("  parallel (%d workers): %s (%.0f variants/sec, %.2fx speedup)",
				runtime.NumCPU(), parDuration.Round(time.Millisecond), parVPS, speedup)

			results = append(results, studyResult{
				name:           name,
				variants:       total,
				categoryCounts: catCounts,
				seqDuration:    seqDuration,
				parDuration:    parDuration,
				conseqMatch:    conseqMatch,
				conseqMismatch: conseqMismatch,
				geneMismatches: geneMismatches,
				geneTotal:      geneTotal,
			})
		})
	}

	return results
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
	seqDuration    time.Duration
	parDuration    time.Duration
	conseqMatch    int
	conseqMismatch int
	geneMismatches map[string]map[string]int // gene → column → count
	geneTotal      map[string]int            // gene → total variants
}

// benchmarkParallel runs the parallel annotation pipeline for a MAF file
// and returns the wall-clock duration of the annotation phase only.
func benchmarkParallel(t *testing.T, mafFile string, ann *annotate.Annotator, workers int) time.Duration {
	t.Helper()

	parser, err := maf.NewParser(mafFile)
	if err != nil {
		t.Fatalf("open MAF for parallel benchmark: %v", err)
	}
	defer parser.Close()

	items := make(chan annotate.WorkItem, 2*workers)

	start := time.Now()

	go func() {
		defer close(items)
		seq := 0
		for {
			v, mafAnn, err := parser.NextWithAnnotation()
			if err != nil {
				t.Errorf("parallel parse error: %v", err)
				return
			}
			if v == nil {
				return
			}
			items <- annotate.WorkItem{Seq: seq, Variant: v, Extra: mafAnn}
			seq++
		}
	}()

	results := ann.ParallelAnnotate(items, workers)

	if err := annotate.OrderedCollect(results, func(r annotate.WorkResult) error {
		// Consume results (simulates writer work) but discard output.
		_ = r.Anns
		return r.Err
	}); err != nil {
		t.Fatalf("parallel annotation: %v", err)
	}

	return time.Since(start)
}

// writeReport generates a markdown validation report.
func writeReport(t *testing.T, path string, results []studyResult, transcriptCount int, cacheDuration time.Duration, loadSource string, cgl oncokb.CancerGeneList, assembly string) {
	t.Helper()

	sys := getSystemInfo()

	var sb strings.Builder
	sb.WriteString(fmt.Sprintf("# %s Validation Report\n\n", assembly))
	sb.WriteString(fmt.Sprintf("Generated: %s  \n", time.Now().UTC().Format("2006-01-02 15:04 UTC")))
	sb.WriteString(fmt.Sprintf("Assembly: %s  \n", assembly))
	sb.WriteString(fmt.Sprintf("GENCODE transcripts: %d (loaded from %s in %s)  \n", transcriptCount, loadSource, cacheDuration.Round(time.Millisecond)))
	sb.WriteString(fmt.Sprintf("Workers: %d (GOMAXPROCS)  \n", runtime.NumCPU()))
	sb.WriteString(fmt.Sprintf("Studies: %d\n\n", len(results)))
	sb.WriteString("## System\n\n")
	sb.WriteString(fmt.Sprintf("| | |\n|---|---|\n"))
	sb.WriteString(fmt.Sprintf("| CPU | %s |\n", sys.CPU))
	sb.WriteString(fmt.Sprintf("| Cores | %d |\n", sys.Cores))
	sb.WriteString(fmt.Sprintf("| Memory | %s |\n", sys.Memory))
	sb.WriteString(fmt.Sprintf("| OS | %s |\n", sys.OS))
	sb.WriteString(fmt.Sprintf("| Kernel | %s |\n", sys.Kernel))
	sb.WriteString(fmt.Sprintf("| Go | %s |\n\n", sys.GoVer))

	// Match rates table (consequence + HGVSp + HGVSc)
	sb.WriteString("## Match Rates\n\n")
	sb.WriteString("| Study | Variants | Conseq Match | Conseq Mismatch | Conseq Rate | HGVSp Match | HGVSp Mismatch | HGVSp Rate | HGVSc Match | HGVSc Mismatch | HGVSc Rate |\n")
	sb.WriteString("|-------|----------|-------------|-----------------|-------------|-------------|----------------|------------|-------------|----------------|------------|\n")

	var totVariants, totConseqMatch, totConseqMismatch int
	var totHGVSpMatch, totHGVSpMismatch, totHGVScMatch, totHGVScMismatch int
	for _, r := range results {
		conseqRate := float64(r.conseqMatch) / float64(r.variants) * 100
		hgvspMatch := r.categoryCounts["HGVSp_Short"][output.CatMatch]
		hgvspMismatch := r.categoryCounts["HGVSp_Short"][output.CatMismatch]
		hgvspRate := float64(hgvspMatch) / float64(r.variants) * 100
		hgvscMatch := r.categoryCounts["HGVSc"][output.CatMatch]
		hgvscMismatch := r.categoryCounts["HGVSc"][output.CatMismatch]
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
	colOrder := []string{"Consequence", "HGVSp_Short", "HGVSc"}
	colNames := map[string]string{
		"Consequence": "Consequence",
		"HGVSp_Short": "HGVSp",
		"HGVSc":       "HGVSc",
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

		genesWithMismatches := len(aggGeneMismatches)
		totalCancerGenes := len(aggGeneTotal)
		perfectGenes := totalCancerGenes - genesWithMismatches

		if genesWithMismatches == 0 {
			sb.WriteString(fmt.Sprintf("No mismatches across all %d cancer genes tested.\n\n", totalCancerGenes))
		} else {
			sb.WriteString(fmt.Sprintf("%d/%d cancer genes have 100%% match across all columns. Mismatches in %d gene(s):\n\n",
				perfectGenes, totalCancerGenes, genesWithMismatches))
			sb.WriteString("| Gene | Variants | Conseq Mismatches | HGVSp Mismatches | HGVSc Mismatches |\n")
			sb.WriteString("|------|----------|-------------------|------------------|------------------|\n")

			// Sort genes by total mismatches descending
			type geneEntry struct {
				gene  string
				total int
			}
			var genes []geneEntry
			for gene, counts := range aggGeneMismatches {
				total := counts["Consequence"] + counts["HGVSp_Short"] + counts["HGVSc"]
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
					counts["Consequence"], counts["HGVSp_Short"], counts["HGVSc"]))
			}
			sb.WriteString("\n")
		}
	}

	// Performance table
	sb.WriteString("## Performance\n\n")
	sb.WriteString(fmt.Sprintf("Transcript load: %s from %s\n\n", cacheDuration.Round(time.Millisecond), loadSource))
	sb.WriteString("| Study | Variants | Sequential | Seq v/s | Parallel | Par v/s | Speedup |\n")
	sb.WriteString("|-------|----------|-----------|---------|----------|---------|--------|\n")

	var totSeqDuration, totParDuration time.Duration
	for _, r := range results {
		seqVPS := float64(r.variants) / r.seqDuration.Seconds()
		parVPS := float64(r.variants) / r.parDuration.Seconds()
		speedup := r.seqDuration.Seconds() / r.parDuration.Seconds()
		sb.WriteString(fmt.Sprintf("| %s | %d | %s | %.0f | %s | %.0f | %.2fx |\n",
			r.name, r.variants,
			r.seqDuration.Round(time.Millisecond), seqVPS,
			r.parDuration.Round(time.Millisecond), parVPS, speedup))
		totSeqDuration += r.seqDuration
		totParDuration += r.parDuration
	}
	totSeqVPS := float64(totVariants) / totSeqDuration.Seconds()
	totParVPS := float64(totVariants) / totParDuration.Seconds()
	totSpeedup := totSeqDuration.Seconds() / totParDuration.Seconds()
	sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%s** | **%.0f** | **%s** | **%.0f** | **%.2fx** |\n",
		totVariants,
		totSeqDuration.Round(time.Millisecond), totSeqVPS,
		totParDuration.Round(time.Millisecond), totParVPS, totSpeedup))

	if err := os.WriteFile(path, []byte(sb.String()), 0644); err != nil {
		t.Fatalf("write report: %v", err)
	}
	t.Logf("report written to %s", path)
}

// reportMeta holds metadata for JSON report generation.
type reportMeta struct {
	Assembly        string
	TranscriptCount int
	CacheDuration   time.Duration
	LoadSource      string
}

// systemInfo collects hardware and OS information for benchmark reports.
type systemInfo struct {
	CPU    string `json:"cpu"`
	Cores  int    `json:"cores"`
	Memory string `json:"memory"`
	OS     string `json:"os"`
	Kernel string `json:"kernel"`
	GoVer  string `json:"go_version"`
}

func getSystemInfo() systemInfo {
	info := systemInfo{
		Cores: runtime.NumCPU(),
		GoVer: runtime.Version(),
	}

	// CPU model
	if data, err := os.ReadFile("/proc/cpuinfo"); err == nil {
		for _, line := range strings.Split(string(data), "\n") {
			if strings.HasPrefix(line, "model name") {
				parts := strings.SplitN(line, ":", 2)
				if len(parts) == 2 {
					info.CPU = strings.TrimSpace(parts[1])
					break
				}
			}
		}
	}
	if info.CPU == "" {
		info.CPU = runtime.GOARCH
	}

	// Memory
	if data, err := os.ReadFile("/proc/meminfo"); err == nil {
		for _, line := range strings.Split(string(data), "\n") {
			if strings.HasPrefix(line, "MemTotal:") {
				parts := strings.Fields(line)
				if len(parts) >= 2 {
					var kb int64
					if _, err := fmt.Sscanf(parts[1], "%d", &kb); err == nil {
						gb := float64(kb) / 1024 / 1024
						info.Memory = fmt.Sprintf("%.0f GB", gb)
					} else {
						info.Memory = parts[1] + " kB"
					}
				}
				break
			}
		}
	}

	// OS
	if data, err := os.ReadFile("/etc/os-release"); err == nil {
		for _, line := range strings.Split(string(data), "\n") {
			if strings.HasPrefix(line, "PRETTY_NAME=") {
				info.OS = strings.Trim(strings.TrimPrefix(line, "PRETTY_NAME="), "\"")
				break
			}
		}
	}

	// Kernel
	info.Kernel = runtime.GOOS + "/" + runtime.GOARCH
	if data, err := os.ReadFile("/proc/version"); err == nil {
		fields := strings.Fields(string(data))
		if len(fields) >= 3 {
			info.Kernel = fields[2]
		}
	}

	return info
}

// findDocsDataPath locates the docs/data directory relative to the test file
// and returns the full path for the given filename.
func findDocsDataPath(t *testing.T, filename string) string {
	t.Helper()
	for _, rel := range []string{
		filepath.Join("docs", "data"),
		filepath.Join("..", "..", "docs", "data"),
	} {
		if info, err := os.Stat(rel); err == nil && info.IsDir() {
			return filepath.Join(rel, filename)
		}
	}
	t.Fatalf("docs/data directory not found")
	return ""
}

// reportTable is a generic table structure for JSON output.
type reportTable struct {
	Columns []string        `json:"columns"`
	Align   []string        `json:"align"`
	Rows    [][]interface{} `json:"rows"`
	Totals  []interface{}   `json:"totals,omitempty"`
}

// reportJSON is the top-level JSON structure for validation reports.
type reportJSON struct {
	Assembly           string                  `json:"assembly"`
	Generated          string                  `json:"generated"`
	Transcripts        int                     `json:"transcripts"`
	LoadSource         string                  `json:"load_source"`
	LoadDuration       string                  `json:"load_duration"`
	Workers            int                     `json:"workers"`
	Studies            int                     `json:"studies"`
	System             *systemInfo             `json:"system,omitempty"`
	MatchRates         reportTable             `json:"match_rates"`
	CategoryBreakdowns map[string]reportTable   `json:"category_breakdowns"`
	CancerGenes        *reportCancerGenes      `json:"cancer_genes,omitempty"`
	Performance        reportTable             `json:"performance"`
}

type reportCancerGenes struct {
	TotalGenes          int             `json:"total_genes"`
	GenesWithMismatches int             `json:"genes_with_mismatches"`
	Columns             []string        `json:"columns,omitempty"`
	Rows                [][]interface{} `json:"rows,omitempty"`
}

// writeReportJSON generates a JSON validation report for Hugo docs.
func writeReportJSON(t *testing.T, path string, results []studyResult, meta reportMeta, cgl oncokb.CancerGeneList) {
	t.Helper()

	sys := getSystemInfo()
	report := reportJSON{
		Assembly:     meta.Assembly,
		Generated:    time.Now().UTC().Format("2006-01-02 15:04 UTC"),
		Transcripts:  meta.TranscriptCount,
		LoadSource:   meta.LoadSource,
		LoadDuration: meta.CacheDuration.Round(time.Millisecond).String(),
		Workers:      runtime.NumCPU(),
		Studies:      len(results),
		System:       &sys,
	}

	// Match rates table
	report.MatchRates = buildMatchRatesTable(results)

	// Category breakdowns
	report.CategoryBreakdowns = buildCategoryBreakdowns(results)

	// Cancer genes
	if cgl != nil {
		report.CancerGenes = buildCancerGenes(results)
	}

	// Performance
	report.Performance = buildPerformanceTable(results)

	data, err := json.MarshalIndent(report, "", "  ")
	if err != nil {
		t.Fatalf("marshal JSON report: %v", err)
	}
	if err := os.WriteFile(path, data, 0644); err != nil {
		t.Fatalf("write JSON report: %v", err)
	}
	t.Logf("JSON report written to %s", path)
}

func buildMatchRatesTable(results []studyResult) reportTable {
	columns := []string{"Study", "Variants", "Conseq Match", "Conseq Mismatch", "Conseq Rate", "HGVSp Match", "HGVSp Mismatch", "HGVSp Rate", "HGVSc Match", "HGVSc Mismatch", "HGVSc Rate"}
	align := []string{"left", "right", "right", "right", "right", "right", "right", "right", "right", "right", "right"}

	var rows [][]interface{}
	var totV, totCM, totCX, totPM, totPX, totHM, totHX int
	for _, r := range results {
		hgvspMatch := r.categoryCounts["HGVSp_Short"][output.CatMatch]
		hgvspMismatch := r.categoryCounts["HGVSp_Short"][output.CatMismatch]
		hgvscMatch := r.categoryCounts["HGVSc"][output.CatMatch]
		hgvscMismatch := r.categoryCounts["HGVSc"][output.CatMismatch]
		conseqRate := fmt.Sprintf("%.1f%%", float64(r.conseqMatch)/float64(r.variants)*100)
		hgvspRate := fmt.Sprintf("%.1f%%", float64(hgvspMatch)/float64(r.variants)*100)
		hgvscRate := fmt.Sprintf("%.1f%%", float64(hgvscMatch)/float64(r.variants)*100)

		rows = append(rows, []interface{}{
			r.name, r.variants,
			r.conseqMatch, r.conseqMismatch, conseqRate,
			hgvspMatch, hgvspMismatch, hgvspRate,
			hgvscMatch, hgvscMismatch, hgvscRate,
		})
		totV += r.variants
		totCM += r.conseqMatch
		totCX += r.conseqMismatch
		totPM += hgvspMatch
		totPX += hgvspMismatch
		totHM += hgvscMatch
		totHX += hgvscMismatch
	}
	totals := []interface{}{
		"Total", totV,
		totCM, totCX, fmt.Sprintf("%.1f%%", float64(totCM)/float64(totV)*100),
		totPM, totPX, fmt.Sprintf("%.1f%%", float64(totPM)/float64(totV)*100),
		totHM, totHX, fmt.Sprintf("%.1f%%", float64(totHM)/float64(totV)*100),
	}

	return reportTable{Columns: columns, Align: align, Rows: rows, Totals: totals}
}

func buildCategoryBreakdowns(results []studyResult) map[string]reportTable {
	breakdowns := make(map[string]reportTable)
	colOrder := []string{"Consequence", "HGVSp_Short", "HGVSc"}
	colNames := map[string]string{
		"Consequence": "Consequence",
		"HGVSp_Short": "HGVSp",
		"HGVSc":       "HGVSc",
	}

	for _, col := range colOrder {
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

		columns := []string{"Study"}
		align := []string{"left"}
		for _, cat := range catList {
			columns = append(columns, string(cat))
			align = append(align, "right")
		}

		var rows [][]interface{}
		totals := make([]int, len(catList))
		for _, r := range results {
			row := []interface{}{r.name}
			for i, cat := range catList {
				count := r.categoryCounts[col][cat]
				totals[i] += count
				row = append(row, count)
			}
			rows = append(rows, row)
		}

		totalRow := []interface{}{"Total"}
		for _, t := range totals {
			totalRow = append(totalRow, t)
		}

		breakdowns[colNames[col]] = reportTable{
			Columns: columns,
			Align:   align,
			Rows:    rows,
			Totals:  totalRow,
		}
	}
	return breakdowns
}

func buildCancerGenes(results []studyResult) *reportCancerGenes {
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

	cg := &reportCancerGenes{
		TotalGenes:          len(aggGeneTotal),
		GenesWithMismatches: len(aggGeneMismatches),
	}

	if len(aggGeneMismatches) > 0 {
		cg.Columns = []string{"Gene", "Variants", "Conseq Mismatches", "HGVSp Mismatches", "HGVSc Mismatches"}

		type geneEntry struct {
			gene  string
			total int
		}
		var genes []geneEntry
		for gene, counts := range aggGeneMismatches {
			total := counts["Consequence"] + counts["HGVSp_Short"] + counts["HGVSc"]
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
			cg.Rows = append(cg.Rows, []interface{}{
				ge.gene, aggGeneTotal[ge.gene],
				counts["Consequence"], counts["HGVSp_Short"], counts["HGVSc"],
			})
		}
	}
	return cg
}

func buildPerformanceTable(results []studyResult) reportTable {
	columns := []string{"Study", "Variants", "Sequential", "Seq v/s", "Parallel", "Par v/s", "Speedup"}
	align := []string{"left", "right", "right", "right", "right", "right", "right"}

	var rows [][]interface{}
	var totV int
	var totSeq, totPar time.Duration
	for _, r := range results {
		seqVPS := float64(r.variants) / r.seqDuration.Seconds()
		parVPS := float64(r.variants) / r.parDuration.Seconds()
		speedup := r.seqDuration.Seconds() / r.parDuration.Seconds()
		rows = append(rows, []interface{}{
			r.name, r.variants,
			r.seqDuration.Round(time.Millisecond).String(), fmt.Sprintf("%.0f", seqVPS),
			r.parDuration.Round(time.Millisecond).String(), fmt.Sprintf("%.0f", parVPS),
			fmt.Sprintf("%.2fx", speedup),
		})
		totV += r.variants
		totSeq += r.seqDuration
		totPar += r.parDuration
	}
	totSeqVPS := float64(totV) / totSeq.Seconds()
	totParVPS := float64(totV) / totPar.Seconds()
	totSpeedup := totSeq.Seconds() / totPar.Seconds()
	totals := []interface{}{
		"Total", totV,
		totSeq.Round(time.Millisecond).String(), fmt.Sprintf("%.0f", totSeqVPS),
		totPar.Round(time.Millisecond).String(), fmt.Sprintf("%.0f", totParVPS),
		fmt.Sprintf("%.2fx", totSpeedup),
	}

	return reportTable{Columns: columns, Align: align, Rows: rows, Totals: totals}
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

// findStudyDir locates a testdata subdirectory (e.g. "tcga", "grch37").
func findStudyDir(t *testing.T, subdir string) string {
	t.Helper()
	for _, p := range []string{
		filepath.Join("testdata", subdir),
		filepath.Join("..", "..", "testdata", subdir),
	} {
		if _, err := os.Stat(p); err == nil {
			return p
		}
	}
	t.Fatalf("testdata/%s directory not found", subdir)
	return ""
}

// findGENCODEFiles locates GENCODE cache files for the given assembly.
// Finds GENCODE files in the default cache directory for the given assembly.
func findGENCODEFiles(t *testing.T, assembly string) (gtfPath, fastaPath, canonicalPath string) {
	t.Helper()

	home, err := os.UserHomeDir()
	if err != nil {
		t.Fatalf("get home dir: %v", err)
	}
	dir := filepath.Join(home, ".vibe-vep", strings.ToLower(assembly))

	// Both assemblies use gencode.v*.annotation.gtf.gz (GRCh37 uses native v19, not liftover).
	gtfPattern := "gencode.v*.annotation.gtf.gz"
	fastaPattern := "gencode.v*.pc_transcripts.fa.gz"

	matches, err := filepath.Glob(filepath.Join(dir, gtfPattern))
	if err != nil || len(matches) == 0 {
		t.Skipf("GENCODE GTF not found in %s — run: vibe-vep download --assembly %s", dir, assembly)
	}
	gtfPath = matches[0]

	matches, err = filepath.Glob(filepath.Join(dir, fastaPattern))
	if err == nil && len(matches) > 0 {
		fastaPath = matches[0]
	}

	cPath := filepath.Join(dir, cache.CanonicalFileName())
	if _, err := os.Stat(cPath); err == nil {
		canonicalPath = cPath
	}

	return
}
