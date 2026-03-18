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
	"github.com/inodb/vibe-vep/internal/datasource/hotspots"
	"github.com/inodb/vibe-vep/internal/genomicindex"
	"github.com/inodb/vibe-vep/internal/maf"
)

// sourcesCtx holds all loaded annotation sources for benchmarking.
type sourcesCtx struct {
	giStore *genomicindex.Store // unified genomic index
	hsStore *hotspots.Store
}

// TestAnnotationSourcesBenchmark runs annotation source coverage and performance benchmarks
// against datahub GDC MAF files and generates testdata/datahub_gdc/annotation_sources_report.md.
//
// Skipped with -short. Run with:
//
//	CGO_ENABLED=1 go test ./internal/output/ -run TestAnnotationSourcesBenchmark -v -count=1 -timeout 30m
func TestAnnotationSourcesBenchmark(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping annotation sources benchmark in short mode")
	}

	// Find datahub GDC MAF files (preferred), fall back to TCGA.
	studyDir := ""
	for _, dir := range []string{"datahub_gdc", "tcga"} {
		d := findStudyDir(t, dir)
		mafFiles, err := filepath.Glob(filepath.Join(d, "*_data_mutations.txt"))
		if err == nil && len(mafFiles) > 0 {
			studyDir = d
			break
		}
	}
	if studyDir == "" {
		t.Skip("no MAF files found in datahub_gdc or tcga")
	}
	mafFiles, _ := filepath.Glob(filepath.Join(studyDir, "*_data_mutations.txt"))
	sort.Strings(mafFiles)
	t.Logf("found %d MAF files in %s", len(mafFiles), studyDir)

	// Load GENCODE cache.
	gtfPath, fastaPath, canonicalPath := findGENCODEFiles(t, "GRCh38")
	c := cache.New()
	loader := cache.NewGENCODELoader(gtfPath, fastaPath)
	if canonicalPath != "" {
		mskOverrides, ensOverrides, err := cache.LoadBiomartCanonicals(canonicalPath)
		if err != nil {
			t.Logf("warning: could not load biomart canonicals: %v", err)
		} else {
			loader.SetCanonicalOverrides(mskOverrides, ensOverrides)
		}
	}
	if err := loader.Load(c); err != nil {
		t.Fatalf("load GENCODE cache: %v", err)
	}
	t.Logf("loaded %d transcripts", c.TranscriptCount())

	var ctx sourcesCtx

	// Open unified genomic index.
	home, err := os.UserHomeDir()
	if err != nil {
		t.Fatalf("get home dir: %v", err)
	}
	dbPath := filepath.Join(home, ".vibe-vep", "grch38", "genomic_annotations.sqlite")
	if _, err := os.Stat(dbPath); os.IsNotExist(err) {
		t.Skipf("genomic index not found at %s — run: vibe-vep prepare", dbPath)
	}
	ctx.giStore, err = genomicindex.Open(dbPath)
	if err != nil {
		t.Fatalf("open genomic index: %v", err)
	}
	defer ctx.giStore.Close()
	t.Logf("opened genomic index at %s", dbPath)

	// Load Hotspots.
	hotspotsPath := "/home/ino/vibe-vep/genome-nexus-importer/data/common_input/hotspots_v2_and_3d.txt"
	if _, err := os.Stat(hotspotsPath); err == nil {
		ctx.hsStore, err = hotspots.Load(hotspotsPath)
		if err != nil {
			t.Fatalf("load hotspots: %v", err)
		}
		t.Logf("loaded hotspots: %d transcripts, %d positions", ctx.hsStore.TranscriptCount(), ctx.hsStore.HotspotCount())
	} else {
		t.Logf("hotspots file not found, skipping")
	}

	var results []sourceStudyResult

	for _, mafFile := range mafFiles {
		name := studyName(mafFile)
		t.Run(name, func(t *testing.T) {
			r := benchmarkStudy(t, mafFile, c, &ctx)
			r.name = name
			results = append(results, r)
		})
	}

	// Write report.
	reportPath := filepath.Join(studyDir, "annotation_sources_report.md")
	writeAnnotationSourcesReport(t, reportPath, results, c.TranscriptCount(), &ctx)
}

// lookupKey identifies a unique variant for deduplication.
type lookupKey struct {
	chrom string
	pos   int64
	ref   string
	alt   string
}

// benchmarkStudy annotates a MAF file and benchmarks all annotation source lookups.
func benchmarkStudy(t *testing.T, mafFile string, c *cache.Cache, ctx *sourcesCtx) sourceStudyResult {
	t.Helper()

	parser, err := maf.NewParser(mafFile)
	if err != nil {
		t.Fatalf("open MAF: %v", err)
	}
	defer parser.Close()

	ann := annotate.NewAnnotator(c)

	var (
		totalVariants int
		missenseCount int
		baseTime      time.Duration
		giLookupTime  time.Duration
		amHits        int
		classCounts   = make(map[string]int)
		hsHits        int
		hsChecked     int
		hsTypeCounts  = make(map[string]int)
		cvHits        int
		cvSigCounts   = make(map[string]int)
		sigHits       int
		gnomadHits    int
		siftHits      int
		siftPredCounts = make(map[string]int)
		pp2Hits       int
		pp2PredCounts = make(map[string]int)
		dbsnpHits     int
	)
	seenMissense := make(map[lookupKey]bool)

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

		// Unified genomic index lookup (all chroms normalized without "chr" prefix).
		if ctx.giStore != nil {
			chrom := v.NormalizeChrom()
			lookupStart := time.Now()
			r, ok := ctx.giStore.Lookup(chrom, v.Pos, v.Ref, v.Alt)
			giLookupTime += time.Since(lookupStart)

			if ok {
				// Track AlphaMissense metrics.
				if r.AMScore > 0 {
					hasMissense := false
					for _, a := range vepAnns {
						if isMissenseConsequence(a.Consequence) {
							hasMissense = true
							break
						}
					}
					if hasMissense {
						k := lookupKey{chrom: chrom, pos: v.Pos, ref: v.Ref, alt: v.Alt}
						if !seenMissense[k] {
							seenMissense[k] = true
							amHits++
							classCounts[r.AMClass]++
						}
					}
				}

				// Track ClinVar metrics.
				if r.CVClnSig != "" {
					cvHits++
					cvSigCounts[r.CVClnSig]++
				}

				// Track SIGNAL metrics.
				if r.SigMutStatus != "" {
					sigHits++
				}

				// Track gnomAD metrics.
				if r.GnomadAF != "" {
					gnomadHits++
				}

				// Track SIFT metrics (missense-relevant).
				if r.SiftScore > 0 || r.SiftPred != "" {
					siftHits++
					if r.SiftPred != "" {
						siftPredCounts[r.SiftPred]++
					}
				}

				// Track PolyPhen-2 metrics (missense-relevant).
				if r.PP2Score > 0 || r.PP2Pred != "" {
					pp2Hits++
					if r.PP2Pred != "" {
						pp2PredCounts[r.PP2Pred]++
					}
				}

				// Track dbSNP metrics.
				if r.DbSnpID != "" {
					dbsnpHits++
				}
			}
		}

		// Count unique missense keys (for coverage calculation).
		for _, a := range vepAnns {
			if isMissenseConsequence(a.Consequence) {
				chrom := v.NormalizeChrom()
				k := lookupKey{chrom: chrom, pos: v.Pos, ref: v.Ref, alt: v.Alt}
				if !seenMissense[k] {
					seenMissense[k] = true
				}
				missenseCount++
				break
			}
		}

		// Hotspot lookup (by transcript + protein position).
		if ctx.hsStore != nil {
			for _, a := range vepAnns {
				if a.ProteinPosition > 0 && a.TranscriptID != "" {
					txID := a.TranscriptID
					if i := strings.IndexByte(txID, '.'); i >= 0 {
						txID = txID[:i]
					}
					hsChecked++
					if h, ok := ctx.hsStore.Lookup(txID, a.ProteinPosition); ok {
						hsHits++
						hsTypeCounts[h.Type]++
					}
					break
				}
			}
		}
	}

	// Deduplicated missense count for AM coverage.
	uniqueMissense := 0
	for range seenMissense {
		uniqueMissense++
	}

	amCoverage := 0.0
	if uniqueMissense > 0 {
		amCoverage = float64(amHits) / float64(uniqueMissense) * 100
	}
	t.Logf("%d variants, %d missense, %d AM hits (%.1f%%), lookup time %s",
		totalVariants, uniqueMissense, amHits, amCoverage, giLookupTime.Round(time.Millisecond))
	if ctx.hsStore != nil {
		t.Logf("  hotspots: %d checked, %d hits", hsChecked, hsHits)
	}
	t.Logf("  ClinVar: %d hits, gnomAD: %d hits, SIGNAL: %d hits",
		cvHits, gnomadHits, sigHits)
	t.Logf("  SIFT: %d hits, PolyPhen-2: %d hits, dbSNP: %d hits",
		siftHits, pp2Hits, dbsnpHits)

	return sourceStudyResult{
		variants:       totalVariants,
		missense:       uniqueMissense,
		amHits:         amHits,
		classCounts:    classCounts,
		baseTime:       baseTime,
		giLookupTime:   giLookupTime,
		hsChecked:      hsChecked,
		hsHits:         hsHits,
		hsTypeCounts:   hsTypeCounts,
		cvHits:         cvHits,
		cvSigCounts:    cvSigCounts,
		sigHits:        sigHits,
		gnomadHits:     gnomadHits,
		siftHits:       siftHits,
		siftPredCounts: siftPredCounts,
		pp2Hits:        pp2Hits,
		pp2PredCounts:  pp2PredCounts,
		dbsnpHits:      dbsnpHits,
	}
}

type sourceStudyResult struct {
	name           string
	variants       int
	missense       int
	amHits         int
	classCounts    map[string]int
	baseTime       time.Duration
	giLookupTime   time.Duration
	hsChecked      int
	hsHits         int
	hsTypeCounts   map[string]int
	cvHits         int
	cvSigCounts    map[string]int
	sigHits        int
	gnomadHits     int
	siftHits       int
	siftPredCounts map[string]int
	pp2Hits        int
	pp2PredCounts  map[string]int
	dbsnpHits      int
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

func writeAnnotationSourcesReport(t *testing.T, path string, results []sourceStudyResult, transcriptCount int, ctx *sourcesCtx) {
	t.Helper()

	var sb strings.Builder
	sb.WriteString("# Annotation Sources Report\n\n")
	sb.WriteString(fmt.Sprintf("Generated: %s  \n", time.Now().UTC().Format("2006-01-02 15:04 UTC")))
	sb.WriteString(fmt.Sprintf("GENCODE transcripts: %d  \n", transcriptCount))
	sb.WriteString("Data: unified genomic index (SQLite)  \n")
	if ctx.hsStore != nil {
		sb.WriteString(fmt.Sprintf("Cancer Hotspots: %d transcripts, %d positions  \n", ctx.hsStore.TranscriptCount(), ctx.hsStore.HotspotCount()))
	}
	sb.WriteString(fmt.Sprintf("Workers: %d (GOMAXPROCS)\n\n", runtime.NumCPU()))

	// Compute totals.
	var totVariants, totMissense, totAMHits int
	var totCVHits, totSigHits, totGnomadHits int
	var totSiftHits, totPP2Hits, totDbsnpHits int
	totClass := make(map[string]int)
	totCVSig := make(map[string]int)
	totSiftPred := make(map[string]int)
	totPP2Pred := make(map[string]int)
	for _, r := range results {
		totVariants += r.variants
		totMissense += r.missense
		totAMHits += r.amHits
		totCVHits += r.cvHits
		totSigHits += r.sigHits
		totGnomadHits += r.gnomadHits
		totSiftHits += r.siftHits
		totPP2Hits += r.pp2Hits
		totDbsnpHits += r.dbsnpHits
		for k, v := range r.classCounts {
			totClass[k] += v
		}
		for k, v := range r.cvSigCounts {
			totCVSig[k] += v
		}
		for k, v := range r.siftPredCounts {
			totSiftPred[k] += v
		}
		for k, v := range r.pp2PredCounts {
			totPP2Pred[k] += v
		}
	}

	// --- AlphaMissense Coverage ---
	sb.WriteString("## AlphaMissense Coverage\n\n")
	sb.WriteString("| Study | Variants | Missense | AM Hits | Coverage | likely_benign | ambiguous | likely_pathogenic |\n")
	sb.WriteString("|-------|----------|----------|---------|----------|---------------|-----------|-------------------|\n")

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
	}
	totCoverage := 0.0
	if totMissense > 0 {
		totCoverage = float64(totAMHits) / float64(totMissense) * 100
	}
	sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%d** | **%d** | **%.1f%%** | **%d** | **%d** | **%d** |\n\n",
		totVariants, totMissense, totAMHits, totCoverage,
		totClass["likely_benign"],
		totClass["ambiguous"],
		totClass["likely_pathogenic"]))

	// --- Cancer Hotspots Coverage ---
	if ctx.hsStore != nil {
		sb.WriteString("## Cancer Hotspots Coverage\n\n")
		sb.WriteString("| Study | Variants | Checked | Hotspot Hits | Hit Rate | single residue | in-frame indel | 3d | splice site |\n")
		sb.WriteString("|-------|----------|---------|--------------|----------|----------------|----------------|----|-------------|\n")

		var totChecked, totHsHits int
		totHsType := make(map[string]int)
		for _, r := range results {
			hitRate := 0.0
			if r.hsChecked > 0 {
				hitRate = float64(r.hsHits) / float64(r.hsChecked) * 100
			}
			sb.WriteString(fmt.Sprintf("| %s | %d | %d | %d | %.2f%% | %d | %d | %d | %d |\n",
				r.name, r.variants, r.hsChecked, r.hsHits, hitRate,
				r.hsTypeCounts["single residue"],
				r.hsTypeCounts["in-frame indel"],
				r.hsTypeCounts["3d"],
				r.hsTypeCounts["splice site"]))
			totChecked += r.hsChecked
			totHsHits += r.hsHits
			for typ, n := range r.hsTypeCounts {
				totHsType[typ] += n
			}
		}
		totHsRate := 0.0
		if totChecked > 0 {
			totHsRate = float64(totHsHits) / float64(totChecked) * 100
		}
		sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%d** | **%d** | **%.2f%%** | **%d** | **%d** | **%d** | **%d** |\n\n",
			totVariants, totChecked, totHsHits, totHsRate,
			totHsType["single residue"],
			totHsType["in-frame indel"],
			totHsType["3d"],
			totHsType["splice site"]))
	}

	// --- ClinVar Coverage ---
	sb.WriteString("## ClinVar Coverage\n\n")
	sb.WriteString("| Study | Variants | ClinVar Hits | Hit Rate | Pathogenic | Likely_pathogenic | Uncertain | Benign | Likely_benign | Other |\n")
	sb.WriteString("|-------|----------|--------------|----------|------------|-------------------|-----------|--------|---------------|-------|\n")

	for _, r := range results {
		hitRate := 0.0
		if r.variants > 0 {
			hitRate = float64(r.cvHits) / float64(r.variants) * 100
		}
		other := r.cvHits
		for _, key := range []string{"Pathogenic", "Likely_pathogenic", "Uncertain_significance", "Benign", "Likely_benign"} {
			other -= r.cvSigCounts[key]
		}
		sb.WriteString(fmt.Sprintf("| %s | %d | %d | %.2f%% | %d | %d | %d | %d | %d | %d |\n",
			r.name, r.variants, r.cvHits, hitRate,
			r.cvSigCounts["Pathogenic"],
			r.cvSigCounts["Likely_pathogenic"],
			r.cvSigCounts["Uncertain_significance"],
			r.cvSigCounts["Benign"],
			r.cvSigCounts["Likely_benign"],
			other))
	}
	totCVRate := 0.0
	if totVariants > 0 {
		totCVRate = float64(totCVHits) / float64(totVariants) * 100
	}
	totOther := totCVHits
	for _, key := range []string{"Pathogenic", "Likely_pathogenic", "Uncertain_significance", "Benign", "Likely_benign"} {
		totOther -= totCVSig[key]
	}
	sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%d** | **%.2f%%** | **%d** | **%d** | **%d** | **%d** | **%d** | **%d** |\n\n",
		totVariants, totCVHits, totCVRate,
		totCVSig["Pathogenic"],
		totCVSig["Likely_pathogenic"],
		totCVSig["Uncertain_significance"],
		totCVSig["Benign"],
		totCVSig["Likely_benign"],
		totOther))

	// --- gnomAD Coverage ---
	sb.WriteString("## gnomAD Coverage\n\n")
	sb.WriteString("| Study | Variants | gnomAD Hits | Hit Rate |\n")
	sb.WriteString("|-------|----------|-------------|----------|\n")

	for _, r := range results {
		hitRate := 0.0
		if r.variants > 0 {
			hitRate = float64(r.gnomadHits) / float64(r.variants) * 100
		}
		sb.WriteString(fmt.Sprintf("| %s | %d | %d | %.2f%% |\n",
			r.name, r.variants, r.gnomadHits, hitRate))
	}
	totGnomadRate := 0.0
	if totVariants > 0 {
		totGnomadRate = float64(totGnomadHits) / float64(totVariants) * 100
	}
	sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%d** | **%.2f%%** |\n\n",
		totVariants, totGnomadHits, totGnomadRate))

	// --- SIFT Coverage ---
	sb.WriteString("## SIFT Coverage\n\n")
	sb.WriteString("| Study | Variants | SIFT Hits | Hit Rate | deleterious | tolerated |\n")
	sb.WriteString("|-------|----------|-----------|----------|-------------|----------|\n")

	for _, r := range results {
		hitRate := 0.0
		if r.variants > 0 {
			hitRate = float64(r.siftHits) / float64(r.variants) * 100
		}
		sb.WriteString(fmt.Sprintf("| %s | %d | %d | %.2f%% | %d | %d |\n",
			r.name, r.variants, r.siftHits, hitRate,
			r.siftPredCounts["deleterious"],
			r.siftPredCounts["tolerated"]))
	}
	totSiftRate := 0.0
	if totVariants > 0 {
		totSiftRate = float64(totSiftHits) / float64(totVariants) * 100
	}
	sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%d** | **%.2f%%** | **%d** | **%d** |\n\n",
		totVariants, totSiftHits, totSiftRate,
		totSiftPred["deleterious"],
		totSiftPred["tolerated"]))

	// --- PolyPhen-2 Coverage ---
	sb.WriteString("## PolyPhen-2 Coverage\n\n")
	sb.WriteString("| Study | Variants | PP2 Hits | Hit Rate | deleterious | possibly_damaging | benign |\n")
	sb.WriteString("|-------|----------|----------|----------|-------------|-------------------|--------|\n")

	for _, r := range results {
		hitRate := 0.0
		if r.variants > 0 {
			hitRate = float64(r.pp2Hits) / float64(r.variants) * 100
		}
		sb.WriteString(fmt.Sprintf("| %s | %d | %d | %.2f%% | %d | %d | %d |\n",
			r.name, r.variants, r.pp2Hits, hitRate,
			r.pp2PredCounts["deleterious"],
			r.pp2PredCounts["possibly_damaging"],
			r.pp2PredCounts["benign"]))
	}
	totPP2Rate := 0.0
	if totVariants > 0 {
		totPP2Rate = float64(totPP2Hits) / float64(totVariants) * 100
	}
	sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%d** | **%.2f%%** | **%d** | **%d** | **%d** |\n\n",
		totVariants, totPP2Hits, totPP2Rate,
		totPP2Pred["deleterious"],
		totPP2Pred["possibly_damaging"],
		totPP2Pred["benign"]))

	// --- dbSNP Coverage ---
	sb.WriteString("## dbSNP Coverage\n\n")
	sb.WriteString("| Study | Variants | dbSNP Hits | Hit Rate |\n")
	sb.WriteString("|-------|----------|------------|----------|\n")

	for _, r := range results {
		hitRate := 0.0
		if r.variants > 0 {
			hitRate = float64(r.dbsnpHits) / float64(r.variants) * 100
		}
		sb.WriteString(fmt.Sprintf("| %s | %d | %d | %.2f%% |\n",
			r.name, r.variants, r.dbsnpHits, hitRate))
	}
	totDbsnpRate := 0.0
	if totVariants > 0 {
		totDbsnpRate = float64(totDbsnpHits) / float64(totVariants) * 100
	}
	sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%d** | **%.2f%%** |\n\n",
		totVariants, totDbsnpHits, totDbsnpRate))

	// --- SIGNAL Coverage ---
	sb.WriteString("## SIGNAL Coverage\n\n")
	sb.WriteString("*Note: SIGNAL data uses GRCh37 coordinates. TCGA GDC MAFs use GRCh38, so few hits are expected.*\n\n")
	sb.WriteString("| Study | Variants | SIGNAL Hits | Hit Rate |\n")
	sb.WriteString("|-------|----------|-------------|----------|\n")

	for _, r := range results {
		hitRate := 0.0
		if r.variants > 0 {
			hitRate = float64(r.sigHits) / float64(r.variants) * 100
		}
		sb.WriteString(fmt.Sprintf("| %s | %d | %d | %.2f%% |\n",
			r.name, r.variants, r.sigHits, hitRate))
	}
	totSigRate := 0.0
	if totVariants > 0 {
		totSigRate = float64(totSigHits) / float64(totVariants) * 100
	}
	sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%d** | **%.2f%%** |\n\n",
		totVariants, totSigHits, totSigRate))

	// --- Genomic Index Performance ---
	sb.WriteString("## Genomic Index Lookup Performance\n\n")
	sb.WriteString("| Study | Variants | Base Time | Index Lookup Time | Overhead | Lookups/sec |\n")
	sb.WriteString("|-------|----------|-----------|-------------------|----------|-------------|\n")

	var totBaseTime, totGILookupTime time.Duration
	for _, r := range results {
		overhead := 0.0
		if r.baseTime > 0 {
			overhead = float64(r.giLookupTime) / float64(r.baseTime) * 100
		}
		lookupsPerSec := 0.0
		if r.giLookupTime > 0 {
			lookupsPerSec = float64(r.variants) / r.giLookupTime.Seconds()
		}
		sb.WriteString(fmt.Sprintf("| %s | %d | %s | %s | %.1f%% | %.0f |\n",
			r.name, r.variants,
			r.baseTime.Round(time.Millisecond),
			r.giLookupTime.Round(time.Millisecond),
			overhead, lookupsPerSec))
		totBaseTime += r.baseTime
		totGILookupTime += r.giLookupTime
	}
	totOverhead := 0.0
	if totBaseTime > 0 {
		totOverhead = float64(totGILookupTime) / float64(totBaseTime) * 100
	}
	totLookupsPerSec := 0.0
	if totGILookupTime > 0 {
		totLookupsPerSec = float64(totVariants) / totGILookupTime.Seconds()
	}
	sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%s** | **%s** | **%.1f%%** | **%.0f** |\n",
		totVariants,
		totBaseTime.Round(time.Millisecond),
		totGILookupTime.Round(time.Millisecond),
		totOverhead, totLookupsPerSec))

	if err := os.WriteFile(path, []byte(sb.String()), 0644); err != nil {
		t.Fatalf("write report: %v", err)
	}
	t.Logf("report written to %s", path)
}
