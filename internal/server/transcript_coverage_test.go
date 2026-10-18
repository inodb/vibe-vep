//go:build apidiff

package server

import (
	"bufio"
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
	"time"
)

// The MAF comparison checks one transcript per variant -- the canonical one,
// which is all the annotation pipeline writes out. A variant reaches many more
// than that, and callers do read them: a second opinion on a gene with two
// plausible isoforms, or the non-coding transcripts a splice variant also
// falls in.
//
// This asks a narrower question than the response diff next door, over a much
// wider sample: for each variant, does vibe-vep report every transcript
// genome-nexus does? A transcript the oracle reports and we do not is a
// failure. One we report and it does not is listed but not failed, since
// answering about more transcripts costs a caller nothing.
//
// The consequence terms of each transcript are carried along and compared
// too, because a transcript is only worth reporting if what it says is right.
// Those differences are counted rather than failed: vibe-vep is not VEP, and
// the size of that gap is what the MAF comparison measures, not this.
//
// It runs in the same two passes as TestAPIDiff, against the same oracle:
//
//	APIDIFF_MODE=capture APIDIFF_BASE=https://vep111.genomenexus.org \
//	  go test -tags apidiff ./internal/server/ -run TestTranscriptCoverage -v
//	APIDIFF_MODE=replay  APIDIFF_BASE=http://localhost:38081 \
//	  go test -tags apidiff ./internal/server/ -run TestTranscriptCoverage -v

const (
	coverageVariantsFile = "testdata/transcript_coverage_variants.txt"
	coverageCaptureFile  = "testdata/transcript_coverage_oracle.json"

	// coverageBatchSize keeps each POST within what the oracle answers
	// comfortably, and keeps a failure from costing the whole run.
	coverageBatchSize = 50
)

// coverageVariant is one variant of the corpus, with the classification it was
// sampled under so a failure says what kind of variant it was.
type coverageVariant struct {
	Location       string // "chromosome,start,end,ref,alt"
	Classification string
}

// coverageCapture is what the oracle answered: per variant, the consequence
// terms it reported for each transcript.
type coverageCapture map[string]map[string]string

func TestTranscriptCoverage(t *testing.T) {
	mode := os.Getenv("APIDIFF_MODE")
	base := os.Getenv("APIDIFF_BASE")
	if mode == "" || base == "" {
		t.Skip("set APIDIFF_MODE=capture|replay and APIDIFF_BASE=<server url>")
	}

	variants := readCoverageVariants(t)
	t.Logf("%d variants, mode=%s, base=%s", len(variants), mode, base)

	delay := 0 * time.Millisecond
	if d := os.Getenv("APIDIFF_DELAY"); d != "" {
		if parsed, err := time.ParseDuration(d); err == nil {
			delay = parsed
		}
	}

	switch mode {
	case "capture":
		captured := make(coverageCapture, len(variants))
		for start := 0; start < len(variants); start += coverageBatchSize {
			end := min(start+coverageBatchSize, len(variants))
			time.Sleep(delay)
			batch := variants[start:end]
			answers, err := transcriptsForBatch(base, batch)
			if err != nil {
				t.Fatalf("capture batch %d-%d: %v", start, end, err)
			}
			for location, terms := range answers {
				captured[location] = terms
			}
			t.Logf("captured %d/%d", end, len(variants))
		}
		encoded, err := json.MarshalIndent(captured, "", " ")
		if err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(coverageCaptureFile, append(encoded, '\n'), 0o644); err != nil {
			t.Fatal(err)
		}
		t.Logf("wrote %d variants to %s", len(captured), coverageCaptureFile)

	case "replay":
		want := readCoverageCapture(t)
		byClass := map[string]*coverageTally{}
		var (
			total, complete, incomplete, unanswered int
			reported                                int
		)
		for start := 0; start < len(variants); start += coverageBatchSize {
			end := min(start+coverageBatchSize, len(variants))
			batch := variants[start:end]
			answers, err := transcriptsForBatch(base, batch)
			if err != nil {
				t.Fatalf("replay batch %d-%d: %v", start, end, err)
			}
			for _, v := range batch {
				expected, ok := want[v.Location]
				if !ok {
					unanswered++
					continue
				}
				total++
				tally := byClass[v.Classification]
				if tally == nil {
					tally = &coverageTally{}
					byClass[v.Classification] = tally
				}
				tally.variants++

				ours := answers[v.Location]
				missing, extra := compareSets(expected, ours)
				tally.oracle += len(expected)
				tally.ours += len(ours)
				tally.missing += len(missing)
				tally.extra += len(extra)
				for id, want := range expected {
					if got, ok := ours[id]; ok && got != want {
						tally.differingTerms++
					}
				}

				if len(missing) == 0 {
					complete++
					continue
				}
				incomplete++
				if reported < 20 {
					reported++
					t.Errorf("%s (%s): %d of the oracle's %d transcripts missing: %s",
						v.Location, v.Classification, len(missing), len(expected),
						strings.Join(missing, " "))
				}
			}
		}
		if incomplete > 20 {
			t.Errorf("... and %d further variants missing transcripts", incomplete-20)
		}
		t.Logf("%d variants: %d complete, %d missing transcripts, %d not captured",
			total, complete, incomplete, unanswered)
		t.Logf("per classification:\n%s", renderCoverage(byClass))

	default:
		t.Fatalf("unknown APIDIFF_MODE %q", mode)
	}
}

// coverageTally counts one classification's variants and transcripts.
type coverageTally struct {
	variants, oracle, ours, missing, extra, differingTerms int
}

// renderCoverage lays the tally out, worst first, so a classification we are
// systematically short on stands out from one with a single odd variant.
func renderCoverage(byClass map[string]*coverageTally) string {
	names := make([]string, 0, len(byClass))
	for name := range byClass {
		names = append(names, name)
	}
	sort.Slice(names, func(i, j int) bool {
		a, b := byClass[names[i]], byClass[names[j]]
		if a.missing != b.missing {
			return a.missing > b.missing
		}
		return names[i] < names[j]
	})

	var b strings.Builder
	fmt.Fprintf(&b, "  %-24s %8s %9s %7s %8s %6s %10s\n",
		"classification", "variants", "oracle tx", "ours", "missing", "extra", "diff terms")
	var total coverageTally
	for _, name := range names {
		c := byClass[name]
		fmt.Fprintf(&b, "  %-24s %8d %9d %7d %8d %6d %10d\n",
			name, c.variants, c.oracle, c.ours, c.missing, c.extra, c.differingTerms)
		total.variants += c.variants
		total.oracle += c.oracle
		total.ours += c.ours
		total.missing += c.missing
		total.extra += c.extra
		total.differingTerms += c.differingTerms
	}
	fmt.Fprintf(&b, "  %-24s %8d %9d %7d %8d %6d %10d\n",
		"TOTAL", total.variants, total.oracle, total.ours,
		total.missing, total.extra, total.differingTerms)
	return strings.TrimRight(b.String(), "\n")
}

// compareSets reports which transcripts the oracle named that ours did not,
// and the other way round.
func compareSets(oracle, ours map[string]string) (missing, extra []string) {
	for id := range oracle {
		if _, ok := ours[id]; !ok {
			missing = append(missing, id)
		}
	}
	for id := range ours {
		if _, ok := oracle[id]; !ok {
			extra = append(extra, id)
		}
	}
	sort.Strings(missing)
	sort.Strings(extra)
	return missing, extra
}

// transcriptsForBatch annotates a batch and returns, per variant, the
// consequence terms reported for each of its transcripts.
func transcriptsForBatch(base string, batch []coverageVariant) (map[string]map[string]string, error) {
	locations := make([]map[string]any, 0, len(batch))
	for _, v := range batch {
		parsed, err := parseCommaGenomicLocation(v.Location)
		if err != nil {
			return nil, fmt.Errorf("corpus entry %q: %w", v.Location, err)
		}
		locations = append(locations, map[string]any{
			"chromosome":      parsed.Chromosome,
			"start":           parsed.Start,
			"end":             parsed.End,
			"referenceAllele": parsed.ReferenceAllele,
			"variantAllele":   parsed.VariantAllele,
		})
	}
	body, err := json.Marshal(locations)
	if err != nil {
		return nil, err
	}

	req, err := http.NewRequest(http.MethodPost,
		strings.TrimRight(base, "/")+"/annotation/genomic", bytes.NewReader(body))
	if err != nil {
		return nil, err
	}
	req.Header.Set("Content-Type", "application/json")

	client := &http.Client{Timeout: 300 * time.Second}
	resp, err := client.Do(req)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()
	payload, err := io.ReadAll(resp.Body)
	if err != nil {
		return nil, err
	}
	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("status %d: %s", resp.StatusCode, firstBytes(payload, 200))
	}

	var answers []struct {
		OriginalVariantQuery   string `json:"originalVariantQuery"`
		TranscriptConsequences []struct {
			TranscriptID     string   `json:"transcript_id"`
			ConsequenceTerms []string `json:"consequence_terms"`
		} `json:"transcript_consequences"`
	}
	if err := json.Unmarshal(payload, &answers); err != nil {
		return nil, fmt.Errorf("decode batch: %w (%s)", err, firstBytes(payload, 200))
	}

	// The answers come back in the order asked, which is how each is matched
	// to its question: originalVariantQuery restates the location only when
	// the caller sent one as text, and this endpoint takes objects.
	if len(answers) != len(batch) {
		return nil, fmt.Errorf("asked about %d variants, got %d answers", len(batch), len(answers))
	}
	out := make(map[string]map[string]string, len(batch))
	for i, answer := range answers {
		terms := make(map[string]string, len(answer.TranscriptConsequences))
		for _, tc := range answer.TranscriptConsequences {
			if tc.TranscriptID == "" {
				continue
			}
			sorted := append([]string(nil), tc.ConsequenceTerms...)
			sort.Strings(sorted)
			terms[tc.TranscriptID] = strings.Join(sorted, ",")
		}
		out[batch[i].Location] = terms
	}
	return out, nil
}

func firstBytes(b []byte, n int) string {
	if len(b) <= n {
		return string(b)
	}
	return string(b[:n]) + "..."
}

// readCoverageVariants reads the checked-in corpus.
func readCoverageVariants(t *testing.T) []coverageVariant {
	t.Helper()
	f, err := os.Open(coverageVariantsFile)
	if err != nil {
		t.Fatalf("open corpus: %v", err)
	}
	defer f.Close()

	var variants []coverageVariant
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		location, classification, _ := strings.Cut(line, "\t")
		variants = append(variants, coverageVariant{
			Location:       location,
			Classification: classification,
		})
	}
	if err := scanner.Err(); err != nil {
		t.Fatalf("read corpus: %v", err)
	}
	if len(variants) == 0 {
		t.Fatalf("%s holds no variants", coverageVariantsFile)
	}
	return variants
}

// readCoverageCapture reads what the oracle answered.
func readCoverageCapture(t *testing.T) coverageCapture {
	t.Helper()
	raw, err := os.ReadFile(coverageCaptureFile)
	if err != nil {
		t.Fatalf("read %s: %v (capture against the oracle first)",
			filepath.Base(coverageCaptureFile), err)
	}
	var captured coverageCapture
	if err := json.Unmarshal(raw, &captured); err != nil {
		t.Fatalf("decode capture: %v", err)
	}
	return captured
}
