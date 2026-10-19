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
	"sort"
	"strconv"
	"strings"
	"sync"
	"testing"
	"time"
)

// The deepest of the three comparisons, and the widest it can be run.
//
// The MAF comparison checks about twenty columns of the one transcript
// genome-nexus picked as canonical. TestTranscriptCoverage checks that every
// transcript comes back, and what its consequence terms say. This checks
// every field of every transcript of every variant.
//
// It drives both servers in the same run rather than replaying a capture,
// because the capture would not fit: a variant's annotation is about 20 KB,
// so the full cohort is several gigabytes. vibe-vep is local and the oracle
// is remote, so there is no contention in asking them at once, and nothing
// has to be stored but the differences.
//
//	ORACLE_BASE=https://vep111.genomenexus.org VIBE_BASE=http://localhost:38081 \
//	  PARITY_LIMIT=1000 go test -tags apidiff ./internal/server/ \
//	  -run TestAnnotationFieldParity -v -timeout 4h
//
// PARITY_VARIANTS names a variant list (one "chrom,start,end,ref,alt" per
// line, anything after a tab ignored); PARITY_LIMIT caps how many are used,
// 0 for all; PARITY_REPORT names a file for the per-field detail.
func TestAnnotationFieldParity(t *testing.T) {
	oracleBase := envOr("ORACLE_BASE", "")
	vibeBase := envOr("VIBE_BASE", "")
	if oracleBase == "" || vibeBase == "" {
		t.Skip("set ORACLE_BASE and VIBE_BASE to the two servers to compare")
	}

	variants := readParityVariants(t,
		envOr("PARITY_VARIANTS", coverageVariantsFile),
		atoiOr(envOr("PARITY_LIMIT", "1000"), 1000))

	delay := parseDelay(envOr("PARITY_DELAY", "0s"))
	batchSize := atoiOr(envOr("PARITY_BATCH", "50"), 50)

	t.Logf("%d variants, batches of %d, oracle=%s vibe=%s",
		len(variants), batchSize, oracleBase, vibeBase)

	report := newParityReport()
	start := time.Now()
	for from := 0; from < len(variants); from += batchSize {
		to := min(from+batchSize, len(variants))
		batch := variants[from:to]

		time.Sleep(delay)
		oracle, ours, err := annotateBothSides(oracleBase, vibeBase, batch)
		if err != nil {
			t.Fatalf("batch %d-%d: %v", from, to, err)
		}
		for i := range batch {
			report.compareVariant(batch[i], oracle[i], ours[i])
		}
		if to%2000 == 0 || to == len(variants) {
			t.Logf("%d/%d variants, %s elapsed", to, len(variants), time.Since(start).Round(time.Second))
		}
	}

	t.Log(report.summary())
	if path := envOr("PARITY_REPORT", ""); path != "" {
		if err := os.WriteFile(path, []byte(report.detail()), 0o644); err != nil {
			t.Fatalf("write report: %v", err)
		}
		t.Logf("per-field detail written to %s", path)
	}

	// A transcript the oracle reports and we do not is the one thing here
	// that is unambiguously wrong. Field differences are the known gap
	// between vibe-vep and VEP, counted rather than failed.
	if report.missingTranscripts > 0 {
		t.Errorf("%d transcripts reported by the oracle are missing from vibe-vep, across %d variants",
			report.missingTranscripts, report.variantsMissingTranscripts)
	}
}

// parityReport accumulates what differed, by field, across every variant.
type parityReport struct {
	variants                   int
	oracleTranscripts          int
	ourTranscripts             int
	missingTranscripts         int
	extraTranscripts           int
	variantsMissingTranscripts int
	comparedFields             int
	differingFields            int
	identicalVariants          int

	byField  map[string]int
	examples map[string]string
}

func newParityReport() *parityReport {
	return &parityReport{
		byField:  map[string]int{},
		examples: map[string]string{},
	}
}

// compareVariant walks one variant's two answers.
func (r *parityReport) compareVariant(v, oracle, ours map[string]any) {
	r.variants++
	before := r.differingFields

	// The transcript-level comparison, which is the point of this test.
	oracleTx := transcriptsByID(oracle)
	ourTx := transcriptsByID(ours)
	r.oracleTranscripts += len(oracleTx)
	r.ourTranscripts += len(ourTx)

	missing := 0
	for id, otc := range oracleTx {
		gtc, ok := ourTx[id]
		if !ok {
			missing++
			r.note("transcript_consequences[] MISSING", v["location"].(string)+" "+id)
			continue
		}
		pruneOracleGaps(otc, gtc)
		r.compareMaps("transcript_consequences[]", otc, gtc, v["location"].(string)+" "+id)
	}
	for id := range ourTx {
		if _, ok := oracleTx[id]; !ok {
			r.extraTranscripts++
		}
	}
	if missing > 0 {
		r.missingTranscripts += missing
		r.variantsMissingTranscripts++
	}

	// The response around them: the scalar fields, and the summary that names
	// the canonical transcript the MAF would have reported.
	pruneOracleGaps(oracle, ours)
	r.compareScalars("", oracle, ours, v["location"].(string))
	r.compareSummaries(oracle, ours, v["location"].(string))

	if r.differingFields == before && missing == 0 {
		r.identicalVariants++
	}
}

// compareSummaries compares annotation_summary, where present on both sides.
func (r *parityReport) compareSummaries(oracle, ours map[string]any, where string) {
	os_, ok1 := oracle["annotation_summary"].(map[string]any)
	gs, ok2 := ours["annotation_summary"].(map[string]any)
	if !ok1 || !ok2 {
		return
	}
	r.compareScalars("annotation_summary", os_, gs, where)

	ots, ok1 := os_["transcriptConsequenceSummary"].(map[string]any)
	gts, ok2 := gs["transcriptConsequenceSummary"].(map[string]any)
	if ok1 && ok2 {
		pruneOracleGaps(ots, gts)
		r.compareMaps("annotation_summary.transcriptConsequenceSummary", ots, gts, where)
	}
}

// compareScalars compares the non-object, non-array fields of two objects,
// which is what is left of a response once its lists are handled on their own.
func (r *parityReport) compareScalars(prefix string, oracle, ours map[string]any, where string) {
	for key, want := range oracle {
		switch want.(type) {
		case map[string]any, []any:
			continue
		}
		r.compareField(prefix, key, want, ours[key], okFlag(ours, key), where)
	}
}

// compareMaps compares every field of two objects, including their lists.
func (r *parityReport) compareMaps(prefix string, oracle, ours map[string]any, where string) {
	for key, want := range oracle {
		r.compareField(prefix, key, want, ours[key], okFlag(ours, key), where)
	}
	for key := range ours {
		if _, present := oracle[key]; !present {
			r.comparedFields++
			r.differingFields++
			r.note(prefix+"."+key+" (extra in ours)", where)
		}
	}
}

func (r *parityReport) compareField(prefix, key string, want, got any, present bool, where string) {
	path := key
	if prefix != "" {
		path = prefix + "." + key
	}
	r.comparedFields++
	if !present {
		r.differingFields++
		r.note(path+" (missing from ours)", where)
		return
	}
	if !sameValue(want, got) {
		r.differingFields++
		r.note(path, fmt.Sprintf("%s  oracle=%v  ours=%v", where, render(want), render(got)))
	}
}

func (r *parityReport) note(field, example string) {
	r.byField[field]++
	if _, seen := r.examples[field]; !seen {
		r.examples[field] = example
	}
}

func (r *parityReport) summary() string {
	var b strings.Builder
	fmt.Fprintf(&b, "\n  variants              %d (%d identical in every field)\n",
		r.variants, r.identicalVariants)
	fmt.Fprintf(&b, "  transcripts           oracle %d, ours %d, missing %d, extra %d\n",
		r.oracleTranscripts, r.ourTranscripts, r.missingTranscripts, r.extraTranscripts)
	fmt.Fprintf(&b, "  fields compared       %d, differing %d (%.4f%%)\n",
		r.comparedFields, r.differingFields, percent(r.differingFields, r.comparedFields))
	fmt.Fprintf(&b, "\n%s", rankPaths(r.byField))
	return b.String()
}

// detail is the report written to disk: every differing field, with how often
// and one example to go and look at.
func (r *parityReport) detail() string {
	type row struct {
		field string
		n     int
	}
	rows := make([]row, 0, len(r.byField))
	for field, n := range r.byField {
		rows = append(rows, row{field, n})
	}
	sort.Slice(rows, func(i, j int) bool {
		if rows[i].n != rows[j].n {
			return rows[i].n > rows[j].n
		}
		return rows[i].field < rows[j].field
	})

	var b strings.Builder
	fmt.Fprintf(&b, "%s\n\n", strings.TrimSpace(r.summary()))
	fmt.Fprintf(&b, "%-8s %-62s %s\n", "count", "field", "first example")
	for _, row := range rows {
		fmt.Fprintf(&b, "%-8d %-62s %s\n", row.n, row.field, r.examples[row.field])
	}
	return b.String()
}

// sameValue compares two decoded JSON values. Numbers are compared as
// numbers, so 600 and 600.0 agree; lists are compared in order, because the
// order VEP lists consequence terms and RefSeq accessions in is its own.
func sameValue(a, b any) bool {
	switch av := a.(type) {
	case []any:
		bv, ok := b.([]any)
		if !ok || len(av) != len(bv) {
			return false
		}
		for i := range av {
			if !sameValue(av[i], bv[i]) {
				return false
			}
		}
		return true
	case map[string]any:
		bv, ok := b.(map[string]any)
		if !ok || len(av) != len(bv) {
			return false
		}
		for k, v := range av {
			other, present := bv[k]
			if !present || !sameValue(v, other) {
				return false
			}
		}
		return true
	case float64:
		bv, ok := b.(float64)
		return ok && av == bv
	default:
		return fmt.Sprint(a) == fmt.Sprint(b)
	}
}

func render(v any) string {
	encoded, err := json.Marshal(v)
	if err != nil {
		return fmt.Sprint(v)
	}
	s := string(encoded)
	if len(s) > 70 {
		return s[:70] + "..."
	}
	return s
}

// transcriptsByID indexes a response's transcript consequences.
func transcriptsByID(annotation map[string]any) map[string]map[string]any {
	out := map[string]map[string]any{}
	list, _ := annotation["transcript_consequences"].([]any)
	for _, entry := range list {
		tc, ok := entry.(map[string]any)
		if !ok {
			continue
		}
		if id, ok := tc["transcript_id"].(string); ok && id != "" {
			out[id] = tc
		}
	}
	return out
}

// annotateBothSides asks the two servers the same batch at the same time and
// returns their answers in the order asked.
func annotateBothSides(oracleBase, vibeBase string, batch []map[string]any) (oracle, ours []map[string]any, err error) {
	locations := make([]map[string]any, 0, len(batch))
	for _, v := range batch {
		locations = append(locations, map[string]any{
			"chromosome":      v["chromosome"],
			"start":           v["start"],
			"end":             v["end"],
			"referenceAllele": v["referenceAllele"],
			"variantAllele":   v["variantAllele"],
		})
	}
	body, err := json.Marshal(locations)
	if err != nil {
		return nil, nil, err
	}

	const query = "/annotation/genomic?fields=annotation_summary&isoformOverrideSource=mskcc"
	var (
		wg                sync.WaitGroup
		oracleErr, ourErr error
	)
	wg.Add(2)
	go func() {
		defer wg.Done()
		oracle, oracleErr = postAnnotations(oracleBase+query, body)
	}()
	go func() {
		defer wg.Done()
		ours, ourErr = postAnnotations(vibeBase+query, body)
	}()
	wg.Wait()

	if oracleErr != nil {
		return nil, nil, fmt.Errorf("oracle: %w", oracleErr)
	}
	if ourErr != nil {
		return nil, nil, fmt.Errorf("vibe-vep: %w", ourErr)
	}
	if len(oracle) != len(batch) || len(ours) != len(batch) {
		return nil, nil, fmt.Errorf("asked about %d variants, oracle gave %d and vibe-vep %d",
			len(batch), len(oracle), len(ours))
	}
	return oracle, ours, nil
}

func postAnnotations(url string, body []byte) ([]map[string]any, error) {
	req, err := http.NewRequest(http.MethodPost, url, bytes.NewReader(body))
	if err != nil {
		return nil, err
	}
	req.Header.Set("Content-Type", "application/json")

	client := &http.Client{Timeout: 600 * time.Second}
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
	var answers []map[string]any
	if err := json.Unmarshal(payload, &answers); err != nil {
		return nil, fmt.Errorf("decode: %w (%s)", err, firstBytes(payload, 200))
	}
	return answers, nil
}

// readParityVariants reads the variant list, keeping at most limit of them.
func readParityVariants(t *testing.T, path string, limit int) []map[string]any {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()

	var variants []map[string]any
	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 0, 64*1024), 4*1024*1024)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		location, _, _ := strings.Cut(line, "\t")
		parsed, err := parseCommaGenomicLocation(location)
		if err != nil {
			continue
		}
		variants = append(variants, map[string]any{
			"location":        location,
			"chromosome":      parsed.Chromosome,
			"start":           parsed.Start,
			"end":             parsed.End,
			"referenceAllele": parsed.ReferenceAllele,
			"variantAllele":   parsed.VariantAllele,
		})
		if limit > 0 && len(variants) >= limit {
			break
		}
	}
	if err := scanner.Err(); err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	if len(variants) == 0 {
		t.Fatalf("%s holds no variants", path)
	}
	return variants
}

func okFlag(m map[string]any, key string) bool {
	_, ok := m[key]
	return ok
}

func percent(part, whole int) float64 {
	if whole == 0 {
		return 0
	}
	return 100 * float64(part) / float64(whole)
}

func envOr(key, fallback string) string {
	if v := os.Getenv(key); v != "" {
		return v
	}
	return fallback
}

func atoiOr(s string, fallback int) int {
	if n, err := strconv.Atoi(s); err == nil {
		return n
	}
	return fallback
}

func parseDelay(s string) time.Duration {
	if d, err := time.ParseDuration(s); err == nil {
		return d
	}
	return 0
}
