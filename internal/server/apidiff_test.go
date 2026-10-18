//go:build apidiff

package server

import (
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

// A differential harness: it drives genome-nexus and vibe-vep with the same
// requests and compares the answers. Hand-written expectations for 38
// endpoints would be both enormous and weaker -- the oracle is the service
// being replaced.
//
// The two servers cannot run at once on a 7 GB box, so this works in two
// passes. Capture against genome-nexus, shut it down, then replay against
// vibe-vep:
//
//	APIDIFF_MODE=capture APIDIFF_BASE=http://localhost:38080 \
//	  go test -tags apidiff ./internal/server/ -run TestAPIDiff -v
//	APIDIFF_MODE=replay  APIDIFF_BASE=http://localhost:38081 \
//	  go test -tags apidiff ./internal/server/ -run TestAPIDiff -v
//
// The captured answers are also reusable, so a later replay needs no Java.

const defaultAPIDiffDir = "testdata/apidiff"

func TestAPIDiff(t *testing.T) {
	mode := os.Getenv("APIDIFF_MODE")
	base := os.Getenv("APIDIFF_BASE")
	if mode == "" || base == "" {
		t.Skip("set APIDIFF_MODE=capture|replay and APIDIFF_BASE=<server url>")
	}
	dir := os.Getenv("APIDIFF_DIR")
	if dir == "" {
		dir = defaultAPIDiffDir
	}

	// Pace requests when the oracle is someone else's service.
	delay := 0 * time.Millisecond
	if d := os.Getenv("APIDIFF_DELAY"); d != "" {
		if parsed, err := time.ParseDuration(d); err == nil {
			delay = parsed
		}
	}

	cases := apiCases()
	t.Logf("%d cases, mode=%s, base=%s, dir=%s", len(cases), mode, base, dir)

	switch mode {
	case "capture":
		if err := os.MkdirAll(dir, 0o755); err != nil {
			t.Fatal(err)
		}
		captured, failed := 0, 0
		for _, c := range cases {
			time.Sleep(delay)
			body, status, err := doRequest(base, c)
			if err != nil {
				t.Errorf("%s: %v", c.Name, err)
				failed++
				continue
			}
			if err := os.WriteFile(filepath.Join(dir, slug(c.Name)+".json"),
				encodeCapture(status, body, c.SortKey), 0o644); err != nil {
				t.Fatal(err)
			}
			captured++
		}
		t.Logf("captured %d, failed %d", captured, failed)

	case "replay":
		same, differ, missing, accepted := 0, 0, 0, 0
		// Counting the paths that differ, rather than only the cases, turns a
		// wall of diffs into a short list of things to go and fix.
		paths := map[string]int{}
		for _, c := range cases {
			want, err := os.ReadFile(filepath.Join(dir, slug(c.Name)+".json"))
			if err != nil {
				missing++
				continue
			}
			body, status, err := doRequest(base, c)
			if err != nil {
				t.Errorf("%s: %v", c.Name, err)
				differ++
				continue
			}
			got := encodeCapture(status, body, c.SortKey)
			want, got = reconcile(want, got, c.SortKey)
			if bytes.Equal(got, want) {
				same++
				continue
			}
			if reason, ok := acceptedDivergences[c.Name]; ok {
				accepted++
				t.Logf("accepted divergence: %s -- %s", c.Name, reason)
				continue
			}
			differ++
			diffs := diffPaths(want, got)
			for _, d := range diffs {
				paths[d]++
			}
			t.Errorf("%s (%s %s)\n%s", c.Name, c.Method, c.Path, summarise(diffs, want, got))
		}
		t.Logf("identical %d, differing %d, accepted %d, not captured %d",
			same, differ, accepted, missing)
		if len(paths) > 0 {
			t.Logf("differing paths, most common first:\n%s", rankPaths(paths))
		}

	default:
		t.Fatalf("unknown APIDIFF_MODE %q", mode)
	}
}

// doRequest issues one case against a server.
func doRequest(base string, c apiCase) ([]byte, int, error) {
	var reqBody io.Reader
	if c.Body != "" {
		reqBody = strings.NewReader(c.Body)
	}
	req, err := http.NewRequest(c.Method, strings.TrimRight(base, "/")+c.Path, reqBody)
	if err != nil {
		return nil, 0, err
	}
	if c.Body != "" {
		req.Header.Set("Content-Type", "application/json")
	}

	client := &http.Client{Timeout: 120 * time.Second}
	resp, err := client.Do(req)
	if err != nil {
		return nil, 0, err
	}
	defer resp.Body.Close()
	body, err := io.ReadAll(resp.Body)
	return body, resp.StatusCode, err
}

// encodeCapture canonicalises a response so that formatting and the order of
// things that have no order cannot be mistaken for a difference. Every rule
// here is a deliberate concession and should be read as one.
func encodeCapture(status int, body []byte, sortKey string) []byte {
	var v any
	if err := json.Unmarshal(body, &v); err != nil {
		// Not JSON: compare the bytes as they came.
		return []byte(fmt.Sprintf("status=%d\n%s", status, body))
	}
	if sortKey != "" {
		v = sortTopLevel(v, sortKey)
	}
	sortUnorderedArrays(v)
	// Marshal sorts object keys, so only genuine differences survive.
	canonical, err := json.MarshalIndent(v, "", "  ")
	if err != nil {
		return []byte(fmt.Sprintf("status=%d\n%s", status, body))
	}
	return []byte(fmt.Sprintf("status=%d\n%s\n", status, canonical))
}

// oracleGapKeys names response fields genome-nexus models but the oracle
// deployment cannot produce, so that its silence is not read as our mistake.
// A gap is only honoured where the oracle actually omits the field: if it
// returns one and we disagree, that stays a difference.
//
// Every entry needs a reason, and the list is meant to stay short.
// acceptedDivergences names cases where vibe-vep deliberately answers
// differently from the oracle, with why. They are counted and listed on every
// run rather than hidden, so the list cannot grow without someone noticing.
//
// Each needs a reason that says why the difference is wanted, not merely that
// it exists. "We have not got to it yet" is not one of these -- that is a
// failing case.
var acceptedDivergences = map[string]string{
	"transcript ENST00000288602.6":             versionedIDReason,
	"hotspots transcript ENST00000288602.6":    versionedIDReason,
	"hotspots transcript v2 ENST00000288602.6": versionedIDReason,
	"transcript filter by ids":                 versionedIDReason + " (one of the ids asked for carries a version)",

	// vep111.genomenexus.org serves its web app for this path rather than a
	// 404 from the API, so there is no oracle answer to compare against.
	"canonical gene NOT_A_GENE": "the oracle answers this path with the frontend's HTML, not an API response",
}

// versionedIDReason is why we answer where genome-nexus does not. Accepting
// a version is a decision taken deliberately, not an oversight, and
// TestVersionedTranscriptIDs holds every endpoint to it.
const versionedIDReason = "genome-nexus keys these on the bare Ensembl id, so a versioned one finds " +
	"nothing and it answers 200 with an empty body. vibe-vep resolves it instead: a caller " +
	"holding ENST00000288602.6 is asking about a transcript we know, and matching the oracle " +
	"here would mean reproducing a defect"

var oracleGapKeys = map[string]string{
	// genome-nexus carries AlphaMissense on TranscriptConsequence,
	// TranscriptConsequenceSummary and VariantAnnotationSummary, but fills it
	// only when VEP runs with the AlphaMissense plugin. Neither
	// vep111.genomenexus.org nor www.genomenexus.org does, so both answer
	// without the field while vibe-vep reads the scores from its own index.
	"alphaMissense": "oracle runs VEP without the AlphaMissense plugin",
}

// reconcile re-encodes both sides so they are comparable: arrays with no
// meaningful order are sorted, and the declared oracle gaps are removed from
// ours wherever the oracle is silent about them. Responses that are not JSON
// are returned untouched.
//
// A capture is normalised here rather than trusted as written, because it was
// recorded by whatever version of this harness ran at the time. Changing a
// normalisation rule must change both sides at once, or the rule reads as a
// difference in every captured case.
func reconcile(want, got []byte, sortKey string) ([]byte, []byte) {
	wantStatus, wantVal, ok := decodeCapture(want)
	if !ok {
		return want, got
	}
	gotStatus, gotVal, ok := decodeCapture(got)
	if !ok {
		return want, got
	}
	if sortKey != "" {
		wantVal = sortTopLevel(wantVal, sortKey)
		gotVal = sortTopLevel(gotVal, sortKey)
	}
	sortUnorderedArrays(wantVal)
	sortUnorderedArrays(gotVal)
	pruneOracleGaps(wantVal, gotVal)
	return reencode(wantStatus, wantVal), reencode(gotStatus, gotVal)
}

// decodeCapture splits a canonical capture back into its status line and the
// JSON value beneath it.
func decodeCapture(b []byte) (string, any, bool) {
	status, rest, found := strings.Cut(string(b), "\n")
	if !found || !strings.HasPrefix(status, "status=") {
		return "", nil, false
	}
	var v any
	if err := json.Unmarshal([]byte(rest), &v); err != nil {
		return "", nil, false
	}
	return status, v, true
}

func reencode(status string, v any) []byte {
	canonical, err := json.MarshalIndent(v, "", "  ")
	if err != nil {
		return []byte(status)
	}
	return []byte(status + "\n" + string(canonical) + "\n")
}

// pruneOracleGaps walks the two responses together, deleting a declared gap
// key from ours wherever the oracle has nothing there. Arrays are walked by
// position, which holds because every list the harness compares is either
// ordered by the server or sorted by the case's SortKey.
func pruneOracleGaps(oracle, ours any) {
	switch o := oracle.(type) {
	case map[string]any:
		m, ok := ours.(map[string]any)
		if !ok {
			return
		}
		for key := range oracleGapKeys {
			if _, present := m[key]; !present {
				continue
			}
			if v, has := o[key]; !has || v == nil {
				delete(m, key)
			}
		}
		for key, ov := range o {
			if mv, has := m[key]; has {
				pruneOracleGaps(ov, mv)
			}
		}
	case []any:
		arr, ok := ours.([]any)
		if !ok {
			return
		}
		for i := 0; i < len(o) && i < len(arr); i++ {
			pruneOracleGaps(o[i], arr[i])
		}
	}
}

// unorderedArrays names response arrays whose order carries no meaning, with
// the field each is sorted by. Without this a different ordering makes every
// field of every element look wrong, burying the real differences. The
// canonical choice, which order does decide, is asserted separately -- the
// summary's own transcriptConsequenceSummary field names it.
var unorderedArrays = map[string]string{
	"transcript_consequences":        "transcript_id",
	"transcriptConsequenceSummaries": "transcriptId",
	"exons":                          "exonId",
	"pfamDomains":                    "pfamDomainId",
}

// unorderedInnerArrays names fields holding an array of arrays whose outer
// order is meaningful but whose inner elements are not ordered. An
// annotation's "hotspots" and "ptms" are both one array per transcript
// consequence, in that order, each holding records in no order at all.
var unorderedInnerArrays = map[string]bool{"annotation": true}

// sortUnorderedArrays walks a decoded response and puts those arrays in a
// deterministic order, wherever they appear.
func sortUnorderedArrays(v any) {
	switch t := v.(type) {
	case map[string]any:
		for key, child := range t {
			if field, ok := unorderedArrays[key]; ok {
				if arr, isArr := child.([]any); isArr {
					sort.SliceStable(arr, func(i, j int) bool {
						return fieldString(arr[i], field) < fieldString(arr[j], field)
					})
				}
			}
			if unorderedInnerArrays[key] {
				if outer, isArr := child.([]any); isArr {
					for _, inner := range outer {
						if arr, isArr := inner.([]any); isArr {
							sort.SliceStable(arr, func(i, j int) bool {
								return sortValue(arr[i], sortAnyElement) < sortValue(arr[j], sortAnyElement)
							})
						}
					}
				}
			}
			sortUnorderedArrays(child)
		}
	case []any:
		for _, child := range t {
			sortUnorderedArrays(child)
		}
	}
}

// sortAnyElement is the SortKey for a top-level array whose elements have no
// one identifying field. Elements are then ordered by their whole encoded
// value, which is deterministic without naming anything.
//
// It is opt-in rather than the default because plenty of arrays here are
// ordered on purpose: a batch response answers in the order it was asked, and
// reordering it would hide exactly the fault that invariant exists to catch.
const sortAnyElement = "*"

// sortTopLevel orders a top-level array by one of its fields, for responses
// whose order carries no meaning.
func sortTopLevel(v any, key string) any {
	arr, ok := v.([]any)
	if !ok {
		return v
	}
	sort.SliceStable(arr, func(i, j int) bool {
		return sortValue(arr[i], key) < sortValue(arr[j], key)
	})
	return arr
}

// sortValue renders the part of an element its order is decided by.
func sortValue(v any, key string) string {
	if key == sortAnyElement {
		encoded, err := json.Marshal(v)
		if err != nil {
			return ""
		}
		return string(encoded)
	}
	return fieldString(v, key)
}

// fieldString reads one field of an element as text. Numbers are included:
// sorting them as text is not numeric order, but it is the same order on both
// sides, which is all this needs.
func fieldString(v any, key string) string {
	m, ok := v.(map[string]any)
	if !ok {
		return ""
	}
	switch value := m[key].(type) {
	case nil:
		return ""
	case string:
		return value
	default:
		return fmt.Sprint(value)
	}
}

// diffPaths lists the JSON paths at which two canonical responses disagree,
// with array indices collapsed to [] so that the same fault across fifteen
// transcript consequences reads as one entry.
func diffPaths(want, got []byte) []string {
	_, wantVal, ok := decodeCapture(want)
	if !ok {
		return []string{"<non-JSON response>"}
	}
	_, gotVal, ok := decodeCapture(got)
	if !ok {
		return []string{"<non-JSON response>"}
	}
	seen := map[string]bool{}
	var out []string
	walkDiff("", wantVal, gotVal, func(path, kind string) {
		entry := path + " (" + kind + ")"
		if !seen[entry] {
			seen[entry] = true
			out = append(out, entry)
		}
	})
	return out
}

// walkDiff reports every path where oracle and ours disagree, saying whether
// the field is missing on one side or merely holds a different value.
func walkDiff(path string, oracle, ours any, report func(path, kind string)) {
	switch o := oracle.(type) {
	case map[string]any:
		m, ok := ours.(map[string]any)
		if !ok {
			report(path, "type")
			return
		}
		for key, ov := range o {
			sub := path + "." + key
			mv, has := m[key]
			if !has {
				report(sub, "missing from ours")
				continue
			}
			walkDiff(sub, ov, mv, report)
		}
		for key := range m {
			if _, has := o[key]; !has {
				report(path+"."+key, "extra in ours")
			}
		}
	case []any:
		arr, ok := ours.([]any)
		if !ok {
			report(path, "type")
			return
		}
		if len(o) != len(arr) {
			report(path+"[]", fmt.Sprintf("length %d vs %d", len(o), len(arr)))
		}
		for i := 0; i < len(o) && i < len(arr); i++ {
			walkDiff(path+"[]", o[i], arr[i], report)
		}
	default:
		if fmt.Sprint(oracle) != fmt.Sprint(ours) {
			report(path, "value")
		}
	}
}

// rankPaths renders the tally, commonest first.
func rankPaths(paths map[string]int) string {
	type row struct {
		path string
		n    int
	}
	rows := make([]row, 0, len(paths))
	for p, n := range paths {
		rows = append(rows, row{p, n})
	}
	sort.Slice(rows, func(i, j int) bool {
		if rows[i].n != rows[j].n {
			return rows[i].n > rows[j].n
		}
		return rows[i].path < rows[j].path
	})
	var b strings.Builder
	for _, r := range rows {
		fmt.Fprintf(&b, "  %4d  %s\n", r.n, r.path)
	}
	return strings.TrimRight(b.String(), "\n")
}

// summarise reports a case's differing paths, falling back to the line diff
// when the paths alone do not say what changed.
func summarise(diffs []string, want, got []byte) string {
	if len(diffs) == 0 {
		return describeDifferences(want, got)
	}
	var b strings.Builder
	for i, d := range diffs {
		if i == maxReportedDifferences {
			fmt.Fprintf(&b, "  ... and %d more paths\n", len(diffs)-i)
			break
		}
		fmt.Fprintf(&b, "  %s\n", d)
	}
	return strings.TrimRight(b.String(), "\n")
}

// maxReportedDifferences bounds the diff printed per case: enough to see what
// kind of difference it is, not so much that one case buries the next.
const maxReportedDifferences = 6

// describeDifferences reports where two canonical responses diverge, since
// printing two whole annotations helps nobody.
func describeDifferences(want, got []byte) string {
	wantLines := strings.Split(string(want), "\n")
	gotLines := strings.Split(string(got), "\n")

	var b strings.Builder
	shown := 0
	for i := 0; i < len(wantLines) && i < len(gotLines); i++ {
		if wantLines[i] == gotLines[i] {
			continue
		}
		if shown == maxReportedDifferences {
			fmt.Fprintf(&b, "  ... more differences past line %d\n", i+1)
			break
		}
		fmt.Fprintf(&b, "  line %d\n    oracle: %s\n    ours:   %s\n",
			i+1, strings.TrimSpace(wantLines[i]), strings.TrimSpace(gotLines[i]))
		shown++
	}
	if len(wantLines) != len(gotLines) {
		fmt.Fprintf(&b, "  oracle has %d lines, ours %d\n", len(wantLines), len(gotLines))
	}
	if b.Len() == 0 {
		return "  responses differ in bytes but not in lines"
	}
	return strings.TrimRight(b.String(), "\n")
}

// slug turns a case name into a filename.
func slug(name string) string {
	var b strings.Builder
	for _, r := range strings.ToLower(name) {
		switch {
		case r >= 'a' && r <= 'z', r >= '0' && r <= '9':
			b.WriteRune(r)
		default:
			b.WriteByte('_')
		}
	}
	return b.String()
}
