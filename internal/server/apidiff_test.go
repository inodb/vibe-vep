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

	cases := apiCases()
	t.Logf("%d cases, mode=%s, base=%s, dir=%s", len(cases), mode, base, dir)

	switch mode {
	case "capture":
		if err := os.MkdirAll(dir, 0o755); err != nil {
			t.Fatal(err)
		}
		captured, failed := 0, 0
		for _, c := range cases {
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
		same, differ, missing := 0, 0, 0
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
			if bytes.Equal(got, want) {
				same++
				continue
			}
			differ++
			t.Errorf("%s (%s %s)\n%s", c.Name, c.Method, c.Path, firstDifference(want, got))
		}
		t.Logf("identical %d, differing %d, not captured %d", same, differ, missing)

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
	// Marshal sorts object keys, so only genuine differences survive.
	canonical, err := json.MarshalIndent(v, "", "  ")
	if err != nil {
		return []byte(fmt.Sprintf("status=%d\n%s", status, body))
	}
	return []byte(fmt.Sprintf("status=%d\n%s\n", status, canonical))
}

// sortTopLevel orders a top-level array by one of its fields, for responses
// whose order carries no meaning.
func sortTopLevel(v any, key string) any {
	arr, ok := v.([]any)
	if !ok {
		return v
	}
	sort.SliceStable(arr, func(i, j int) bool {
		return fieldString(arr[i], key) < fieldString(arr[j], key)
	})
	return arr
}

func fieldString(v any, key string) string {
	m, ok := v.(map[string]any)
	if !ok {
		return ""
	}
	s, _ := m[key].(string)
	return s
}

// firstDifference reports where two canonical responses first diverge, since
// printing two whole annotations helps nobody.
func firstDifference(want, got []byte) string {
	wantLines := strings.Split(string(want), "\n")
	gotLines := strings.Split(string(got), "\n")
	for i := 0; i < len(wantLines) && i < len(gotLines); i++ {
		if wantLines[i] != gotLines[i] {
			return fmt.Sprintf("line %d\n  oracle: %s\n  ours:   %s",
				i+1, strings.TrimSpace(wantLines[i]), strings.TrimSpace(gotLines[i]))
		}
	}
	return fmt.Sprintf("one response has %d lines, the other %d",
		len(wantLines), len(gotLines))
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
