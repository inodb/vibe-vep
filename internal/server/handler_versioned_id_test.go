package server

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"
)

// Ensembl ids are often written with their version, and every endpoint that
// takes one accepts both spellings. genome-nexus does not: it keys these
// lookups on the bare id, so a versioned one finds nothing and it answers 200
// with an empty body. Answering anyway is a deliberate divergence, recorded
// in the differential harness, and this is what holds the endpoints to it.
//
// The two spellings must also give the same answer, since a caller should not
// have to know which one a given endpoint prefers.
func TestVersionedTranscriptIDs(t *testing.T) {
	srv := newTestServerWithKRAS(t)
	handler := srv.Handler()

	const (
		bare      = "ENST00000311936"
		versioned = "ENST00000311936.8"
	)

	for _, c := range []struct {
		name          string
		bare, withVer string
	}{
		{
			name:    "GET /ensembl/transcript",
			bare:    "/genome-nexus/grch38/ensembl/transcript/" + bare,
			withVer: "/genome-nexus/grch38/ensembl/transcript/" + versioned,
		},
		{
			name:    "GET /cancer_hotspots/transcript",
			bare:    "/genome-nexus/grch38/cancer_hotspots/transcript/" + bare,
			withVer: "/genome-nexus/grch38/cancer_hotspots/transcript/" + versioned,
		},
		{
			name:    "GET /ptm/experimental",
			bare:    "/genome-nexus/grch38/ptm/experimental?ensemblTranscriptId=" + bare,
			withVer: "/genome-nexus/grch38/ptm/experimental?ensemblTranscriptId=" + versioned,
		},
	} {
		t.Run(c.name, func(t *testing.T) {
			wantStatus, wantBody := get(t, handler, c.bare)
			gotStatus, gotBody := get(t, handler, c.withVer)

			if wantStatus != http.StatusOK {
				t.Fatalf("the bare id did not resolve: status %d, %s", wantStatus, wantBody)
			}
			if gotStatus != wantStatus {
				t.Errorf("versioned id: status %d, bare id: %d", gotStatus, wantStatus)
			}
			if gotBody != wantBody {
				t.Errorf("versioned and bare ids gave different answers\n  bare:      %s\n  versioned: %s",
					truncate(wantBody), truncate(gotBody))
			}
		})
	}

	// The batch filter resolves them too, and reports the transcript under
	// its bare id either way.
	t.Run("POST /ensembl/transcript transcriptIds", func(t *testing.T) {
		status, body := post(t, handler, "/genome-nexus/grch38/ensembl/transcript",
			`{"transcriptIds":["`+versioned+`"]}`)
		if status != http.StatusOK {
			t.Fatalf("status %d: %s", status, body)
		}
		var found []struct {
			TranscriptID string `json:"transcriptId"`
		}
		if err := json.Unmarshal([]byte(body), &found); err != nil {
			t.Fatalf("decode: %v (%s)", err, body)
		}
		if len(found) != 1 || found[0].TranscriptID != bare {
			t.Errorf("got %v, want one entry for %s", found, bare)
		}
	})
}

func get(t *testing.T, handler http.Handler, path string) (int, string) {
	t.Helper()
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, httptest.NewRequest(http.MethodGet, path, nil))
	return w.Code, w.Body.String()
}

func post(t *testing.T, handler http.Handler, path, body string) (int, string) {
	t.Helper()
	req := httptest.NewRequest(http.MethodPost, path, strings.NewReader(body))
	req.Header.Set("Content-Type", "application/json")
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, req)
	return w.Code, w.Body.String()
}

func truncate(s string) string {
	const limit = 160
	if len(s) <= limit {
		return s
	}
	return s[:limit] + "..."
}
