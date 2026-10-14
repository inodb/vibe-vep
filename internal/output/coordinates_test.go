package output

import (
	"bytes"
	"encoding/json"
	"strings"
	"testing"

	"github.com/inodb/vibe-vep/internal/vcf"
)

// Ensembl VEP reports an insertion between bases N and N+1 as start=N+1, end=N.
// genome-nexus recovers the MAF convention with min(start,end)/max(start,end),
// so emitting (N, N) leaves nothing to swap and every insertion ends up with a
// MAF End_Position one base short.
func TestVEPCoordinates(t *testing.T) {
	tests := []struct {
		name               string
		ref, alt           string
		pos                int64
		wantStart, wantEnd int64
	}{
		{"SNV", "C", "T", 100, 100, 100},
		{"insertion", "", "AA", 100, 101, 100},
		{"single-base deletion", "G", "", 100, 100, 100},
		{"multi-base deletion", "CTGTCTG", "", 100, 100, 106},
		{"delins", "GG", "AA", 100, 100, 101},
	}
	for _, tc := range tests {
		v := &vcf.Variant{Chrom: "1", Pos: tc.pos, Ref: tc.ref, Alt: tc.alt}
		start, end := vepCoordinates(v)
		if start != tc.wantStart || end != tc.wantEnd {
			t.Errorf("%s: vepCoordinates = (%d,%d), want (%d,%d)", tc.name, start, end, tc.wantStart, tc.wantEnd)
		}
	}
}

// The MAF convention is the un-inverted one: an insertion spans (N, N+1).
func TestMAFCoordinates(t *testing.T) {
	tests := []struct {
		name               string
		ref, alt           string
		pos                int64
		wantStart, wantEnd int64
	}{
		{"SNV", "C", "T", 100, 100, 100},
		{"insertion", "", "AA", 100, 100, 101},
		{"multi-base deletion", "CTGTCTG", "", 100, 100, 106},
	}
	for _, tc := range tests {
		v := &vcf.Variant{Chrom: "1", Pos: tc.pos, Ref: tc.ref, Alt: tc.alt}
		start, end := mafCoordinates(v)
		if start != tc.wantStart || end != tc.wantEnd {
			t.Errorf("%s: mafCoordinates = (%d,%d), want (%d,%d)", tc.name, start, end, tc.wantStart, tc.wantEnd)
		}
	}
}

// End-to-end through the writer genome-nexus actually consumes.
func TestJSONLWriterInsertionCoordinates(t *testing.T) {
	v := &vcf.Variant{Chrom: "13", Pos: 32910578, Ref: "", Alt: "AA"}
	_, ann := krasVariantAndAnnotation()
	ann.TranscriptID = "ENST00000544455.1"
	line := vepLineFromWriter(t, v, ann)

	var result VEPVariantAnnotation
	if err := json.Unmarshal([]byte(line), &result); err != nil {
		t.Fatalf("invalid JSON: %v", err)
	}
	if result.Start != 32910579 || result.End != 32910578 {
		t.Errorf("start=%d end=%d, want start=32910579 end=32910578 (VEP inverts insertions)",
			result.Start, result.End)
	}
	// What genome-nexus derives from it.
	mafStart, mafEnd := min64(result.Start, result.End), max64(result.Start, result.End)
	if mafStart != 32910578 || mafEnd != 32910579 {
		t.Errorf("MAF coords = (%d,%d), want (32910578,32910579)", mafStart, mafEnd)
	}
}

func min64(a, b int64) int64 {
	if a < b {
		return a
	}
	return b
}
func max64(a, b int64) int64 {
	if a > b {
		return a
	}
	return b
}

// VEP reports how far it 3'-shifted an indel as hgvs_offset, and omits the
// field when there was no shift. genome-nexus reads it straight off the
// transcript consequence for the MAF HGVS_Offset column.
func TestHGVSOffsetEmission(t *testing.T) {
	v, ann := krasVariantAndAnnotation()
	ann.HGVSOffset = 5

	var result VEPVariantAnnotation
	if err := json.Unmarshal([]byte(vepLineFromWriter(t, v, ann)), &result); err != nil {
		t.Fatalf("invalid JSON: %v", err)
	}
	if got := result.TranscriptConsequences[0].HGVSOffset; got != 5 {
		t.Errorf("hgvs_offset=%d, want 5", got)
	}
}

func TestHGVSOffsetOmittedWhenZero(t *testing.T) {
	v, ann := krasVariantAndAnnotation()
	ann.HGVSOffset = 0

	if line := vepLineFromWriter(t, v, ann); strings.Contains(line, "hgvs_offset") {
		t.Errorf("expected hgvs_offset omitted when zero, got: %s", line)
	}
}

// genome-nexus breaks a tie between equally-severe canonical candidates with
// bestCandidates.get(0), so the order transcripts are emitted in decides which
// one lands in the MAF. Ensembl VEP orders them by ascending transcript ID
// (22/22 real VEP fixtures, no exceptions).
func TestVEPTranscriptsSortedAscending(t *testing.T) {
	v, ann := krasVariantAndAnnotation()
	ann.TranscriptID = "ENST00000579755.1"
	_, ann2 := krasVariantAndAnnotation()
	ann2.TranscriptID = "ENST00000304494.1"

	var buf bytes.Buffer
	w := NewJSONLWriter(&buf, "ensembl-vep-jsonl", "GRCh37")
	w.SetInput("9,21971091,21971092,-,C")
	if err := w.Write(v, ann); err != nil {
		t.Fatal(err)
	}
	if err := w.Write(v, ann2); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}
	var result VEPVariantAnnotation
	if err := json.Unmarshal([]byte(buf.String()), &result); err != nil {
		t.Fatalf("invalid JSON: %v", err)
	}
	ids := make([]string, len(result.TranscriptConsequences))
	for i, tc := range result.TranscriptConsequences {
		ids[i] = tc.TranscriptID
	}
	if len(ids) != 2 || ids[0] != "ENST00000304494" || ids[1] != "ENST00000579755" {
		t.Errorf("transcript order = %v, want ascending [ENST00000304494 ENST00000579755]", ids)
	}
}
