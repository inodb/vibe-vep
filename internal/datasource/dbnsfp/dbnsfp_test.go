package dbnsfp

import (
	"strings"
	"testing"
)

func TestIndexColumns(t *testing.T) {
	header := strings.Split("#chr\tpos(1-based)\tref\talt\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred", "\t")
	col := IndexColumns(header,
		ColChrom, ColPos, ColRef, ColAlt,
		ColSIFTScore, ColSIFTPred, ColPP2Score, ColPP2Pred,
	)

	if col[ColChrom] != 0 {
		t.Errorf("ColChrom=%d, want 0", col[ColChrom])
	}
	if col[ColSIFTScore] != 4 {
		t.Errorf("ColSIFTScore=%d, want 4", col[ColSIFTScore])
	}
	if col[ColPP2Pred] != 7 {
		t.Errorf("ColPP2Pred=%d, want 7", col[ColPP2Pred])
	}
}

func TestParseLine(t *testing.T) {
	header := strings.Split("#chr\tpos(1-based)\tref\talt\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred", "\t")
	col := IndexColumns(header,
		ColChrom, ColPos, ColRef, ColAlt,
		ColSIFTScore, ColSIFTPred, ColPP2Score, ColPP2Pred,
	)

	tests := []struct {
		name       string
		line       string
		wantCount  int
		wantSIFT   float32
		wantSPred  string
		wantPP2    float32
		wantPPPred string
		wantOK     bool
	}{
		{
			name:       "single alt",
			line:       "chr1\t69091\tA\tC\t0.123\tD\t0.987\tD",
			wantCount:  1,
			wantSIFT:   0.123,
			wantSPred:  "deleterious",
			wantPP2:    0.987,
			wantPPPred: "deleterious",
			wantOK:     true,
		},
		{
			name:       "multi-allelic",
			line:       "chr1\t69091\tA\tC;T\t0.1;0.2\tD;T\t0.9;0.3\tD;B",
			wantCount:  2,
			wantSIFT:   0.1,
			wantSPred:  "deleterious",
			wantPP2:    0.9,
			wantPPPred: "deleterious",
			wantOK:     true,
		},
		{
			name:   "all dots",
			line:   "chr1\t69091\tA\tC\t.\t.\t.\t.",
			wantOK: false,
		},
		{
			name:      "partial data",
			line:      "chr1\t69091\tA\tC\t0.5\tT\t.\t.",
			wantCount: 1,
			wantSIFT:  0.5,
			wantSPred: "tolerated",
			wantOK:    true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			fields := strings.Split(tt.line, "\t")
			entries, _, ok := ParseLine(fields, col)
			if ok != tt.wantOK {
				t.Fatalf("ok=%v, want %v", ok, tt.wantOK)
			}
			if !ok {
				return
			}
			if len(entries) != tt.wantCount {
				t.Fatalf("got %d entries, want %d", len(entries), tt.wantCount)
			}
			e := entries[0]
			if diff := e.SIFTScore - tt.wantSIFT; diff < -0.001 || diff > 0.001 {
				t.Errorf("SIFTScore=%v, want %v", e.SIFTScore, tt.wantSIFT)
			}
			if e.SIFTPred != tt.wantSPred {
				t.Errorf("SIFTPred=%q, want %q", e.SIFTPred, tt.wantSPred)
			}
			if diff := e.PP2Score - tt.wantPP2; diff < -0.001 || diff > 0.001 {
				t.Errorf("PP2Score=%v, want %v", e.PP2Score, tt.wantPP2)
			}
			if e.PP2Pred != tt.wantPPPred {
				t.Errorf("PP2Pred=%q, want %q", e.PP2Pred, tt.wantPPPred)
			}
		})
	}
}

func TestExpandPrediction(t *testing.T) {
	tests := []struct {
		code string
		want string
	}{
		{"D", "deleterious"},
		{"T", "tolerated"},
		{"P", "possibly_damaging"},
		{"B", "benign"},
		{"unknown", "unknown"},
	}
	for _, tt := range tests {
		if got := ExpandPrediction(tt.code); got != tt.want {
			t.Errorf("ExpandPrediction(%q) = %q, want %q", tt.code, got, tt.want)
		}
	}
}
