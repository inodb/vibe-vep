package dbsnp

import "testing"

func TestParseVCFLine(t *testing.T) {
	tests := []struct {
		name    string
		line    string
		wantID  string
		wantRef string
		wantAlt string
		wantOK  bool
	}{
		{
			name:    "simple SNP",
			line:    "1\t10177\trs367896724\tA\tAC\t.\t.\tRS=367896724",
			wantID:  "rs367896724",
			wantRef: "A",
			wantAlt: "AC",
			wantOK:  true,
		},
		{
			name:    "multi-allelic takes first",
			line:    "1\t10177\trs367896724\tA\tAC,T\t.\t.\tRS=367896724",
			wantID:  "rs367896724",
			wantRef: "A",
			wantAlt: "AC",
			wantOK:  true,
		},
		{
			name:   "no RS ID",
			line:   "1\t10177\t.\tA\tAC\t.\t.\t.",
			wantOK: false,
		},
		{
			name:   "too few fields",
			line:   "1\t10177",
			wantOK: false,
		},
		{
			name:    "RefSeq accession",
			line:    "NC_000001.11\t10177\trs367896724\tA\tAC\t.\t.\tRS=367896724",
			wantID:  "rs367896724",
			wantRef: "A",
			wantAlt: "AC",
			wantOK:  true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			entry, _, ok := ParseVCFLine(tt.line)
			if ok != tt.wantOK {
				t.Fatalf("ok=%v, want %v", ok, tt.wantOK)
			}
			if !ok {
				return
			}
			if entry.ID != tt.wantID {
				t.Errorf("ID=%q, want %q", entry.ID, tt.wantID)
			}
			if entry.Ref != tt.wantRef {
				t.Errorf("Ref=%q, want %q", entry.Ref, tt.wantRef)
			}
			if entry.Alt != tt.wantAlt {
				t.Errorf("Alt=%q, want %q", entry.Alt, tt.wantAlt)
			}
		})
	}
}

func TestNormalizeChrom(t *testing.T) {
	tests := []struct {
		input string
		want  string
	}{
		{"chr1", "1"},
		{"1", "1"},
		{"NC_000001.11", "1"},
		{"NC_000022.11", "22"},
		{"NC_000023.11", "X"},
		{"NC_000024.10", "Y"},
		{"chrX", "X"},
	}

	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			got := NormalizeChrom(tt.input)
			if got != tt.want {
				t.Errorf("NormalizeChrom(%q) = %q, want %q", tt.input, got, tt.want)
			}
		})
	}
}
