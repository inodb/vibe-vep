package cache

import (
	"testing"

	"github.com/Sereal/Sereal/Go/sereal"
)

func TestIsSereal(t *testing.T) {
	tests := []struct {
		name     string
		data     []byte
		expected bool
	}{
		{
			name:     "standard magic",
			data:     []byte{0x3D, 0x73, 0x72, 0x6C, 0x01, 0x02, 0x03},
			expected: true,
		},
		{
			name:     "high-bit magic",
			data:     []byte{0x3D, 0xF3, 0x72, 0x6C, 0x01, 0x02, 0x03},
			expected: true,
		},
		{
			name:     "JSON data",
			data:     []byte(`[{"id": "test"}]`),
			expected: false,
		},
		{
			name:     "empty data",
			data:     []byte{},
			expected: false,
		},
		{
			name:     "short data",
			data:     []byte{0x3D, 0x73, 0x72},
			expected: false,
		},
		{
			name:     "random data",
			data:     []byte{0x00, 0x01, 0x02, 0x03, 0x04},
			expected: false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := IsSereal(tt.data)
			if got != tt.expected {
				t.Errorf("IsSereal() = %v, want %v", got, tt.expected)
			}
		})
	}
}

func TestConvertTranscript(t *testing.T) {
	tests := []struct {
		name    string
		raw     interface{}
		chrom   string
		wantID  string
		wantErr bool
	}{
		{
			name: "basic transcript",
			raw: map[string]interface{}{
				"stable_id":        "ENST00000311936",
				"_gene_stable_id":  "ENSG00000133703",
				"_gene_symbol":     "KRAS",
				"start":            int64(25205246),
				"end":              int64(25250929),
				"strand":           int64(-1),
				"biotype":          "protein_coding",
				"is_canonical":     int64(1),
				"translateable_seq": "ATG",
			},
			chrom:   "12",
			wantID:  "ENST00000311936",
			wantErr: false,
		},
		{
			name: "transcript with bool canonical",
			raw: map[string]interface{}{
				"stable_id":       "ENST00000556131",
				"_gene_stable_id": "ENSG00000133703",
				"_gene_symbol":    "KRAS",
				"start":           int64(25205246),
				"end":             int64(25250929),
				"strand":          int64(-1),
				"biotype":         "protein_coding",
				"is_canonical":    false,
			},
			chrom:   "12",
			wantID:  "ENST00000556131",
			wantErr: false,
		},
		{
			name: "transcript with nested translation",
			raw: map[string]interface{}{
				"stable_id": "ENST00000311936",
				"start":     int64(25205246),
				"end":       int64(25250929),
				"strand":    int64(-1),
				"biotype":   "protein_coding",
				"translation": map[string]interface{}{
					"seq": "MTEYKLVVVGAGGVGKSALTI",
				},
			},
			chrom:   "12",
			wantID:  "ENST00000311936",
			wantErr: false,
		},
		{
			name:    "nil input",
			raw:     nil,
			chrom:   "12",
			wantID:  "",
			wantErr: true,
		},
		{
			name:    "non-map input",
			raw:     "invalid",
			chrom:   "12",
			wantID:  "",
			wantErr: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := convertTranscript(tt.raw, tt.chrom)
			if (err != nil) != tt.wantErr {
				t.Errorf("convertTranscript() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if tt.wantErr {
				return
			}
			if got == nil {
				if tt.wantID != "" {
					t.Errorf("convertTranscript() returned nil, want ID %v", tt.wantID)
				}
				return
			}
			if got.ID != tt.wantID {
				t.Errorf("convertTranscript() ID = %v, want %v", got.ID, tt.wantID)
			}
			if got.Chrom != tt.chrom {
				t.Errorf("convertTranscript() Chrom = %v, want %v", got.Chrom, tt.chrom)
			}
		})
	}
}

func TestConvertTranscriptFields(t *testing.T) {
	raw := map[string]interface{}{
		"stable_id":           "ENST00000311936",
		"_gene_stable_id":     "ENSG00000133703",
		"_gene_symbol":        "KRAS",
		"start":               int64(25205246),
		"end":                 int64(25250929),
		"strand":              int64(-1),
		"biotype":             "protein_coding",
		"is_canonical":        int64(1),
		"coding_region_start": int64(25209798),
		"coding_region_end":   int64(25245384),
		"translateable_seq":   "ATGACTGAATATAAACTTGTG",
		"translation": map[string]interface{}{
			"seq": "MTEYKLVVVGAGGVGKSALTI",
		},
	}

	tr, err := convertTranscript(raw, "12")
	if err != nil {
		t.Fatalf("convertTranscript() error = %v", err)
	}

	// Verify all fields
	if tr.ID != "ENST00000311936" {
		t.Errorf("ID = %v, want ENST00000311936", tr.ID)
	}
	if tr.GeneID != "ENSG00000133703" {
		t.Errorf("GeneID = %v, want ENSG00000133703", tr.GeneID)
	}
	if tr.GeneName != "KRAS" {
		t.Errorf("GeneName = %v, want KRAS", tr.GeneName)
	}
	if tr.Start != 25205246 {
		t.Errorf("Start = %v, want 25205246", tr.Start)
	}
	if tr.End != 25250929 {
		t.Errorf("End = %v, want 25250929", tr.End)
	}
	if tr.Strand != -1 {
		t.Errorf("Strand = %v, want -1", tr.Strand)
	}
	if tr.Biotype != "protein_coding" {
		t.Errorf("Biotype = %v, want protein_coding", tr.Biotype)
	}
	if !tr.IsCanonical {
		t.Error("IsCanonical = false, want true")
	}
	if tr.CDSStart != 25209798 {
		t.Errorf("CDSStart = %v, want 25209798", tr.CDSStart)
	}
	if tr.CDSEnd != 25245384 {
		t.Errorf("CDSEnd = %v, want 25245384", tr.CDSEnd)
	}
	if tr.CDSSequence != "ATGACTGAATATAAACTTGTG" {
		t.Errorf("CDSSequence = %v, want ATGACTGAATATAAACTTGTG", tr.CDSSequence)
	}
	if tr.ProteinSequence != "MTEYKLVVVGAGGVGKSALTI" {
		t.Errorf("ProteinSequence = %v, want MTEYKLVVVGAGGVGKSALTI", tr.ProteinSequence)
	}
}

func TestConvertExons(t *testing.T) {
	raw := []interface{}{
		map[string]interface{}{
			"start": int64(25245274),
			"end":   int64(25245395),
			"phase": int64(0),
		},
		map[string]interface{}{
			"start": int64(25227234),
			"end":   int64(25227412),
			"phase": int64(0),
		},
		map[string]interface{}{
			"start": int64(25225614),
			"end":   int64(25225773),
			"phase": int64(2),
		},
	}

	exons, err := convertExons(raw, 25209798, 25245384)
	if err != nil {
		t.Fatalf("convertExons() error = %v", err)
	}

	if len(exons) != 3 {
		t.Fatalf("got %d exons, want 3", len(exons))
	}

	// First exon
	if exons[0].Number != 1 {
		t.Errorf("exon[0].Number = %d, want 1", exons[0].Number)
	}
	if exons[0].Start != 25245274 {
		t.Errorf("exon[0].Start = %d, want 25245274", exons[0].Start)
	}
	if exons[0].End != 25245395 {
		t.Errorf("exon[0].End = %d, want 25245395", exons[0].End)
	}
	// CDSEnd should be capped at transcript CDSEnd
	if exons[0].CDSStart != 25245274 {
		t.Errorf("exon[0].CDSStart = %d, want 25245274", exons[0].CDSStart)
	}
	if exons[0].CDSEnd != 25245384 {
		t.Errorf("exon[0].CDSEnd = %d, want 25245384 (capped at CDS end)", exons[0].CDSEnd)
	}
	if exons[0].Frame != 0 {
		t.Errorf("exon[0].Frame = %d, want 0", exons[0].Frame)
	}

	// Third exon
	if exons[2].Frame != 2 {
		t.Errorf("exon[2].Frame = %d, want 2", exons[2].Frame)
	}
}

func TestDecodeSereal(t *testing.T) {
	// Create test data by encoding transcripts with sereal
	testTranscripts := []interface{}{
		map[string]interface{}{
			"stable_id":       "ENST00000311936",
			"_gene_stable_id": "ENSG00000133703",
			"_gene_symbol":    "KRAS",
			"start":           int64(25205246),
			"end":             int64(25250929),
			"strand":          int64(-1),
			"biotype":         "protein_coding",
			"is_canonical":    int64(1),
		},
		map[string]interface{}{
			"stable_id":       "ENST00000556131",
			"_gene_stable_id": "ENSG00000133703",
			"_gene_symbol":    "KRAS",
			"start":           int64(25205246),
			"end":             int64(25250929),
			"strand":          int64(-1),
			"biotype":         "protein_coding",
			"is_canonical":    int64(0),
		},
	}

	// Encode to sereal
	encoder := sereal.NewEncoderV3()
	data, err := encoder.Marshal(testTranscripts)
	if err != nil {
		t.Fatalf("failed to encode test data: %v", err)
	}

	// Verify it's detected as sereal
	if !IsSereal(data) {
		t.Fatal("encoded data not detected as sereal")
	}

	// Decode
	transcripts, err := DecodeSereal(data, "12")
	if err != nil {
		t.Fatalf("DecodeSereal() error = %v", err)
	}

	if len(transcripts) != 2 {
		t.Fatalf("got %d transcripts, want 2", len(transcripts))
	}

	// Check first transcript
	if transcripts[0].ID != "ENST00000311936" {
		t.Errorf("transcripts[0].ID = %v, want ENST00000311936", transcripts[0].ID)
	}
	if transcripts[0].GeneName != "KRAS" {
		t.Errorf("transcripts[0].GeneName = %v, want KRAS", transcripts[0].GeneName)
	}
	if transcripts[0].Chrom != "12" {
		t.Errorf("transcripts[0].Chrom = %v, want 12", transcripts[0].Chrom)
	}
	if !transcripts[0].IsCanonical {
		t.Error("transcripts[0].IsCanonical = false, want true")
	}

	// Check second transcript
	if transcripts[1].ID != "ENST00000556131" {
		t.Errorf("transcripts[1].ID = %v, want ENST00000556131", transcripts[1].ID)
	}
	if transcripts[1].IsCanonical {
		t.Error("transcripts[1].IsCanonical = true, want false")
	}
}

func TestGetString(t *testing.T) {
	m := map[string]interface{}{
		"string": "value",
		"bytes":  []byte("bytes_value"),
		"int":    int64(123),
		"nil":    nil,
	}

	if got := getString(m, "string"); got != "value" {
		t.Errorf("getString(string) = %v, want value", got)
	}
	if got := getString(m, "bytes"); got != "bytes_value" {
		t.Errorf("getString(bytes) = %v, want bytes_value", got)
	}
	if got := getString(m, "int"); got != "" {
		t.Errorf("getString(int) = %v, want empty", got)
	}
	if got := getString(m, "missing"); got != "" {
		t.Errorf("getString(missing) = %v, want empty", got)
	}
}

func TestGetInt64(t *testing.T) {
	m := map[string]interface{}{
		"int64":   int64(123),
		"int":     int(456),
		"float64": float64(789.5),
		"uint64":  uint64(999),
		"string":  "not a number",
	}

	if got := getInt64(m, "int64"); got != 123 {
		t.Errorf("getInt64(int64) = %v, want 123", got)
	}
	if got := getInt64(m, "int"); got != 456 {
		t.Errorf("getInt64(int) = %v, want 456", got)
	}
	if got := getInt64(m, "float64"); got != 789 {
		t.Errorf("getInt64(float64) = %v, want 789", got)
	}
	if got := getInt64(m, "uint64"); got != 999 {
		t.Errorf("getInt64(uint64) = %v, want 999", got)
	}
	if got := getInt64(m, "string"); got != 0 {
		t.Errorf("getInt64(string) = %v, want 0", got)
	}
	if got := getInt64(m, "missing"); got != 0 {
		t.Errorf("getInt64(missing) = %v, want 0", got)
	}
}

func TestGetBool(t *testing.T) {
	m := map[string]interface{}{
		"true_bool":   true,
		"false_bool":  false,
		"one_int":     int64(1),
		"zero_int":    int64(0),
		"one_float":   float64(1.0),
		"zero_float":  float64(0.0),
		"one_str":     "1",
		"true_str":    "true",
		"zero_str":    "0",
		"false_str":   "false",
		"random_str":  "random",
	}

	tests := []struct {
		key  string
		want bool
	}{
		{"true_bool", true},
		{"false_bool", false},
		{"one_int", true},
		{"zero_int", false},
		{"one_float", true},
		{"zero_float", false},
		{"one_str", true},
		{"true_str", true},
		{"zero_str", false},
		{"false_str", false},
		{"random_str", false},
		{"missing", false},
	}

	for _, tt := range tests {
		t.Run(tt.key, func(t *testing.T) {
			if got := getBool(m, tt.key); got != tt.want {
				t.Errorf("getBool(%s) = %v, want %v", tt.key, got, tt.want)
			}
		})
	}
}
