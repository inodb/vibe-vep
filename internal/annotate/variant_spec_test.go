package annotate

import "testing"

func TestParseVariantSpec(t *testing.T) {
	tests := []struct {
		input    string
		wantType VariantSpecType
		wantErr  bool
		check    func(*testing.T, *VariantSpec)
	}{
		// Genomic formats
		{
			input:    "12:25245350:C:A",
			wantType: SpecGenomic,
			check: func(t *testing.T, s *VariantSpec) {
				t.Helper()
				if s.Chrom != "12" || s.Pos != 25245350 || s.Ref != "C" || s.Alt != "A" {
					t.Errorf("got %+v", s)
				}
			},
		},
		{
			input:    "chr12:25245350:C>A",
			wantType: SpecGenomic,
			check: func(t *testing.T, s *VariantSpec) {
				t.Helper()
				if s.Chrom != "12" || s.Pos != 25245350 || s.Ref != "C" || s.Alt != "A" {
					t.Errorf("got %+v", s)
				}
			},
		},
		{
			input:    "12-25245350-C-A",
			wantType: SpecGenomic,
			check: func(t *testing.T, s *VariantSpec) {
				t.Helper()
				if s.Chrom != "12" || s.Pos != 25245350 || s.Ref != "C" || s.Alt != "A" {
					t.Errorf("got %+v", s)
				}
			},
		},
		// Protein formats (single-letter)
		{
			input:    "KRAS G12C",
			wantType: SpecProtein,
			check: func(t *testing.T, s *VariantSpec) {
				t.Helper()
				if s.GeneName != "KRAS" || s.RefAA != 'G' || s.Position != 12 || s.AltAA != 'C' {
					t.Errorf("got %+v", s)
				}
			},
		},
		{
			input:    "KRAS p.G12C",
			wantType: SpecProtein,
			check: func(t *testing.T, s *VariantSpec) {
				t.Helper()
				if s.GeneName != "KRAS" || s.RefAA != 'G' || s.Position != 12 || s.AltAA != 'C' {
					t.Errorf("got %+v", s)
				}
			},
		},
		// Protein formats (three-letter)
		{
			input:    "KRAS p.Gly12Cys",
			wantType: SpecProtein,
			check: func(t *testing.T, s *VariantSpec) {
				t.Helper()
				if s.GeneName != "KRAS" || s.RefAA != 'G' || s.Position != 12 || s.AltAA != 'C' {
					t.Errorf("got %+v", s)
				}
			},
		},
		// HGVSc formats
		{
			input:    "KRAS c.35G>T",
			wantType: SpecHGVSc,
			check: func(t *testing.T, s *VariantSpec) {
				t.Helper()
				if s.TranscriptID != "KRAS" || s.CDSChange != "35G>T" {
					t.Errorf("got %+v", s)
				}
			},
		},
		{
			input:    "ENST00000311936:c.35G>T",
			wantType: SpecHGVSc,
			check: func(t *testing.T, s *VariantSpec) {
				t.Helper()
				if s.TranscriptID != "ENST00000311936" || s.CDSChange != "35G>T" {
					t.Errorf("got %+v", s)
				}
			},
		},
		// Stop codon in protein
		{
			input:    "TP53 R196*",
			wantType: SpecProtein,
			check: func(t *testing.T, s *VariantSpec) {
				t.Helper()
				if s.GeneName != "TP53" || s.RefAA != 'R' || s.Position != 196 || s.AltAA != '*' {
					t.Errorf("got %+v", s)
				}
			},
		},
		// Errors
		{input: "", wantErr: true},
		{input: "not a variant", wantErr: true},
	}

	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			spec, err := ParseVariantSpec(tt.input)
			if tt.wantErr {
				if err == nil {
					t.Fatalf("expected error, got %+v", spec)
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if spec.Type != tt.wantType {
				t.Errorf("type = %d, want %d", spec.Type, tt.wantType)
			}
			if tt.check != nil {
				tt.check(t, spec)
			}
		})
	}
}
