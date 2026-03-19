package input

import "testing"

func TestParseGenomicLocation(t *testing.T) {
	tests := []struct {
		name    string
		line    string
		wantGL  GenomicLocation
		wantErr bool
	}{
		{
			name: "SNP",
			line: `{"chromosome":"7","start":140453136,"end":140453136,"referenceAllele":"A","variantAllele":"T"}`,
			wantGL: GenomicLocation{Chromosome: "7", Start: 140453136, End: 140453136, ReferenceAllele: "A", VariantAllele: "T"},
		},
		{
			name: "deletion",
			line: `{"chromosome":"1","start":100,"end":102,"referenceAllele":"ATG","variantAllele":"-"}`,
			wantGL: GenomicLocation{Chromosome: "1", Start: 100, End: 102, ReferenceAllele: "ATG", VariantAllele: "-"},
		},
		{
			name: "insertion",
			line: `{"chromosome":"1","start":100,"end":101,"referenceAllele":"-","variantAllele":"CGA"}`,
			wantGL: GenomicLocation{Chromosome: "1", Start: 100, End: 101, ReferenceAllele: "-", VariantAllele: "CGA"},
		},
		{
			name:    "missing chromosome",
			line:    `{"start":100,"end":100,"referenceAllele":"A","variantAllele":"T"}`,
			wantErr: true,
		},
		{
			name:    "invalid JSON",
			line:    `not json`,
			wantErr: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gl, err := ParseGenomicLocation([]byte(tt.line))
			if (err != nil) != tt.wantErr {
				t.Fatalf("err=%v, wantErr=%v", err, tt.wantErr)
			}
			if tt.wantErr {
				return
			}
			if gl != tt.wantGL {
				t.Errorf("got %+v, want %+v", gl, tt.wantGL)
			}
		})
	}
}

func TestToVariant(t *testing.T) {
	tests := []struct {
		name    string
		gl      GenomicLocation
		wantRef string
		wantAlt string
		wantPos int64
	}{
		{
			name:    "SNP",
			gl:      GenomicLocation{Chromosome: "7", Start: 100, ReferenceAllele: "A", VariantAllele: "T"},
			wantRef: "A", wantAlt: "T", wantPos: 100,
		},
		{
			name:    "deletion dash becomes empty",
			gl:      GenomicLocation{Chromosome: "1", Start: 100, ReferenceAllele: "ATG", VariantAllele: "-"},
			wantRef: "ATG", wantAlt: "", wantPos: 100,
		},
		{
			name:    "insertion dash becomes empty",
			gl:      GenomicLocation{Chromosome: "1", Start: 100, ReferenceAllele: "-", VariantAllele: "CGA"},
			wantRef: "", wantAlt: "CGA", wantPos: 100,
		},
		{
			name:    "chr prefix stripped",
			gl:      GenomicLocation{Chromosome: "chr7", Start: 100, ReferenceAllele: "A", VariantAllele: "T"},
			wantRef: "A", wantAlt: "T", wantPos: 100,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := tt.gl.ToVariant()
			if v.Ref != tt.wantRef {
				t.Errorf("Ref=%q, want %q", v.Ref, tt.wantRef)
			}
			if v.Alt != tt.wantAlt {
				t.Errorf("Alt=%q, want %q", v.Alt, tt.wantAlt)
			}
			if v.Pos != tt.wantPos {
				t.Errorf("Pos=%d, want %d", v.Pos, tt.wantPos)
			}
		})
	}
}

func TestFormatInput(t *testing.T) {
	gl := GenomicLocation{Chromosome: "7", Start: 140453136, End: 140453136, ReferenceAllele: "A", VariantAllele: "T"}
	want := "7,140453136,140453136,A,T"
	if got := gl.FormatInput(); got != want {
		t.Errorf("FormatInput()=%q, want %q", got, want)
	}
}
