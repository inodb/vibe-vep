package vcf

import "testing"

func TestVariant_IsSNV(t *testing.T) {
	tests := []struct {
		name string
		ref  string
		alt  string
		want bool
	}{
		{"A to G", "A", "G", true},
		{"G to C (KRAS G12C)", "G", "C", true},
		{"deletion", "AT", "A", false},
		{"insertion", "A", "AT", false},
		{"MNV", "AT", "GC", false},
		{"complex indel", "ATG", "A", false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := &Variant{Ref: tt.ref, Alt: tt.alt}
			if got := v.IsSNV(); got != tt.want {
				t.Errorf("IsSNV() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestVariant_IsIndel(t *testing.T) {
	tests := []struct {
		name string
		ref  string
		alt  string
		want bool
	}{
		{"SNV", "A", "G", false},
		{"deletion", "AT", "A", true},
		{"insertion", "A", "AT", true},
		{"complex deletion", "ATGC", "A", true},
		{"MNV same length", "AT", "GC", false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := &Variant{Ref: tt.ref, Alt: tt.alt}
			if got := v.IsIndel(); got != tt.want {
				t.Errorf("IsIndel() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestVariant_IsInsertion(t *testing.T) {
	tests := []struct {
		name string
		ref  string
		alt  string
		want bool
	}{
		{"SNV", "A", "G", false},
		{"deletion", "AT", "A", false},
		{"insertion", "A", "AT", true},
		{"larger insertion", "A", "ATGC", true},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := &Variant{Ref: tt.ref, Alt: tt.alt}
			if got := v.IsInsertion(); got != tt.want {
				t.Errorf("IsInsertion() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestVariant_IsDeletion(t *testing.T) {
	tests := []struct {
		name string
		ref  string
		alt  string
		want bool
	}{
		{"SNV", "A", "G", false},
		{"deletion", "AT", "A", true},
		{"insertion", "A", "AT", false},
		{"larger deletion", "ATGC", "A", true},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := &Variant{Ref: tt.ref, Alt: tt.alt}
			if got := v.IsDeletion(); got != tt.want {
				t.Errorf("IsDeletion() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestVariant_NormalizeChrom(t *testing.T) {
	tests := []struct {
		name  string
		chrom string
		want  string
	}{
		{"with chr prefix", "chr12", "12"},
		{"without chr prefix", "12", "12"},
		{"chrX", "chrX", "X"},
		{"X", "X", "X"},
		{"chrM", "chrM", "M"},
		{"MT", "MT", "MT"},
		{"chr1", "chr1", "1"},
		{"empty", "", ""},
		{"short chr", "ch", "ch"}, // too short for "chr" prefix
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := &Variant{Chrom: tt.chrom}
			if got := v.NormalizeChrom(); got != tt.want {
				t.Errorf("NormalizeChrom() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestVariant_KRASG12C(t *testing.T) {
	// Test the specific KRAS G12C variant (c.34G>T p.G12C)
	// KRAS is on reverse strand: coding G->T = genomic C->A
	v := &Variant{
		Chrom: "12",
		Pos:   25245351,
		Ref:   "C",
		Alt:   "A",
	}

	if !v.IsSNV() {
		t.Error("KRAS G12C should be classified as SNV")
	}

	if v.IsIndel() {
		t.Error("KRAS G12C should not be classified as indel")
	}

	if v.NormalizeChrom() != "12" {
		t.Errorf("Expected chromosome 12, got %s", v.NormalizeChrom())
	}
}
