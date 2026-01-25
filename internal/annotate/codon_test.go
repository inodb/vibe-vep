package annotate

import "testing"

func TestTranslateCodon(t *testing.T) {
	tests := []struct {
		name   string
		codon  string
		want   byte
	}{
		// Standard amino acids
		{"ATG -> Met (start)", "ATG", 'M'},
		{"GGT -> Gly", "GGT", 'G'},
		{"TGT -> Cys", "TGT", 'C'},
		{"TTT -> Phe", "TTT", 'F'},
		{"AAA -> Lys", "AAA", 'K'},

		// Stop codons
		{"TAA -> Stop", "TAA", '*'},
		{"TAG -> Stop", "TAG", '*'},
		{"TGA -> Stop", "TGA", '*'},

		// Case insensitivity
		{"lowercase atg", "atg", 'M'},
		{"mixed case AtG", "AtG", 'M'},

		// Invalid codons
		{"too short", "AT", 'X'},
		{"too long", "ATGG", 'X'},
		{"invalid bases", "XYZ", 'X'},
		{"empty", "", 'X'},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := TranslateCodon(tt.codon)
			if got != tt.want {
				t.Errorf("TranslateCodon(%q) = %c, want %c", tt.codon, got, tt.want)
			}
		})
	}
}

func TestReverseComplement(t *testing.T) {
	tests := []struct {
		name string
		seq  string
		want string
	}{
		{"simple", "ATGC", "GCAT"},
		{"single base", "A", "T"},
		{"palindrome", "ATAT", "ATAT"},
		{"poly-A", "AAAA", "TTTT"},
		{"GC rich", "GCGC", "GCGC"},
		{"lowercase", "atgc", "gcat"},
		{"mixed case", "AtGc", "gCaT"},
		{"empty", "", ""},

		// KRAS G12C specific: The mutation is G>C at a specific position
		// On reverse strand, we need to complement
		{"KRAS codon 12 ref (GGT)", "GGT", "ACC"},
		{"KRAS codon 12 alt (TGT)", "TGT", "ACA"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := ReverseComplement(tt.seq)
			if got != tt.want {
				t.Errorf("ReverseComplement(%q) = %q, want %q", tt.seq, got, tt.want)
			}
		})
	}
}

func TestIsStopCodon(t *testing.T) {
	tests := []struct {
		codon string
		want  bool
	}{
		{"TAA", true},
		{"TAG", true},
		{"TGA", true},
		{"ATG", false},
		{"GGT", false},
		{"taa", true}, // case insensitive
	}

	for _, tt := range tests {
		t.Run(tt.codon, func(t *testing.T) {
			got := IsStopCodon(tt.codon)
			if got != tt.want {
				t.Errorf("IsStopCodon(%q) = %v, want %v", tt.codon, got, tt.want)
			}
		})
	}
}

func TestIsStartCodon(t *testing.T) {
	tests := []struct {
		codon string
		want  bool
	}{
		{"ATG", true},
		{"atg", true},
		{"TAA", false},
		{"GGT", false},
	}

	for _, tt := range tests {
		t.Run(tt.codon, func(t *testing.T) {
			got := IsStartCodon(tt.codon)
			if got != tt.want {
				t.Errorf("IsStartCodon(%q) = %v, want %v", tt.codon, got, tt.want)
			}
		})
	}
}

func TestGetCodon(t *testing.T) {
	cds := "ATGGGTCGATAA" // ATG GGT CGA TAA (Met Gly Arg Stop)

	tests := []struct {
		name        string
		cds         string
		codonNumber int64
		want        string
	}{
		{"codon 1", cds, 1, "ATG"},
		{"codon 2", cds, 2, "GGT"},
		{"codon 3", cds, 3, "CGA"},
		{"codon 4 (stop)", cds, 4, "TAA"},
		{"codon 0 (invalid)", cds, 0, ""},
		{"codon 5 (beyond)", cds, 5, ""},
		{"negative codon", cds, -1, ""},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := GetCodon(tt.cds, tt.codonNumber)
			if got != tt.want {
				t.Errorf("GetCodon(%q, %d) = %q, want %q", tt.cds, tt.codonNumber, got, tt.want)
			}
		})
	}
}

func TestMutateCodon(t *testing.T) {
	tests := []struct {
		name           string
		codon          string
		positionInCodon int
		newBase        byte
		want           string
	}{
		{"first position", "GGT", 0, 'T', "TGT"},   // G12C: GGT -> TGT
		{"second position", "GGT", 1, 'A', "GAT"},
		{"third position", "GGT", 2, 'A', "GGA"},
		{"invalid position", "GGT", 3, 'A', "GGT"},
		{"negative position", "GGT", -1, 'A', "GGT"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := MutateCodon(tt.codon, tt.positionInCodon, tt.newBase)
			if got != tt.want {
				t.Errorf("MutateCodon(%q, %d, %c) = %q, want %q",
					tt.codon, tt.positionInCodon, tt.newBase, got, tt.want)
			}
		})
	}
}

func TestTranslateSequence(t *testing.T) {
	tests := []struct {
		name string
		seq  string
		want string
	}{
		{"simple protein", "ATGGGTCGA", "MGR"},
		{"with stop", "ATGGGTCGATAA", "MGR*"},
		{"incomplete codon truncated", "ATGGGTCGAT", "MGR"},
		{"empty", "", ""},
		{"single codon", "ATG", "M"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := TranslateSequence(tt.seq)
			if got != tt.want {
				t.Errorf("TranslateSequence(%q) = %q, want %q", tt.seq, got, tt.want)
			}
		})
	}
}

// TestKRASG12CMutation verifies the KRAS G12C mutation logic
func TestKRASG12CMutation(t *testing.T) {
	// KRAS codon 12 on the coding strand is GGT (Glycine)
	// G12C mutation: GGT -> TGT (Glycine -> Cysteine)
	// Note: KRAS is on reverse strand, so the VCF shows G>C which
	// becomes C>G on the coding strand after complement
	// But the codon itself changes from GGT to TGT

	refCodon := "GGT"
	refAA := TranslateCodon(refCodon)
	if refAA != 'G' {
		t.Errorf("Reference codon GGT should translate to G (Gly), got %c", refAA)
	}

	// Mutate first position: G -> T gives TGT
	altCodon := MutateCodon(refCodon, 0, 'T')
	if altCodon != "TGT" {
		t.Errorf("Expected TGT after mutation, got %s", altCodon)
	}

	altAA := TranslateCodon(altCodon)
	if altAA != 'C' {
		t.Errorf("Alternate codon TGT should translate to C (Cys), got %c", altAA)
	}

	// The amino acid change notation
	aaChange := string(refAA) + "12" + string(altAA)
	if aaChange != "G12C" {
		t.Errorf("Expected amino acid change G12C, got %s", aaChange)
	}
}
