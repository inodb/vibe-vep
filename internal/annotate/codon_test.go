package annotate

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestTranslateCodon(t *testing.T) {
	tests := []struct {
		name  string
		codon string
		want  byte
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

		// Lowercase not supported (CDS data is always uppercase)
		{"lowercase atg", "atg", 'X'},
		{"mixed case AtG", "AtG", 'X'},

		// Invalid codons
		{"too short", "AT", 'X'},
		{"too long", "ATGG", 'X'},
		{"invalid bases", "XYZ", 'X'},
		{"empty", "", 'X'},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := TranslateCodon(tt.codon)
			assert.Equal(t, tt.want, got, "TranslateCodon(%q)", tt.codon)
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
			assert.Equal(t, tt.want, got, "ReverseComplement(%q)", tt.seq)
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
		{"taa", false}, // lowercase not supported (CDS data is always uppercase)
	}

	for _, tt := range tests {
		t.Run(tt.codon, func(t *testing.T) {
			got := IsStopCodon(tt.codon)
			assert.Equal(t, tt.want, got)
		})
	}
}

func TestIsStartCodon(t *testing.T) {
	tests := []struct {
		codon string
		want  bool
	}{
		{"ATG", true},
		{"atg", false}, // lowercase not supported
		{"TAA", false},
		{"GGT", false},
	}

	for _, tt := range tests {
		t.Run(tt.codon, func(t *testing.T) {
			got := IsStartCodon(tt.codon)
			assert.Equal(t, tt.want, got)
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
			assert.Equal(t, tt.want, got)
		})
	}
}

func TestMutateCodon(t *testing.T) {
	tests := []struct {
		name            string
		codon           string
		positionInCodon int
		newBase         byte
		want            string
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
			assert.Equal(t, tt.want, got)
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
			assert.Equal(t, tt.want, got)
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
	assert.Equal(t, byte('G'), refAA, "Reference codon GGT should translate to G (Gly)")

	// Mutate first position: G -> T gives TGT
	altCodon := MutateCodon(refCodon, 0, 'T')
	assert.Equal(t, "TGT", altCodon, "Expected TGT after mutation")

	altAA := TranslateCodon(altCodon)
	assert.Equal(t, byte('C'), altAA, "Alternate codon TGT should translate to C (Cys)")

	// The amino acid change notation
	aaChange := string(refAA) + "12" + string(altAA)
	assert.Equal(t, "G12C", aaChange, "Expected amino acid change G12C")
}

func BenchmarkTranslateCodon(b *testing.B) {
	codons := []string{"ATG", "GGT", "TGT", "TAA", "TAG", "TGA", "CCC", "AAA"}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for _, c := range codons {
			TranslateCodon(c)
		}
	}
}

func BenchmarkComplement(b *testing.B) {
	bases := []byte{'A', 'T', 'G', 'C', 'N'}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for _, base := range bases {
			Complement(base)
		}
	}
}

func BenchmarkReverseComplement(b *testing.B) {
	seqs := []string{"ATGCGTACGTAGCTAG", "GGT", "CCCAAATTTGGG"}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for _, s := range seqs {
			ReverseComplement(s)
		}
	}
}

// TestAllocRegression_TranslateCodon verifies zero allocations for codon translation.
func TestAllocRegression_TranslateCodon(t *testing.T) {
	allocs := testing.AllocsPerRun(100, func() {
		TranslateCodon("ATG")
		TranslateCodon("TAA")
		TranslateCodon("GGT")
	})
	if allocs > 0 {
		t.Errorf("TranslateCodon allocates: got %.0f, want 0", allocs)
	}
}

// TestAllocRegression_Complement verifies zero allocations for complement.
func TestAllocRegression_Complement(t *testing.T) {
	allocs := testing.AllocsPerRun(100, func() {
		Complement('A')
		Complement('T')
		Complement('G')
		Complement('C')
	})
	if allocs > 0 {
		t.Errorf("Complement allocates: got %.0f, want 0", allocs)
	}
}

// TestAllocRegression_ReverseComplement verifies minimal allocations for small sequences.
func TestAllocRegression_ReverseComplement(t *testing.T) {
	// Typical codon (3 bases) - should be 1 alloc (string return)
	allocs := testing.AllocsPerRun(100, func() {
		ReverseComplement("ATG")
	})
	if allocs > 1 {
		t.Errorf("ReverseComplement(3bp) allocs: %.0f, want <= 1", allocs)
	}

	// Short variant ref/alt (10 bases)
	allocs = testing.AllocsPerRun(100, func() {
		ReverseComplement("ATGCATGCAT")
	})
	if allocs > 1 {
		t.Errorf("ReverseComplement(10bp) allocs: %.0f, want <= 1", allocs)
	}
}

// TestAllocRegression_MutateCodon verifies minimal allocations.
func TestAllocRegression_MutateCodon(t *testing.T) {
	allocs := testing.AllocsPerRun(100, func() {
		MutateCodon("GGT", 2, 'A')
	})
	if allocs > 1 {
		t.Errorf("MutateCodon allocs: %.0f, want <= 1", allocs)
	}
}
