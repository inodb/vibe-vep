// Package annotate provides variant effect prediction functionality.
package annotate

import "strings"

// Standard genetic code: DNA codon to amino acid (single letter).
var codonTable = map[string]byte{
	"TTT": 'F', "TTC": 'F', "TTA": 'L', "TTG": 'L',
	"TCT": 'S', "TCC": 'S', "TCA": 'S', "TCG": 'S',
	"TAT": 'Y', "TAC": 'Y', "TAA": '*', "TAG": '*',
	"TGT": 'C', "TGC": 'C', "TGA": '*', "TGG": 'W',

	"CTT": 'L', "CTC": 'L', "CTA": 'L', "CTG": 'L',
	"CCT": 'P', "CCC": 'P', "CCA": 'P', "CCG": 'P',
	"CAT": 'H', "CAC": 'H', "CAA": 'Q', "CAG": 'Q',
	"CGT": 'R', "CGC": 'R', "CGA": 'R', "CGG": 'R',

	"ATT": 'I', "ATC": 'I', "ATA": 'I', "ATG": 'M',
	"ACT": 'T', "ACC": 'T', "ACA": 'T', "ACG": 'T',
	"AAT": 'N', "AAC": 'N', "AAA": 'K', "AAG": 'K',
	"AGT": 'S', "AGC": 'S', "AGA": 'R', "AGG": 'R',

	"GTT": 'V', "GTC": 'V', "GTA": 'V', "GTG": 'V',
	"GCT": 'A', "GCC": 'A', "GCA": 'A', "GCG": 'A',
	"GAT": 'D', "GAC": 'D', "GAA": 'E', "GAG": 'E',
	"GGT": 'G', "GGC": 'G', "GGA": 'G', "GGG": 'G',
}

// TranslateCodon translates a DNA codon to its amino acid.
// Returns 'X' for unknown codons and '*' for stop codons.
// CDS data is already uppercase, so no ToUpper conversion is needed.
func TranslateCodon(codon string) byte {
	if len(codon) != 3 {
		return 'X'
	}
	if aa, ok := codonTable[codon]; ok {
		return aa
	}
	return 'X'
}

// IsStopCodon returns true if the codon is a stop codon (TAA, TAG, TGA).
func IsStopCodon(codon string) bool {
	return TranslateCodon(codon) == '*'
}

// IsStartCodon returns true if the codon is the start codon (ATG).
func IsStartCodon(codon string) bool {
	return codon == "ATG"
}

// ReverseComplement returns the reverse complement of a DNA sequence.
// For strand handling: reverse strand genes need their coding sequence
// reverse complemented to get the template strand.
func ReverseComplement(seq string) string {
	n := len(seq)
	// Stack-allocate for typical variant ref/alt lengths (≤64 bases).
	var buf [64]byte
	var result []byte
	if n <= len(buf) {
		result = buf[:n]
	} else {
		result = make([]byte, n)
	}
	for i := 0; i < n; i++ {
		result[i] = Complement(seq[n-1-i])
	}
	return string(result)
}

// Complement returns the complement of a single base.
func Complement(base byte) byte {
	switch base {
	case 'A':
		return 'T'
	case 'T':
		return 'A'
	case 'G':
		return 'C'
	case 'C':
		return 'G'
	case 'a':
		return 't'
	case 't':
		return 'a'
	case 'g':
		return 'c'
	case 'c':
		return 'g'
	default:
		return 'N'
	}
}

// TranslateSequence translates a DNA sequence to amino acids.
// Sequence length must be divisible by 3.
func TranslateSequence(seq string) string {
	n := len(seq)
	if n%3 != 0 {
		// Truncate to complete codons
		n = (n / 3) * 3
	}

	var result strings.Builder
	result.Grow(n / 3)

	for i := 0; i < n; i += 3 {
		codon := seq[i : i+3]
		aa := TranslateCodon(codon)
		result.WriteByte(aa)
	}

	return result.String()
}

// GetCodon extracts a codon from a CDS sequence at a given codon number.
// Codon numbers are 1-based (codon 1 is positions 1-3).
func GetCodon(cdsSequence string, codonNumber int64) string {
	if codonNumber < 1 {
		return ""
	}
	startIdx := (codonNumber - 1) * 3
	endIdx := startIdx + 3
	if endIdx > int64(len(cdsSequence)) {
		return ""
	}
	return cdsSequence[startIdx:endIdx]
}

// MutateCodon applies a mutation to a codon at a specific position.
// positionInCodon is 0, 1, or 2 (first, second, or third base).
func MutateCodon(codon string, positionInCodon int, newBase byte) string {
	if len(codon) != 3 || positionInCodon < 0 || positionInCodon > 2 {
		return codon
	}
	var buf [3]byte
	copy(buf[:], codon)
	buf[positionInCodon] = newBase
	return string(buf[:])
}

// AminoAcidSingleToThree converts single letter amino acid to three letter code.
var AminoAcidSingleToThree = map[byte]string{
	'A': "Ala", 'C': "Cys", 'D': "Asp", 'E': "Glu",
	'F': "Phe", 'G': "Gly", 'H': "His", 'I': "Ile",
	'K': "Lys", 'L': "Leu", 'M': "Met", 'N': "Asn",
	'P': "Pro", 'Q': "Gln", 'R': "Arg", 'S': "Ser",
	'T': "Thr", 'V': "Val", 'W': "Trp", 'Y': "Tyr",
	'*': "Ter", 'X': "Xaa",
}
