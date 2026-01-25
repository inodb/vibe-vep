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

// Complement map for reverse complementing DNA sequences.
var complementMap = map[byte]byte{
	'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
	'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
	'N': 'N', 'n': 'n',
}

// TranslateCodon translates a DNA codon to its amino acid.
// Returns 'X' for unknown codons and '*' for stop codons.
func TranslateCodon(codon string) byte {
	if len(codon) != 3 {
		return 'X'
	}
	codon = strings.ToUpper(codon)
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
	return strings.ToUpper(codon) == "ATG"
}

// ReverseComplement returns the reverse complement of a DNA sequence.
// For strand handling: reverse strand genes need their coding sequence
// reverse complemented to get the template strand.
func ReverseComplement(seq string) string {
	n := len(seq)
	result := make([]byte, n)
	for i := 0; i < n; i++ {
		base := seq[n-1-i]
		if comp, ok := complementMap[base]; ok {
			result[i] = comp
		} else {
			result[i] = 'N' // Unknown base
		}
	}
	return string(result)
}

// Complement returns the complement of a single base.
func Complement(base byte) byte {
	if comp, ok := complementMap[base]; ok {
		return comp
	}
	return 'N'
}

// TranslateSequence translates a DNA sequence to amino acids.
// Sequence length must be divisible by 3.
func TranslateSequence(seq string) string {
	seq = strings.ToUpper(seq)
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
	codonBytes := []byte(codon)
	codonBytes[positionInCodon] = newBase
	return string(codonBytes)
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
