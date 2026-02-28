package annotate

import (
	"fmt"
	"regexp"
	"strconv"
	"strings"
)

// VariantSpecType identifies the kind of variant specification.
type VariantSpecType int

const (
	SpecGenomic VariantSpecType = iota
	SpecProtein
	SpecHGVSc
)

// VariantSpec holds a parsed variant specification.
type VariantSpec struct {
	Type VariantSpecType
	// Genomic fields
	Chrom string
	Pos   int64
	Ref   string
	Alt   string
	// Protein fields
	GeneName string
	RefAA    byte  // single-letter
	Position int64 // protein position
	AltAA    byte  // single-letter
	// HGVSc fields
	TranscriptID string // transcript ID or gene name
	CDSChange    string // e.g. "35G>T"
}

// AminoAcidThreeToSingle maps three-letter amino acid codes to single-letter.
var AminoAcidThreeToSingle map[string]byte

func init() {
	AminoAcidThreeToSingle = make(map[string]byte, len(AminoAcidSingleToThree))
	for single, three := range AminoAcidSingleToThree {
		AminoAcidThreeToSingle[three] = single
	}
}

// Regexes for variant spec parsing.
var (
	// Genomic: chr12:25245350:C:A  or  12-25245350-C-A  or  chr12:25245350:C>A
	reGenomic = regexp.MustCompile(`^(chr)?(\w+)[:\-](\d+)[:\-]([ACGTNacgtn]+)[>:\-/]([ACGTNacgtn]+)$`)
	// HGVSc: KRAS c.35G>T  or  ENST00000311936:c.35G>T
	reHGVSc = regexp.MustCompile(`^(\S+)\s+c\.(.+)$`)
	// HGVSc with colon separator: ENST00000311936:c.35G>T
	reHGVScColon = regexp.MustCompile(`^(\S+):c\.(.+)$`)
	// Protein with three-letter codes: KRAS p.Gly12Cys
	reProteinThree = regexp.MustCompile(`^(\w+)\s+p?\.?([A-Z][a-z]{2}|\*)(\d+)([A-Z][a-z]{2}|\*)$`)
	// Protein with single-letter codes: KRAS G12C  or  KRAS p.G12C
	reProteinSingle = regexp.MustCompile(`^(\w+)\s+p?\.?([ACDEFGHIKLMNPQRSTVWY*])(\d+)([ACDEFGHIKLMNPQRSTVWY*])$`)
)

// ParseVariantSpec parses a variant specification string into a VariantSpec.
// It tries genomic, then HGVSc, then protein format.
func ParseVariantSpec(input string) (*VariantSpec, error) {
	input = strings.TrimSpace(input)
	if input == "" {
		return nil, fmt.Errorf("empty variant specification")
	}

	// Try genomic format
	if spec, ok := parseGenomic(input); ok {
		return spec, nil
	}

	// Try HGVSc with colon separator (must be before generic HGVSc)
	if spec, ok := parseHGVScColon(input); ok {
		return spec, nil
	}

	// Try HGVSc with space separator
	if spec, ok := parseHGVSc(input); ok {
		return spec, nil
	}

	// Try protein (three-letter codes first, then single-letter)
	if spec, ok := parseProteinThree(input); ok {
		return spec, nil
	}
	if spec, ok := parseProteinSingle(input); ok {
		return spec, nil
	}

	return nil, fmt.Errorf("cannot parse variant specification %q (expected genomic, protein, or HGVSc format)", input)
}

func parseGenomic(input string) (*VariantSpec, bool) {
	m := reGenomic.FindStringSubmatch(input)
	if m == nil {
		return nil, false
	}
	pos, err := strconv.ParseInt(m[3], 10, 64)
	if err != nil {
		return nil, false
	}
	return &VariantSpec{
		Type:  SpecGenomic,
		Chrom: m[2],
		Pos:   pos,
		Ref:   strings.ToUpper(m[4]),
		Alt:   strings.ToUpper(m[5]),
	}, true
}

func parseHGVSc(input string) (*VariantSpec, bool) {
	m := reHGVSc.FindStringSubmatch(input)
	if m == nil {
		return nil, false
	}
	// Verify the CDS change looks valid (contains digits)
	if !strings.ContainsAny(m[2], "0123456789") {
		return nil, false
	}
	return &VariantSpec{
		Type:         SpecHGVSc,
		TranscriptID: m[1],
		CDSChange:    m[2],
	}, true
}

func parseHGVScColon(input string) (*VariantSpec, bool) {
	m := reHGVScColon.FindStringSubmatch(input)
	if m == nil {
		return nil, false
	}
	if !strings.ContainsAny(m[2], "0123456789") {
		return nil, false
	}
	return &VariantSpec{
		Type:         SpecHGVSc,
		TranscriptID: m[1],
		CDSChange:    m[2],
	}, true
}

func parseProteinThree(input string) (*VariantSpec, bool) {
	m := reProteinThree.FindStringSubmatch(input)
	if m == nil {
		return nil, false
	}
	pos, err := strconv.ParseInt(m[3], 10, 64)
	if err != nil {
		return nil, false
	}
	refAA := threeLetterToSingle(m[2])
	altAA := threeLetterToSingle(m[4])
	if refAA == 0 || altAA == 0 {
		return nil, false
	}
	return &VariantSpec{
		Type:     SpecProtein,
		GeneName: m[1],
		RefAA:    refAA,
		Position: pos,
		AltAA:    altAA,
	}, true
}

func parseProteinSingle(input string) (*VariantSpec, bool) {
	m := reProteinSingle.FindStringSubmatch(input)
	if m == nil {
		return nil, false
	}
	pos, err := strconv.ParseInt(m[3], 10, 64)
	if err != nil {
		return nil, false
	}
	return &VariantSpec{
		Type:     SpecProtein,
		GeneName: m[1],
		RefAA:    m[2][0],
		Position: pos,
		AltAA:    m[4][0],
	}, true
}

func threeLetterToSingle(code string) byte {
	if code == "*" {
		return '*'
	}
	if aa, ok := AminoAcidThreeToSingle[code]; ok {
		return aa
	}
	return 0
}
