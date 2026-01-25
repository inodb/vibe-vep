// Package cache provides VEP cache loading functionality.
package cache

import (
	"fmt"

	"github.com/Sereal/Sereal/Go/sereal"
)

// Sereal magic bytes for format detection.
// Standard magic is "=srl" (0x3D 0x73 0x72 0x6C)
// High-bit variant is "=\xF3rl" (0x3D 0xF3 0x72 0x6C)
var (
	serealMagicStandard = []byte{0x3D, 0x73, 0x72, 0x6C} // =srl
	serealMagicHighBit  = []byte{0x3D, 0xF3, 0x72, 0x6C} // =\xF3rl
)

// IsSereal checks if the data starts with Sereal magic bytes.
func IsSereal(data []byte) bool {
	if len(data) < 4 {
		return false
	}
	return matchMagic(data[:4], serealMagicStandard) ||
		matchMagic(data[:4], serealMagicHighBit)
}

// matchMagic compares two byte slices for equality.
func matchMagic(a, b []byte) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// DecodeSereal decodes Sereal-encoded VEP cache data into transcripts.
func DecodeSereal(data []byte, chrom string) ([]*Transcript, error) {
	var raw interface{}
	if err := sereal.Unmarshal(data, &raw); err != nil {
		return nil, fmt.Errorf("sereal unmarshal: %w", err)
	}

	// VEP cache files contain an array of transcript hashes
	rawArray, ok := raw.([]interface{})
	if !ok {
		// Try to unwrap from map if needed
		if rawMap, isMap := raw.(map[string]interface{}); isMap {
			if items, hasItems := rawMap["transcripts"]; hasItems {
				rawArray, ok = items.([]interface{})
			}
		}
		if !ok {
			return nil, fmt.Errorf("expected array of transcripts, got %T", raw)
		}
	}

	transcripts := make([]*Transcript, 0, len(rawArray))
	for i, item := range rawArray {
		t, err := convertTranscript(item, chrom)
		if err != nil {
			return nil, fmt.Errorf("convert transcript %d: %w", i, err)
		}
		if t != nil {
			transcripts = append(transcripts, t)
		}
	}

	return transcripts, nil
}

// convertTranscript converts a raw Perl hash to a Go Transcript struct.
func convertTranscript(raw interface{}, chrom string) (*Transcript, error) {
	// Handle Perl blessed objects
	m := unwrapPerlObject(raw)
	if m == nil {
		return nil, fmt.Errorf("expected hash, got %T", raw)
	}

	t := &Transcript{
		Chrom: chrom,
	}

	// Extract basic fields
	t.ID = getString(m, "stable_id")
	t.GeneID = getString(m, "_gene_stable_id")
	t.GeneName = getString(m, "_gene_symbol")
	t.Start = getInt64(m, "start")
	t.End = getInt64(m, "end")
	t.Strand = int8(getInt64(m, "strand"))
	t.Biotype = getString(m, "biotype")
	t.IsCanonical = getBool(m, "is_canonical")

	// Extract CDS boundaries
	t.CDSStart = getInt64(m, "coding_region_start")
	t.CDSEnd = getInt64(m, "coding_region_end")

	// Extract sequences
	t.CDSSequence = getString(m, "translateable_seq")

	// Protein sequence may be nested in translation object
	if translation := getMap(m, "translation"); translation != nil {
		t.ProteinSequence = getString(translation, "seq")
	}

	// Extract exons
	if exonArray := getArray(m, "_trans_exon_array"); exonArray != nil {
		exons, err := convertExons(exonArray, t.CDSStart, t.CDSEnd)
		if err != nil {
			return nil, fmt.Errorf("convert exons: %w", err)
		}
		t.Exons = exons
	}

	// Skip invalid transcripts
	if t.ID == "" {
		return nil, nil
	}

	return t, nil
}

// convertExons converts a raw Perl array to Go Exon structs.
func convertExons(raw []interface{}, cdsStart, cdsEnd int64) ([]Exon, error) {
	exons := make([]Exon, 0, len(raw))

	for i, item := range raw {
		m := unwrapPerlObject(item)
		if m == nil {
			continue
		}

		exon := Exon{
			Number: i + 1,
			Start:  getInt64(m, "start"),
			End:    getInt64(m, "end"),
			Frame:  -1, // Default to non-coding
		}

		// Calculate CDS portion of exon
		if cdsStart > 0 && cdsEnd > 0 {
			// Check if exon overlaps CDS
			if exon.End >= cdsStart && exon.Start <= cdsEnd {
				exon.CDSStart = max(exon.Start, cdsStart)
				exon.CDSEnd = min(exon.End, cdsEnd)

				// Get phase/frame if available
				if phase := getInt64(m, "phase"); phase >= 0 {
					exon.Frame = int(phase)
				}
			}
		}

		exons = append(exons, exon)
	}

	return exons, nil
}

// unwrapPerlObject extracts the underlying map from a Perl object or returns
// the map directly if it's already a map.
func unwrapPerlObject(raw interface{}) map[string]interface{} {
	// Direct map
	if m, ok := raw.(map[string]interface{}); ok {
		return m
	}

	// Sereal PerlObject (blessed reference)
	if obj, ok := raw.(sereal.PerlObject); ok {
		if m, ok := obj.Reference.(map[string]interface{}); ok {
			return m
		}
	}

	return nil
}

// getString extracts a string value from a map.
func getString(m map[string]interface{}, key string) string {
	if v, ok := m[key]; ok {
		switch s := v.(type) {
		case string:
			return s
		case []byte:
			return string(s)
		}
	}
	return ""
}

// getInt64 extracts an integer value from a map.
func getInt64(m map[string]interface{}, key string) int64 {
	if v, ok := m[key]; ok {
		switch n := v.(type) {
		case int64:
			return n
		case int:
			return int64(n)
		case float64:
			return int64(n)
		case uint64:
			return int64(n)
		}
	}
	return 0
}

// getBool extracts a boolean value from a map.
// Handles Perl's loose boolean typing (0/1, true/false, etc.).
func getBool(m map[string]interface{}, key string) bool {
	if v, ok := m[key]; ok {
		switch b := v.(type) {
		case bool:
			return b
		case int64:
			return b != 0
		case int:
			return b != 0
		case float64:
			return b != 0
		case string:
			return b == "1" || b == "true"
		}
	}
	return false
}

// getMap extracts a nested map from a map.
func getMap(m map[string]interface{}, key string) map[string]interface{} {
	if v, ok := m[key]; ok {
		return unwrapPerlObject(v)
	}
	return nil
}

// getArray extracts an array from a map.
func getArray(m map[string]interface{}, key string) []interface{} {
	if v, ok := m[key]; ok {
		if arr, ok := v.([]interface{}); ok {
			return arr
		}
	}
	return nil
}
