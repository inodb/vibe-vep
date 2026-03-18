// Package dbnsfp provides dbNSFP TSV parsing helpers.
// The actual store is in the genomicindex package (unified SQLite).
package dbnsfp

import (
	"strconv"
	"strings"
)

// Column names in dbNSFP TSV files.
const (
	ColChrom     = "#chr"
	ColPos       = "pos(1-based)"
	ColRef       = "ref"
	ColAlt       = "alt"
	ColSIFTScore = "SIFT_score"
	ColSIFTPred  = "SIFT_pred"
	ColPP2Score  = "Polyphen2_HDIV_score"
	ColPP2Pred   = "Polyphen2_HDIV_pred"
)

// Entry represents a single dbNSFP variant annotation.
type Entry struct {
	Pos       int64
	Ref       string
	Alt       string
	SIFTScore float32
	SIFTPred  string
	PP2Score  float32
	PP2Pred   string
}

// IndexColumns maps column names to their indices in a TSV header.
// Returns a map; missing columns have value -1.
func IndexColumns(header []string, names ...string) map[string]int {
	m := make(map[string]int, len(names))
	for _, name := range names {
		m[name] = -1
	}
	for i, col := range header {
		if _, ok := m[col]; ok {
			m[col] = i
		}
	}
	return m
}

// ParseLine parses a single dbNSFP TSV data line using pre-indexed column positions.
// altIndex is the 0-based index of the desired alt allele in multi-allelic rows.
// Returns entries (one per alt allele), normalized chromosome, and whether parsing succeeded.
func ParseLine(fields []string, col map[string]int) ([]Entry, string, bool) {
	chromIdx := col[ColChrom]
	posIdx := col[ColPos]
	refIdx := col[ColRef]
	altIdx := col[ColAlt]

	if chromIdx < 0 || posIdx < 0 || refIdx < 0 || altIdx < 0 {
		return nil, "", false
	}
	if chromIdx >= len(fields) || posIdx >= len(fields) || refIdx >= len(fields) || altIdx >= len(fields) {
		return nil, "", false
	}

	chrom := NormalizeChrom(fields[chromIdx])
	pos, err := strconv.ParseInt(fields[posIdx], 10, 64)
	if err != nil {
		return nil, "", false
	}

	ref := fields[refIdx]
	alts := strings.Split(fields[altIdx], ";")

	siftScores := getMultiField(fields, col[ColSIFTScore])
	siftPreds := getMultiField(fields, col[ColSIFTPred])
	pp2Scores := getMultiField(fields, col[ColPP2Score])
	pp2Preds := getMultiField(fields, col[ColPP2Pred])

	var entries []Entry
	for i, alt := range alts {
		if alt == "." || alt == "" {
			continue
		}
		e := Entry{
			Pos: pos,
			Ref: ref,
			Alt: alt,
		}

		if s := getIndex(siftScores, i); s != "" && s != "." {
			if v, err := strconv.ParseFloat(s, 32); err == nil {
				e.SIFTScore = float32(v)
			}
		}
		if p := getIndex(siftPreds, i); p != "" && p != "." {
			e.SIFTPred = ExpandPrediction(p)
		}
		if s := getIndex(pp2Scores, i); s != "" && s != "." {
			if v, err := strconv.ParseFloat(s, 32); err == nil {
				e.PP2Score = float32(v)
			}
		}
		if p := getIndex(pp2Preds, i); p != "" && p != "." {
			e.PP2Pred = ExpandPrediction(p)
		}

		// Skip entries with no data.
		if e.SIFTScore == 0 && e.SIFTPred == "" && e.PP2Score == 0 && e.PP2Pred == "" {
			continue
		}

		entries = append(entries, e)
	}

	if len(entries) == 0 {
		return nil, "", false
	}
	return entries, chrom, true
}

// ExpandPrediction maps single-letter prediction codes to full words.
func ExpandPrediction(code string) string {
	switch code {
	case "D":
		return "deleterious"
	case "T":
		return "tolerated"
	case "P":
		return "possibly_damaging"
	case "B":
		return "benign"
	default:
		return code
	}
}

// NormalizeChrom removes "chr" prefix for consistent lookups.
func NormalizeChrom(chrom string) string {
	return strings.TrimPrefix(chrom, "chr")
}

// getMultiField splits a semicolon-separated field value.
func getMultiField(fields []string, idx int) []string {
	if idx < 0 || idx >= len(fields) {
		return nil
	}
	return strings.Split(fields[idx], ";")
}

// getIndex returns the i-th element of a slice, or "" if out of range.
// If the slice has only one element, returns that (dbNSFP sometimes uses a single
// value for all alt alleles).
func getIndex(s []string, i int) string {
	if len(s) == 0 {
		return ""
	}
	if i < len(s) {
		return s[i]
	}
	if len(s) == 1 {
		return s[0]
	}
	return ""
}
