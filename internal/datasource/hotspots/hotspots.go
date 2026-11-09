// Package hotspots provides cancer hotspot mutation lookups.
// Data comes from the Cancer Hotspots database (cancerhotspots.org).
package hotspots

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
)

// Hotspot represents a single hotspot entry at a protein position.
type Hotspot struct {
	Position     int64   // first amino acid position the hotspot covers
	PositionEnd  int64   // last position covered; equal to Position for a single residue
	TranscriptID string  // Ensembl transcript (unversioned, e.g. "ENST00000311936")
	Type         string  // "single residue", "in-frame indel", "3d", "splice site"
	QValue       float64 // statistical significance (q-value)

	// The rest are reported by the /cancer_hotspots endpoints, which answer
	// from this table directly rather than through a variant annotation.
	HugoSymbol          string
	Residue             string // e.g. "Q61"
	TranscriptIDVersion string
	Version             string // the hotspots release this row came from
	TumorCount          int
	// The per-class counts are absent on 3d rows, where the file leaves the
	// columns empty. genome-nexus reports them as null rather than zero, so
	// they are pointers to keep the distinction.
	MissenseCount   *int
	TruncatingCount *int
	InframeCount    *int
	SpliceCount     *int
}

// Store holds hotspot data keyed by transcript ID.
// Positions are only meaningful for the specific transcript they were defined on,
// so lookups require matching the transcript, not just the gene.
type Store struct {
	data  map[string][]Hotspot // transcript ID → every row, sorted by position
	byPos map[string][]Hotspot // transcript ID → one single-residue row per position
}

// Load parses a hotspots TSV file (hotspots_v2_and_3d.txt format).
// Expected columns: hugo_symbol, amino_acid_position, type, q_value (among others).
func Load(path string) (*Store, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open hotspots file: %w", err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	// Read header
	if !scanner.Scan() {
		return nil, fmt.Errorf("empty hotspots file")
	}
	header := strings.Split(scanner.Text(), "\t")
	colIdx := indexColumns(header,
		"hugo_symbol", "amino_acid_position", "type", "q_value", "transcript_id",
		"residue", "transcript_id_version", "version",
		"tumor_count", "missense_count", "trunc_count", "inframe_count", "splice_count")
	if colIdx["amino_acid_position"] < 0 || colIdx["transcript_id"] < 0 {
		return nil, fmt.Errorf("missing required columns: amino_acid_position, transcript_id")
	}

	data := make(map[string][]Hotspot)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")

		txIdx := colIdx["transcript_id"]
		posIdx := colIdx["amino_acid_position"]
		if len(fields) <= txIdx || len(fields) <= posIdx {
			continue
		}

		txID := fields[txIdx]
		if txID == "" {
			continue // skip rows without a transcript
		}
		start, end, ok := parsePositionRange(fields[posIdx])
		if !ok {
			continue // skip rows with an unparseable position
		}

		h := Hotspot{Position: start, PositionEnd: end, TranscriptID: txID}

		if idx := colIdx["type"]; idx >= 0 && idx < len(fields) {
			h.Type = fields[idx]
		}
		if idx := colIdx["q_value"]; idx >= 0 && idx < len(fields) && fields[idx] != "" {
			if v, err := strconv.ParseFloat(fields[idx], 64); err == nil {
				h.QValue = v
			}
		}

		text := func(name string) string {
			if idx := colIdx[name]; idx >= 0 && idx < len(fields) {
				return fields[idx]
			}
			return ""
		}
		count := func(name string) *int {
			raw := strings.TrimSpace(text(name))
			if raw == "" {
				return nil
			}
			// The counts are written as floats ("897.0"), so Atoi will not do.
			f, err := strconv.ParseFloat(raw, 64)
			if err != nil {
				return nil
			}
			n := int(f)
			return &n
		}
		h.HugoSymbol = text("hugo_symbol")
		h.Residue = text("residue")
		h.TranscriptIDVersion = text("transcript_id_version")
		h.Version = text("version")
		if n := count("tumor_count"); n != nil {
			h.TumorCount = *n
		}
		h.MissenseCount = count("missense_count")
		h.TruncatingCount = count("trunc_count")
		h.InframeCount = count("inframe_count")
		h.SpliceCount = count("splice_count")

		data[txID] = append(data[txID], h)
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("read hotspots file: %w", err)
	}

	byPos := make(map[string][]Hotspot, len(data))
	for tx, spots := range data {
		sort.SliceStable(spots, func(i, j int) bool {
			return spots[i].Position < spots[j].Position
		})
		data[tx] = dedupRows(spots)
		byPos[tx] = collapseByPosition(data[tx])
	}

	return &Store{data: data, byPos: byPos}, nil
}

// LookupAll returns every hotspot on a transcript, with or without a version
// suffix on the id. The /cancer_hotspots endpoints report them all, where
// annotation only ever asks about one position.
func (s *Store) LookupAll(transcriptID string) []Hotspot {
	if s == nil {
		return nil
	}
	return s.data[unversioned(transcriptID)]
}

// unversioned drops the ".N" suffix an Ensembl id may carry. The table is
// keyed on the bare id, so a caller supplying a version still matches.
func unversioned(id string) string {
	if i := strings.IndexByte(id, '.'); i >= 0 {
		return id[:i]
	}
	return id
}

// Lookup checks if a transcript+position is a known hotspot.
// The transcriptID should be unversioned (e.g. "ENST00000311936").
func (s *Store) Lookup(transcriptID string, proteinPosition int64) (Hotspot, bool) {
	spots := s.byPos[transcriptID]
	if len(spots) == 0 {
		return Hotspot{}, false
	}

	// Binary search for position
	i := sort.Search(len(spots), func(i int) bool {
		return spots[i].Position >= proteinPosition
	})
	if i < len(spots) && spots[i].Position == proteinPosition {
		return spots[i], true
	}
	return Hotspot{}, false
}

// TranscriptCount returns the number of transcripts with hotspots.
func (s *Store) TranscriptCount() int {
	return len(s.data)
}

// HotspotCount returns the total number of hotspot entries.
func (s *Store) HotspotCount() int {
	var n int
	for _, spots := range s.data {
		n += len(spots)
	}
	return n
}

// typePriority returns a priority for hotspot types (lower = preferred).
// "single residue" is direct evidence, preferred over structural prediction ("3d").
func typePriority(t string) int {
	switch t {
	case "single residue":
		return 0
	case "in-frame indel":
		return 1
	case "splice site":
		return 2
	case "3d":
		return 3
	default:
		return 4
	}
}

// parsePositionRange reads an amino_acid_position, which is a single number
// for most rows but a span ("486-494") on in-frame indel hotspots.
func parsePositionRange(s string) (start, end int64, ok bool) {
	s = strings.TrimSpace(s)
	lo, hi, isRange := strings.Cut(s, "-")
	start, err := strconv.ParseInt(lo, 10, 64)
	if err != nil {
		return 0, 0, false
	}
	if !isRange {
		return start, start, true
	}
	end, err = strconv.ParseInt(hi, 10, 64)
	if err != nil {
		return 0, 0, false
	}
	return start, end, true
}

// dedupRows collapses rows that the importer would have merged. It loads this
// same file into MongoDB with
//
//	--mode upsert --upsertFields hugo_symbol,residue,type,tumor_count
//
// so rows agreeing on those four are one document, the last one read winning.
// 3d hotspots repeat per structural cluster and are the usual case.
func dedupRows(spots []Hotspot) []Hotspot {
	type key struct {
		hugo, residue, hsType string
		tumorCount            int
	}
	at := make(map[key]int, len(spots))
	result := spots[:0]
	for _, h := range spots {
		k := key{h.HugoSymbol, h.Residue, h.Type, h.TumorCount}
		if i, seen := at[k]; seen {
			result[i] = h // last row wins, as the upsert does
			continue
		}
		at[k] = len(result)
		result = append(result, h)
	}
	return result
}

// collapseByPosition is the view annotation uses: at most one hotspot per
// residue, preferring the lowest q-value and then the strongest evidence.
// Hotspots spanning a range are left out, as they were when the loader could
// not parse their position at all.
func collapseByPosition(spots []Hotspot) []Hotspot {
	result := make([]Hotspot, 0, len(spots))
	for _, h := range spots {
		if h.PositionEnd != h.Position {
			continue
		}
		if n := len(result); n > 0 && result[n-1].Position == h.Position {
			cur := &result[n-1]
			if h.QValue < cur.QValue ||
				(h.QValue == cur.QValue && typePriority(h.Type) < typePriority(cur.Type)) {
				*cur = h
			}
			continue
		}
		result = append(result, h)
	}
	return result
}

// indexColumns maps column names to their indices in the header.
func indexColumns(header []string, names ...string) map[string]int {
	idx := make(map[string]int, len(names))
	for _, name := range names {
		idx[name] = -1
	}
	for i, col := range header {
		for _, name := range names {
			if col == name {
				idx[name] = i
			}
		}
	}
	return idx
}

// FormatQValue formats a q-value for output.
func FormatQValue(v float64) string {
	if v == 0 {
		return "0"
	}
	if math.Abs(v) < 0.001 {
		return strconv.FormatFloat(v, 'e', 3, 64)
	}
	return strconv.FormatFloat(v, 'f', -1, 64)
}
