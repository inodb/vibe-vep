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
	Position     int64   // amino acid position
	TranscriptID string  // Ensembl transcript (unversioned, e.g. "ENST00000311936")
	Type         string  // "single residue", "in-frame indel", "3d", "splice"
	QValue       float64 // statistical significance (q-value)

	// The rest are reported by the /cancer_hotspots endpoints, which answer
	// from this table directly rather than through a variant annotation.
	HugoSymbol          string
	Residue             string // e.g. "Q61"
	TranscriptIDVersion string
	Version             string // the hotspots release this row came from
	TumorCount          int
	MissenseCount       int
	TruncatingCount     int
	InframeCount        int
	SpliceCount         int
}

// Store holds hotspot data keyed by transcript ID.
// Positions are only meaningful for the specific transcript they were defined on,
// so lookups require matching the transcript, not just the gene.
type Store struct {
	data map[string][]Hotspot // transcript ID → sorted hotspots
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
		posStr := fields[posIdx]
		pos, err := strconv.ParseInt(posStr, 10, 64)
		if err != nil {
			continue // skip rows with non-numeric positions
		}

		h := Hotspot{Position: pos, TranscriptID: txID}

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
		count := func(name string) int {
			n, _ := strconv.Atoi(text(name))
			return n
		}
		h.HugoSymbol = text("hugo_symbol")
		h.Residue = text("residue")
		h.TranscriptIDVersion = text("transcript_id_version")
		h.Version = text("version")
		h.TumorCount = count("tumor_count")
		h.MissenseCount = count("missense_count")
		h.TruncatingCount = count("trunc_count")
		h.InframeCount = count("inframe_count")
		h.SpliceCount = count("splice_count")

		data[txID] = append(data[txID], h)
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("read hotspots file: %w", err)
	}

	// Sort each gene's hotspots by position and deduplicate
	for gene, spots := range data {
		sort.Slice(spots, func(i, j int) bool {
			return spots[i].Position < spots[j].Position
		})
		data[gene] = dedup(spots)
	}

	return &Store{data: data}, nil
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
	spots := s.data[transcriptID]
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
	case "splice":
		return 2
	case "3d":
		return 3
	default:
		return 4
	}
}

// dedup removes duplicate positions, keeping the entry with the best evidence.
// Prefers lower q-value, then "single residue" type over "3d".
func dedup(spots []Hotspot) []Hotspot {
	if len(spots) <= 1 {
		return spots
	}
	result := spots[:1]
	for i := 1; i < len(spots); i++ {
		if spots[i].Position == result[len(result)-1].Position {
			cur := &result[len(result)-1]
			// Prefer lower q-value, then higher-priority type
			if spots[i].QValue < cur.QValue ||
				(spots[i].QValue == cur.QValue && typePriority(spots[i].Type) < typePriority(cur.Type)) {
				*cur = spots[i]
			}
		} else {
			result = append(result, spots[i])
		}
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
