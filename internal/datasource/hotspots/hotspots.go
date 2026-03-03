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
	Position int64   // amino acid position
	Type     string  // "single residue", "in-frame indel", "3d", "splice"
	QValue   float64 // statistical significance (q-value)
}

// Store holds hotspot data keyed by gene symbol.
// Each gene's hotspots are sorted by amino acid position for binary search.
type Store struct {
	data map[string][]Hotspot // gene symbol → sorted hotspots
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
	colIdx := indexColumns(header, "hugo_symbol", "amino_acid_position", "type", "q_value")
	if colIdx["hugo_symbol"] < 0 || colIdx["amino_acid_position"] < 0 {
		return nil, fmt.Errorf("missing required columns: hugo_symbol, amino_acid_position")
	}

	data := make(map[string][]Hotspot)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) <= colIdx["hugo_symbol"] || len(fields) <= colIdx["amino_acid_position"] {
			continue
		}

		gene := fields[colIdx["hugo_symbol"]]
		posStr := fields[colIdx["amino_acid_position"]]
		pos, err := strconv.ParseInt(posStr, 10, 64)
		if err != nil {
			continue // skip rows with non-numeric positions
		}

		h := Hotspot{Position: pos}

		if idx := colIdx["type"]; idx >= 0 && idx < len(fields) {
			h.Type = fields[idx]
		}
		if idx := colIdx["q_value"]; idx >= 0 && idx < len(fields) && fields[idx] != "" {
			if v, err := strconv.ParseFloat(fields[idx], 64); err == nil {
				h.QValue = v
			}
		}

		data[gene] = append(data[gene], h)
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

// Lookup checks if a gene+position is a known hotspot.
func (s *Store) Lookup(gene string, proteinPosition int64) (Hotspot, bool) {
	spots := s.data[gene]
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

// GeneCount returns the number of genes with hotspots.
func (s *Store) GeneCount() int {
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
