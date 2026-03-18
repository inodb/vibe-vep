// Package ensemblpred provides SIFT and PolyPhen-2 predictions from Ensembl's
// pre-computed SQLite database. The database stores compressed prediction matrices
// keyed by peptide MD5 hash, covering every possible amino acid substitution.
//
// Database: https://ftp.ensembl.org/pub/current_variation/pangenomes/Human/
// Format: SQLite with table predictions(md5, analysis, matrix)
// Binary: gzip-compressed, 3-byte header, then (position * 20 + aa_index) * 2 bytes
// Each prediction is 16-bit LE: top 2 bits = qualitative, bottom 10 bits = score/1000
package ensemblpred

import (
	"bytes"
	"compress/gzip"
	"database/sql"
	"encoding/binary"
	"fmt"
	"io"
	"strconv"
	"strings"
	"sync"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
	_ "modernc.org/sqlite"
)

// Amino acids in alphabetical order (Ensembl convention).
var aaOrder = [20]byte{'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}

// aaIndex returns the index of an amino acid in the Ensembl alphabetical order.
// Returns -1 if the amino acid is not found.
func aaIndex(aa byte) int {
	for i, a := range aaOrder {
		if a == aa {
			return i
		}
	}
	return -1
}

// Prediction holds a decoded SIFT or PolyPhen prediction.
type Prediction struct {
	Score float64
	Pred  string
}

// SIFT qualitative predictions (top 2 bits of the 16-bit value).
var siftPreds = [4]string{"tolerated", "deleterious", "tolerated_low_confidence", "deleterious_low_confidence"}

// PolyPhen HDIV qualitative predictions.
var pp2Preds = [4]string{"probably_damaging", "possibly_damaging", "benign", "unknown"}

// headerSize is the number of bytes at the start of the decompressed matrix to skip.
const headerSize = 3

// Store provides lookups against the Ensembl SIFT/PolyPhen SQLite database.
type Store struct {
	db       *sql.DB
	lookupPS *sql.Stmt

	// LRU cache: variants often hit the same transcript repeatedly.
	mu    sync.Mutex
	cache map[cacheKey][]byte // decompressed matrix
}

type cacheKey struct {
	md5      string
	analysis string
}

// Open opens the Ensembl predictions database.
func Open(dbPath string) (*Store, error) {
	db, err := sql.Open("sqlite", dbPath+"?mode=ro&_pragma=mmap_size%3D2147483648")
	if err != nil {
		return nil, fmt.Errorf("open ensembl predictions: %w", err)
	}

	ps, err := db.Prepare(`SELECT matrix FROM predictions WHERE md5 = ? AND analysis = ?`)
	if err != nil {
		db.Close()
		return nil, fmt.Errorf("prepare lookup: %w", err)
	}

	return &Store{
		db:       db,
		lookupPS: ps,
		cache:    make(map[cacheKey][]byte, 128),
	}, nil
}

// Close closes the database.
func (s *Store) Close() error {
	if s.lookupPS != nil {
		s.lookupPS.Close()
	}
	return s.db.Close()
}

// Lookup retrieves a prediction for a specific amino acid change.
// position is 1-based. altAA is the single-letter alternate amino acid.
// analysis is "sift" or "polyphen_humdiv".
func (s *Store) Lookup(md5, analysis string, position int, altAA byte) (Prediction, bool) {
	idx := aaIndex(altAA)
	if idx < 0 || position < 1 {
		return Prediction{}, false
	}

	matrix, ok := s.getMatrix(md5, analysis)
	if !ok {
		return Prediction{}, false
	}

	// Offset: 3-byte header, then (pos-1)*20 + aa_index, each entry is 2 bytes.
	offset := headerSize + ((position-1)*20+idx)*2
	if offset+2 > len(matrix) {
		return Prediction{}, false
	}

	raw := binary.LittleEndian.Uint16(matrix[offset : offset+2])

	// 0xFFFF means "reference amino acid" (no prediction).
	if raw == 0xFFFF {
		return Prediction{}, false
	}

	pred := int(raw >> 14)          // top 2 bits
	scoreRaw := int(raw & 0x03FF)   // bottom 10 bits
	score := float64(scoreRaw) / 1000.0

	var predStr string
	switch analysis {
	case "sift":
		predStr = siftPreds[pred]
	case "polyphen_humdiv":
		predStr = pp2Preds[pred]
	default:
		predStr = ""
	}

	return Prediction{Score: score, Pred: predStr}, true
}

// getMatrix retrieves and decompresses a prediction matrix, using a cache.
func (s *Store) getMatrix(md5, analysis string) ([]byte, bool) {
	key := cacheKey{md5, analysis}

	s.mu.Lock()
	if m, ok := s.cache[key]; ok {
		s.mu.Unlock()
		return m, true
	}
	s.mu.Unlock()

	var compressed []byte
	err := s.lookupPS.QueryRow(md5, analysis).Scan(&compressed)
	if err != nil {
		return nil, false
	}

	matrix, err := decompress(compressed)
	if err != nil {
		return nil, false
	}

	s.mu.Lock()
	// Simple bounded cache: evict all if too large.
	if len(s.cache) > 200 {
		s.cache = make(map[cacheKey][]byte, 128)
	}
	s.cache[key] = matrix
	s.mu.Unlock()

	return matrix, true
}

// decompress decompresses a gzip-compressed prediction matrix.
func decompress(data []byte) ([]byte, error) {
	r, err := gzip.NewReader(bytes.NewReader(data))
	if err != nil {
		return nil, err
	}
	defer r.Close()
	return io.ReadAll(r)
}

// Pre-built keys for Extra map.
const (
	extraKeySiftScore = "sift.score"
	extraKeySiftPred  = "sift.prediction"
	extraKeyPP2Score  = "polyphen.score"
	extraKeyPP2Pred   = "polyphen.prediction"
)

// Source implements annotate.AnnotationSource for Ensembl SIFT/PolyPhen predictions.
type Source struct {
	store *Store
}

// NewSource creates an AnnotationSource backed by the given Store.
func NewSource(store *Store) *Source {
	return &Source{store: store}
}

func (s *Source) Name() string                   { return "ensembl_predictions" }
func (s *Source) Version() string                 { return "115" }
func (s *Source) MatchLevel() annotate.MatchLevel { return annotate.MatchProteinPosition }
func (s *Source) Store() *Store                   { return s.store }

func (s *Source) Columns() []annotate.ColumnDef {
	return []annotate.ColumnDef{
		{Name: "sift.score", Description: "SIFT score (0-1, lower = more damaging)"},
		{Name: "sift.prediction", Description: "SIFT prediction (deleterious/tolerated)"},
		{Name: "polyphen.score", Description: "PolyPhen-2 HDIV score (0-1, higher = more damaging)"},
		{Name: "polyphen.prediction", Description: "PolyPhen-2 HDIV prediction (probably_damaging/possibly_damaging/benign)"},
	}
}

// Annotate adds SIFT and PolyPhen-2 predictions to missense annotations.
func (s *Source) Annotate(_ *vcf.Variant, anns []*annotate.Annotation) {
	for _, ann := range anns {
		if ann.PeptideMD5 == "" || ann.ProteinPosition == 0 || ann.AminoAcidChange == "" {
			continue
		}
		if !isMissense(ann.Consequence) {
			continue
		}

		altAA := extractAltAA(ann.AminoAcidChange)
		if altAA == 0 {
			continue
		}

		pos := int(ann.ProteinPosition)

		if p, ok := s.store.Lookup(ann.PeptideMD5, "sift", pos, altAA); ok {
			ann.SetExtraKey(extraKeySiftScore, formatScore(p.Score))
			ann.SetExtraKey(extraKeySiftPred, p.Pred)
		}

		if p, ok := s.store.Lookup(ann.PeptideMD5, "polyphen_humdiv", pos, altAA); ok {
			ann.SetExtraKey(extraKeyPP2Score, formatScore(p.Score))
			ann.SetExtraKey(extraKeyPP2Pred, p.Pred)
		}
	}
}

// extractAltAA extracts the alternate amino acid from an AminoAcidChange string.
// Format is like "G12C" — the last character is the alt AA.
func extractAltAA(change string) byte {
	if len(change) < 2 {
		return 0
	}
	aa := change[len(change)-1]
	// Must be an uppercase letter.
	if aa < 'A' || aa > 'Z' {
		return 0
	}
	return aa
}

// formatScore formats a float64 score as a 3-decimal string.
func formatScore(score float64) string {
	return strconv.FormatFloat(score, 'f', 3, 64)
}

// isMissense returns true if the consequence includes missense_variant.
func isMissense(consequence string) bool {
	for rest := consequence; rest != ""; {
		term := rest
		if i := strings.IndexByte(rest, ','); i >= 0 {
			term = rest[:i]
			rest = rest[i+1:]
		} else {
			rest = ""
		}
		if term == annotate.ConsequenceMissenseVariant {
			return true
		}
	}
	return false
}
