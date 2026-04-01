// Package ensemblpred provides SIFT and PolyPhen-2 predictions from Ensembl's
// variation database. Prediction matrices are downloaded as MySQL dump TSVs and
// built into a local SQLite database keyed by (peptide_md5, analysis).
//
// Source: https://ftp.ensembl.org/pub/release-115/mysql/homo_sapiens_variation_115_38/
// Binary: gzip-compressed, 3-byte header, then (position * 20 + aa_index) * 2 bytes
// Each prediction is 16-bit LE: top 2 bits = qualitative, bottom 10 bits = score/1000
package ensemblpred

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"database/sql"
	"encoding/binary"
	"fmt"
	"io"
	"os"
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

	mu         sync.Mutex
	cache      map[cacheKey][]byte // decompressed matrix (LRU)
	compressed map[cacheKey][]byte // preloaded compressed data (nil if not preloaded)
	preloaded  bool                // true if compressed data is in memory
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

// Preload reads all prediction matrices from SQLite into memory so that
// subsequent lookups don't hit the database. This trades ~200-500MB RAM
// for instant lookups, which is critical on networked filesystems like EFS
// where each SQLite read has high latency.
// Preload reads all compressed prediction matrices from SQLite into memory.
// Subsequent lookups decompress on demand (fast) instead of hitting SQLite (slow on EFS).
// Compressed data is ~100-200MB; decompression happens per-lookup with an LRU cache.
func (s *Store) Preload() (int, error) {
	rows, err := s.db.Query(`SELECT md5, analysis, matrix FROM predictions`)
	if err != nil {
		return 0, fmt.Errorf("query predictions: %w", err)
	}
	defer rows.Close()

	compressed := make(map[cacheKey][]byte, 200000)
	n := 0
	for rows.Next() {
		var md5, analysis string
		var data []byte
		if err := rows.Scan(&md5, &analysis, &data); err != nil {
			return n, fmt.Errorf("scan row %d: %w", n, err)
		}
		// Store compressed data as-is (no decompression yet).
		compressed[cacheKey{md5, analysis}] = data
		n++
	}
	if err := rows.Err(); err != nil {
		return n, fmt.Errorf("iterate predictions: %w", err)
	}

	s.mu.Lock()
	s.compressed = compressed
	s.preloaded = true
	s.mu.Unlock()

	// Close the database — no longer needed.
	s.lookupPS.Close()
	s.db.Close()
	s.db = nil
	s.lookupPS = nil

	return n, nil
}

// Close closes the database (no-op if preloaded).
func (s *Store) Close() error {
	if s.lookupPS != nil {
		s.lookupPS.Close()
	}
	if s.db != nil {
		return s.db.Close()
	}
	return nil
}

// BuildSources holds paths to the Ensembl MySQL dump files for building the predictions DB.
type BuildSources struct {
	TranslationMD5TSV  string // translation_md5.txt.gz
	PredictionsTSV     string // protein_function_predictions.txt.gz
}

// Ready returns true if the SQLite DB exists and is newer than source files.
func Ready(dbPath string, sources BuildSources) bool {
	dbInfo, err := os.Stat(dbPath)
	if err != nil || dbInfo.Size() == 0 {
		return false
	}
	dbMod := dbInfo.ModTime()

	for _, src := range []string{sources.TranslationMD5TSV, sources.PredictionsTSV} {
		if src == "" {
			continue
		}
		srcInfo, err := os.Stat(src)
		if err != nil {
			continue
		}
		if srcInfo.ModTime().After(dbMod) {
			return false
		}
	}

	// Quick integrity check.
	db, err := sql.Open("sqlite", dbPath+"?mode=ro")
	if err != nil {
		return false
	}
	defer db.Close()
	var n int
	return db.QueryRow("SELECT 1 FROM predictions LIMIT 1").Scan(&n) == nil
}

// Analysis attrib IDs from the Ensembl variation database.
const (
	attribSIFT          = 267
	attribPolyPhenHumVar = 268
	attribPolyPhenHumDiv = 269
)

// analysisName maps Ensembl attrib IDs to analysis names.
var analysisName = map[int]string{
	attribSIFT:          "sift",
	attribPolyPhenHumVar: "polyphen_humvar",
	attribPolyPhenHumDiv: "polyphen_humdiv",
}

// Build creates the SQLite predictions database from Ensembl MySQL dump files.
func Build(dbPath string, sources BuildSources, logf func(string, ...any)) error {
	if sources.TranslationMD5TSV == "" || sources.PredictionsTSV == "" {
		return fmt.Errorf("both TranslationMD5TSV and PredictionsTSV are required")
	}

	// Load translation_md5 mapping: id → md5 hex string.
	logf("loading translation MD5 mapping from %s", sources.TranslationMD5TSV)
	md5Map, err := loadTranslationMD5(sources.TranslationMD5TSV)
	if err != nil {
		return fmt.Errorf("load translation MD5: %w", err)
	}
	logf("loaded %d translation MD5 mappings", len(md5Map))

	// Create SQLite.
	os.Remove(dbPath)
	db, err := sql.Open("sqlite", dbPath)
	if err != nil {
		return fmt.Errorf("create sqlite: %w", err)
	}
	defer db.Close()

	for _, pragma := range []string{
		"PRAGMA journal_mode = OFF",
		"PRAGMA synchronous = OFF",
		"PRAGMA temp_store = MEMORY",
		"PRAGMA cache_size = -64000",
		"PRAGMA page_size = 8192",
	} {
		if _, err := db.Exec(pragma); err != nil {
			return fmt.Errorf("set pragma %q: %w", pragma, err)
		}
	}

	if _, err := db.Exec(`CREATE TABLE predictions (
		md5 TEXT NOT NULL,
		analysis TEXT NOT NULL,
		matrix BLOB NOT NULL,
		PRIMARY KEY (md5, analysis)
	) WITHOUT ROWID`); err != nil {
		return fmt.Errorf("create table: %w", err)
	}

	// Load predictions.
	logf("loading predictions from %s", sources.PredictionsTSV)
	n, err := loadPredictions(db, sources.PredictionsTSV, md5Map)
	if err != nil {
		return fmt.Errorf("load predictions: %w", err)
	}
	logf("loaded %d prediction matrices", n)

	return nil
}

// loadTranslationMD5 reads the translation_md5.txt.gz MySQL dump.
// Format: translation_md5_id\ttranslation_md5
func loadTranslationMD5(path string) (map[int]string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var reader io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return nil, err
		}
		defer gz.Close()
		reader = gz
	}

	m := make(map[int]string, 210000)
	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	for scanner.Scan() {
		line := scanner.Text()
		tab := strings.IndexByte(line, '\t')
		if tab < 0 {
			continue
		}
		id, err := strconv.Atoi(line[:tab])
		if err != nil {
			continue
		}
		m[id] = line[tab+1:]
	}
	return m, scanner.Err()
}

// loadPredictions reads protein_function_predictions.txt.gz and inserts into SQLite.
// Format: translation_md5_id\tanalysis_attrib_id\tprediction_matrix_blob
// The blob is MySQL-escaped binary (backslash sequences: \0, \n, \t, \\).
func loadPredictions(db *sql.DB, path string, md5Map map[int]string) (int64, error) {
	f, err := os.Open(path)
	if err != nil {
		return 0, err
	}
	defer f.Close()

	var reader io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return 0, err
		}
		defer gz.Close()
		reader = gz
	}

	tx, err := db.Begin()
	if err != nil {
		return 0, err
	}
	defer tx.Rollback()

	stmt, err := tx.Prepare(`INSERT OR IGNORE INTO predictions (md5, analysis, matrix) VALUES (?, ?, ?)`)
	if err != nil {
		return 0, err
	}
	defer stmt.Close()

	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 16*1024*1024), 16*1024*1024) // large buffer for binary blobs

	var count int64
	for scanner.Scan() {
		line := scanner.Bytes()

		// Parse: translation_md5_id \t analysis_attrib_id \t blob
		tab1 := bytes.IndexByte(line, '\t')
		if tab1 < 0 {
			continue
		}
		tab2 := bytes.IndexByte(line[tab1+1:], '\t')
		if tab2 < 0 {
			continue
		}
		tab2 += tab1 + 1

		mdID, err := strconv.Atoi(string(line[:tab1]))
		if err != nil {
			continue
		}
		attribID, err := strconv.Atoi(string(line[tab1+1 : tab2]))
		if err != nil {
			continue
		}

		md5hex, ok := md5Map[mdID]
		if !ok {
			continue
		}
		analysis, ok := analysisName[attribID]
		if !ok {
			continue // skip polyphen_humvar and unknown analyses
		}
		// Only keep sift and polyphen_humdiv.
		if analysis != "sift" && analysis != "polyphen_humdiv" {
			continue
		}

		blob := unescapeMySQL(line[tab2+1:])

		if _, err := stmt.Exec(md5hex, analysis, blob); err != nil {
			return 0, fmt.Errorf("insert prediction: %w", err)
		}
		count++
	}
	if err := scanner.Err(); err != nil {
		return 0, err
	}

	if err := tx.Commit(); err != nil {
		return 0, err
	}
	return count, nil
}

// unescapeMySQL decodes MySQL's backslash-escaped binary format in TSV dumps.
// Sequences: \0 → 0x00, \n → 0x0A, \t → 0x09, \\ → 0x5C, \r → 0x0D
func unescapeMySQL(data []byte) []byte {
	if bytes.IndexByte(data, '\\') < 0 {
		return data // fast path: no escapes
	}

	out := make([]byte, 0, len(data))
	for i := 0; i < len(data); i++ {
		if data[i] == '\\' && i+1 < len(data) {
			i++
			switch data[i] {
			case '0':
				out = append(out, 0)
			case 'n':
				out = append(out, '\n')
			case 't':
				out = append(out, '\t')
			case '\\':
				out = append(out, '\\')
			case 'r':
				out = append(out, '\r')
			default:
				out = append(out, data[i])
			}
		} else {
			out = append(out, data[i])
		}
	}
	return out
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

	// Check decompressed cache first.
	s.mu.Lock()
	m, ok := s.cache[key]
	s.mu.Unlock()
	if ok {
		return m, true
	}

	// Try preloaded compressed data (fast — in-memory decompress).
	if s.preloaded {
		s.mu.Lock()
		comp, ok := s.compressed[key]
		s.mu.Unlock()
		if !ok {
			return nil, false
		}
		matrix, err := decompress(comp)
		if err != nil {
			return nil, false
		}
		s.mu.Lock()
		if len(s.cache) > 10000 {
			s.cache = make(map[cacheKey][]byte, 1000)
		}
		s.cache[key] = matrix
		s.mu.Unlock()
		return matrix, true
	}

	// Fallback: query SQLite (slow on EFS).
	if s.lookupPS == nil {
		return nil, false
	}
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
