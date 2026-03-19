package ensemblpred

import (
	"bytes"
	"compress/gzip"
	"crypto/md5"
	"database/sql"
	"encoding/binary"
	"encoding/hex"
	"os"
	"path/filepath"
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	_ "modernc.org/sqlite"
)

func TestAAIndex(t *testing.T) {
	tests := []struct {
		aa   byte
		want int
	}{
		{'A', 0}, {'C', 1}, {'D', 2}, {'Y', 19},
		{'B', -1}, {'Z', -1}, {'a', -1},
	}
	for _, tt := range tests {
		if got := aaIndex(tt.aa); got != tt.want {
			t.Errorf("aaIndex(%q) = %d, want %d", tt.aa, got, tt.want)
		}
	}
}

func TestExtractAltAA(t *testing.T) {
	tests := []struct {
		change string
		want   byte
	}{
		{"G12C", 'C'},
		{"V600E", 'E'},
		{"", 0},
		{"A", 0},
	}
	for _, tt := range tests {
		if got := extractAltAA(tt.change); got != tt.want {
			t.Errorf("extractAltAA(%q) = %q, want %q", tt.change, got, tt.want)
		}
	}
}

// buildTestMatrix creates a gzip-compressed prediction matrix with a known value
// at the given position and amino acid index.
func buildTestMatrix(position, aaIdx int, predBits uint16, scoreRaw uint16) []byte {
	// Calculate total size needed.
	totalPositions := position + 1
	matrixSize := headerSize + totalPositions*20*2

	raw := make([]byte, matrixSize)
	// Set header bytes (3 bytes, content doesn't matter for our decoding).
	raw[0] = 0
	raw[1] = 0
	raw[2] = 0

	// Fill with 0xFFFF (no prediction) everywhere.
	for i := headerSize; i < len(raw)-1; i += 2 {
		raw[i] = 0xFF
		raw[i+1] = 0xFF
	}

	// Encode the prediction at the target position.
	value := (predBits << 14) | scoreRaw
	offset := headerSize + (position*20+aaIdx)*2
	binary.LittleEndian.PutUint16(raw[offset:], value)

	// Gzip compress.
	var buf bytes.Buffer
	w := gzip.NewWriter(&buf)
	w.Write(raw)
	w.Close()
	return buf.Bytes()
}

func setupTestDB(t *testing.T) *Store {
	t.Helper()

	db, err := sql.Open("sqlite", ":memory:")
	if err != nil {
		t.Fatal(err)
	}

	_, err = db.Exec(`CREATE TABLE predictions (md5 TEXT, analysis TEXT, matrix BLOB)`)
	if err != nil {
		t.Fatal(err)
	}
	_, err = db.Exec(`CREATE INDEX md5_idx ON predictions(md5)`)
	if err != nil {
		t.Fatal(err)
	}

	// Insert a SIFT prediction: position 12 (0-based=11), alt AA = C (index 1).
	// SIFT pred=1 (deleterious), score=50 (0.050).
	siftMatrix := buildTestMatrix(11, 1, 1, 50)
	_, err = db.Exec(`INSERT INTO predictions (md5, analysis, matrix) VALUES (?, ?, ?)`,
		"abc123", "sift", siftMatrix)
	if err != nil {
		t.Fatal(err)
	}

	// Insert a PolyPhen prediction at same position.
	// PP2 pred=0 (probably_damaging), score=998 (0.998).
	pp2Matrix := buildTestMatrix(11, 1, 0, 998)
	_, err = db.Exec(`INSERT INTO predictions (md5, analysis, matrix) VALUES (?, ?, ?)`,
		"abc123", "polyphen_humdiv", pp2Matrix)
	if err != nil {
		t.Fatal(err)
	}

	ps, err := db.Prepare(`SELECT matrix FROM predictions WHERE md5 = ? AND analysis = ?`)
	if err != nil {
		t.Fatal(err)
	}

	return &Store{
		db:       db,
		lookupPS: ps,
		cache:    make(map[cacheKey][]byte, 128),
	}
}

func TestLookupSIFT(t *testing.T) {
	store := setupTestDB(t)
	defer store.Close()

	// Position 12 (1-based), alt AA = C.
	p, ok := store.Lookup("abc123", "sift", 12, 'C')
	if !ok {
		t.Fatal("expected SIFT hit")
	}
	if p.Pred != "deleterious" {
		t.Errorf("SIFT pred=%q, want %q", p.Pred, "deleterious")
	}
	if diff := p.Score - 0.050; diff < -0.001 || diff > 0.001 {
		t.Errorf("SIFT score=%v, want 0.050", p.Score)
	}
}

func TestLookupPolyPhen(t *testing.T) {
	store := setupTestDB(t)
	defer store.Close()

	p, ok := store.Lookup("abc123", "polyphen_humdiv", 12, 'C')
	if !ok {
		t.Fatal("expected PP2 hit")
	}
	if p.Pred != "probably_damaging" {
		t.Errorf("PP2 pred=%q, want %q", p.Pred, "probably_damaging")
	}
	if diff := p.Score - 0.998; diff < -0.001 || diff > 0.001 {
		t.Errorf("PP2 score=%v, want 0.998", p.Score)
	}
}

func TestLookupMiss(t *testing.T) {
	store := setupTestDB(t)
	defer store.Close()

	// Wrong MD5.
	_, ok := store.Lookup("nonexistent", "sift", 12, 'C')
	if ok {
		t.Error("expected miss for wrong MD5")
	}

	// Wrong position (too large).
	_, ok = store.Lookup("abc123", "sift", 9999, 'C')
	if ok {
		t.Error("expected miss for out-of-range position")
	}

	// Reference AA position returns 0xFFFF.
	_, ok = store.Lookup("abc123", "sift", 12, 'A')
	if ok {
		t.Error("expected miss for reference AA (0xFFFF)")
	}
}

func TestSourceAnnotate(t *testing.T) {
	store := setupTestDB(t)
	defer store.Close()

	src := NewSource(store)

	anns := []*annotate.Annotation{
		{
			Consequence:     "missense_variant",
			ProteinPosition: 12,
			AminoAcidChange: "G12C",
			PeptideMD5:      "abc123",
		},
		{
			Consequence:     "synonymous_variant",
			ProteinPosition: 12,
			AminoAcidChange: "",
			PeptideMD5:      "abc123",
		},
	}

	src.Annotate(nil, anns)

	// Missense should have SIFT + PP2.
	if got := anns[0].GetExtraKey("sift.prediction"); got != "deleterious" {
		t.Errorf("missense sift.prediction=%q, want %q", got, "deleterious")
	}
	if got := anns[0].GetExtraKey("polyphen.prediction"); got != "probably_damaging" {
		t.Errorf("missense polyphen.prediction=%q, want %q", got, "probably_damaging")
	}
	if got := anns[0].GetExtraKey("sift.score"); got != "0.050" {
		t.Errorf("missense sift.score=%q, want %q", got, "0.050")
	}

	// Synonymous should have nothing.
	if got := anns[1].GetExtraKey("sift.prediction"); got != "" {
		t.Errorf("synonymous should not have sift.prediction, got %q", got)
	}
}

// TestRealEnsemblDB tests against the real Ensembl SIFT/PolyPhen database.
// Skipped if the database is not present.
func TestRealEnsemblDB(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping real DB test in short mode")
	}

	home, err := os.UserHomeDir()
	if err != nil {
		t.Skip("cannot determine home dir")
	}
	dbPath := filepath.Join(home, ".vibe-vep", "grch38", "ensembl_sift_polyphen.db")
	if _, err := os.Stat(dbPath); err != nil {
		t.Skipf("Ensembl predictions DB not found at %s", dbPath)
	}

	store, err := Open(dbPath)
	if err != nil {
		t.Fatalf("open: %v", err)
	}
	defer store.Close()

	// KRAS protein (Ensembl canonical, 189 AA).
	krasSeq := "MTEYKLVVVGAVGVGKSALTIQLIQNHFVDEYDPTIEDSYNTKQFVQTVDSSSYRKKVDKDLSALLQHSLPEESFRKSVDHAKDDPMVEGSGGAALEEDTAVDFLHRIEMRQKELDKFNSHECSHFQPINNLKGESIPSEPASPAKDCAKDSNRGCTFNSSSSAKRDDSSPFLDSHQPKPNKKNTDCSTRILDTAGQEEFGVEQSGDDNYTEDEESTF"
	h := md5.Sum([]byte(krasSeq))
	krasMD5 := hex.EncodeToString(h[:])
	t.Logf("KRAS MD5: %s", krasMD5)

	// G12C: position 12, alt AA = C — should be deleterious/damaging.
	sift, ok := store.Lookup(krasMD5, "sift", 12, 'C')
	if !ok {
		t.Fatal("expected SIFT hit for KRAS G12C")
	}
	t.Logf("SIFT G12C: score=%.3f pred=%s", sift.Score, sift.Pred)
	if sift.Pred != "deleterious" {
		t.Errorf("SIFT G12C pred=%q, want deleterious", sift.Pred)
	}

	pp2, ok := store.Lookup(krasMD5, "polyphen_humdiv", 12, 'C')
	if !ok {
		t.Fatal("expected PP2 hit for KRAS G12C")
	}
	t.Logf("PP2  G12C: score=%.3f pred=%s", pp2.Score, pp2.Pred)
	if pp2.Pred != "probably_damaging" {
		t.Errorf("PP2 G12C pred=%q, want probably_damaging", pp2.Pred)
	}

	// Miss: nonexistent MD5.
	_, ok = store.Lookup("0000000000000000000000000000000", "sift", 1, 'A')
	if ok {
		t.Error("expected miss for nonexistent MD5")
	}
}
