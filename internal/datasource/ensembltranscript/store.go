package ensembltranscript

import (
	"database/sql"
	"encoding/json"
	"errors"
	"fmt"
	"os"
	"sort"
	"strings"
)

// Store answers the /ensembl endpoints from the built index. Records come back
// already encoded, so a handler can write one straight out.
type Store struct {
	db *sql.DB
}

// Open attaches to an index built by Build. It reports whether the file is
// there and usable rather than failing, so a bundle without this source still
// serves everything else.
func Open(dbPath string) (*Store, error) {
	if _, err := os.Stat(dbPath); err != nil {
		return nil, err
	}
	db, err := sql.Open("sqlite", dbPath+"?_pragma=mmap_size%3D268435456&mode=ro")
	if err != nil {
		return nil, fmt.Errorf("open ensembl transcript index: %w", err)
	}
	var version string
	if err := db.QueryRow(`SELECT value FROM meta WHERE key = 'schema_version'`).Scan(&version); err != nil {
		db.Close()
		return nil, fmt.Errorf("read ensembl transcript index metadata: %w", err)
	}
	if version != schemaVersion {
		db.Close()
		return nil, fmt.Errorf("ensembl transcript index is schema %s, this binary reads %s", version, schemaVersion)
	}
	return &Store{db: db}, nil
}

// Close releases the index.
func (s *Store) Close() error {
	if s == nil || s.db == nil {
		return nil
	}
	return s.db.Close()
}

// Count reports how many transcripts the index holds.
func (s *Store) Count() int {
	if s == nil {
		return 0
	}
	var n int
	if err := s.db.QueryRow(`SELECT value FROM meta WHERE key = 'transcripts'`).Scan(&n); err != nil {
		return 0
	}
	return n
}

// stripVersion drops the ".N" suffix an Ensembl id may carry, since the index
// is keyed on the bare id.
func stripVersion(id string) string {
	if i := strings.IndexByte(id, '.'); i >= 0 {
		return id[:i]
	}
	return id
}

// ByTranscriptID returns one transcript, with or without a version suffix on
// the id, reporting false when the index has no such transcript.
func (s *Store) ByTranscriptID(id string) (json.RawMessage, bool) {
	if s == nil {
		return nil, false
	}
	var doc []byte
	err := s.db.QueryRow(`SELECT doc FROM transcript WHERE transcript_id = ?`, stripVersion(id)).Scan(&doc)
	if errors.Is(err, sql.ErrNoRows) || err != nil {
		return nil, false
	}
	return doc, true
}

// ByTranscriptIDs returns the transcripts named, in the order asked for,
// skipping those the index does not hold.
func (s *Store) ByTranscriptIDs(ids []string) []json.RawMessage {
	out := make([]json.RawMessage, 0, len(ids))
	for _, id := range ids {
		if doc, ok := s.ByTranscriptID(id); ok {
			out = append(out, doc)
		}
	}
	return out
}

// query runs one of the indexed lookups and returns the documents, ordered by
// transcript id so that repeating a request gives the same answer.
func (s *Store) query(sqlText string, args ...any) []json.RawMessage {
	if s == nil {
		return nil
	}
	rows, err := s.db.Query(sqlText, args...)
	if err != nil {
		return nil
	}
	defer rows.Close()

	type entry struct {
		id  string
		doc []byte
	}
	var found []entry
	for rows.Next() {
		var e entry
		if err := rows.Scan(&e.id, &e.doc); err != nil {
			return nil
		}
		found = append(found, e)
	}
	sort.Slice(found, func(i, j int) bool { return found[i].id < found[j].id })
	out := make([]json.RawMessage, 0, len(found))
	for _, e := range found {
		out = append(out, e.doc)
	}
	return out
}

// ByGeneID returns every transcript of a gene.
func (s *Store) ByGeneID(geneID string) []json.RawMessage {
	return s.query(`SELECT transcript_id, doc FROM transcript WHERE gene_id = ?`, stripVersion(geneID))
}

// ByProteinID returns the transcripts translating to a protein.
func (s *Store) ByProteinID(proteinID string) []json.RawMessage {
	return s.query(`SELECT transcript_id, doc FROM transcript WHERE protein_id = ?`, stripVersion(proteinID))
}

// BySymbol returns every transcript a gene symbol claims.
func (s *Store) BySymbol(symbol string) []json.RawMessage {
	return s.query(`SELECT t.transcript_id, t.doc FROM transcript t
		JOIN transcript_symbol s ON s.transcript_id = t.transcript_id
		WHERE s.symbol = ?`, symbol)
}

// HugoSymbol returns the HGNC symbol a transcript's gene is approved under,
// or "" when the index does not name one.
//
// This is the symbol Ensembl and VEP report. A GENCODE release predating an
// HGNC assignment still carries the clone-based name -- AC005152.2 where the
// gene is now SOX9-AS1 -- and that name appears in no alias table, so it can
// only be resolved through the transcript.
func (s *Store) HugoSymbol(transcriptID string) string {
	doc, ok := s.ByTranscriptID(transcriptID)
	if !ok {
		return ""
	}
	var record struct {
		HugoSymbols []string `json:"hugoSymbols"`
	}
	if err := json.Unmarshal(doc, &record); err != nil || len(record.HugoSymbols) == 0 {
		return ""
	}
	return record.HugoSymbols[0]
}
