// Package ensembltranscript serves the Ensembl transcript records the
// /ensembl endpoints report.
//
// These cannot be derived from the GENCODE cache: CCDS and RefSeq accessions,
// exon identifiers and their versions, and Pfam domain ranges are not in the
// GTF. genome-nexus answers them from its ensembl.transcript collection, which
// genome-nexus-importer builds as ensembl_biomart_transcripts.json.gz, so that
// file is the source here too and the records match by construction.
//
// Each record is stored already in its wire shape, so serving one is a copy
// rather than a reconstruction.
package ensembltranscript

import (
	"bufio"
	"compress/gzip"
	"database/sql"
	"encoding/json"
	"fmt"
	"os"
	"strconv"
	"strings"

	_ "modernc.org/sqlite"
)

// Transcript is one Ensembl transcript as the /ensembl endpoints report it,
// mirroring genome-nexus's EnsemblTranscript. Absent values are omitted
// rather than zeroed, which is what genome-nexus does and what callers
// distinguishing "no CCDS" from "CCDS 0" need.
type Transcript struct {
	TranscriptID        string       `json:"transcriptId"`
	TranscriptIDVersion string       `json:"transcriptIdVersion,omitempty"`
	GeneID              string       `json:"geneId,omitempty"`
	RefseqMrnaID        string       `json:"refseqMrnaId,omitempty"`
	CcdsID              string       `json:"ccdsId,omitempty"`
	HugoSymbols         []string     `json:"hugoSymbols,omitempty"`
	ProteinID           string       `json:"proteinId,omitempty"`
	ProteinLength       *int         `json:"proteinLength,omitempty"`
	PfamDomains         []PfamDomain `json:"pfamDomains,omitempty"`
	Exons               []Exon       `json:"exons,omitempty"`
	UTRs                []UTR        `json:"utrs,omitempty"`
	UniprotID           string       `json:"uniprotId,omitempty"`
}

// PfamDomain is one Pfam domain's span in protein coordinates.
type PfamDomain struct {
	PfamDomainID    string `json:"pfamDomainId,omitempty"`
	PfamDomainStart *int   `json:"pfamDomainStart,omitempty"`
	PfamDomainEnd   *int   `json:"pfamDomainEnd,omitempty"`
}

// Exon is one exon of a transcript, in genomic coordinates.
type Exon struct {
	ExonID    string `json:"exonId,omitempty"`
	ExonStart *int   `json:"exonStart,omitempty"`
	ExonEnd   *int   `json:"exonEnd,omitempty"`
	Rank      *int   `json:"rank,omitempty"`
	Strand    *int   `json:"strand,omitempty"`
	Version   *int   `json:"version,omitempty"`
}

// UTR is one untranslated region. genome-nexus's model keeps only these four
// fields, although the source record carries an id, a rank and a version.
type UTR struct {
	Type   string `json:"type,omitempty"`
	Start  *int   `json:"start,omitempty"`
	End    *int   `json:"end,omitempty"`
	Strand *int   `json:"strand,omitempty"`
}

// sourceRecord is one line of ensembl_biomart_transcripts.json.gz. Its numbers
// are written as strings, and a field with no value is null rather than absent.
type sourceRecord struct {
	TranscriptStableID  string   `json:"transcript_stable_id"`
	TranscriptIDVersion string   `json:"transcript_id_version"`
	GeneStableID        string   `json:"gene_stable_id"`
	ProteinStableID     string   `json:"protein_stable_id"`
	ProteinLength       *int     `json:"protein_length"`
	RefseqMrnaID        string   `json:"refseq_mrna_id"`
	CcdsID              string   `json:"ccds_id"`
	UniprotID           string   `json:"uniprot_id"`
	HgncSymbols         []string `json:"hgnc_symbols"`
	Domains             []struct {
		PfamDomainID    string `json:"pfam_domain_id"`
		PfamDomainStart string `json:"pfam_domain_start"`
		PfamDomainEnd   string `json:"pfam_domain_end"`
	} `json:"domains"`
	Exons []struct {
		ID      string `json:"id"`
		Start   string `json:"start"`
		End     string `json:"end"`
		Rank    string `json:"rank"`
		Strand  string `json:"strand"`
		Version string `json:"version"`
	} `json:"exons"`
	UTRs []struct {
		Type   string `json:"type"`
		Start  string `json:"start"`
		End    string `json:"end"`
		Strand string `json:"strand"`
	} `json:"utrs"`
}

// atoi reads one of the source file's string-encoded numbers, reporting
// nothing where the field is empty or unparseable so the value is omitted
// rather than reported as zero.
func atoi(s string) *int {
	n, err := strconv.Atoi(strings.TrimSpace(s))
	if err != nil {
		return nil
	}
	return &n
}

// toWire converts a source record into the shape the endpoints serve.
func (r sourceRecord) toWire() Transcript {
	t := Transcript{
		TranscriptID:        r.TranscriptStableID,
		TranscriptIDVersion: r.TranscriptIDVersion,
		GeneID:              r.GeneStableID,
		RefseqMrnaID:        r.RefseqMrnaID,
		CcdsID:              r.CcdsID,
		HugoSymbols:         r.HgncSymbols,
		ProteinID:           r.ProteinStableID,
		ProteinLength:       r.ProteinLength,
		UniprotID:           r.UniprotID,
	}
	for _, d := range r.Domains {
		t.PfamDomains = append(t.PfamDomains, PfamDomain{
			PfamDomainID:    d.PfamDomainID,
			PfamDomainStart: atoi(d.PfamDomainStart),
			PfamDomainEnd:   atoi(d.PfamDomainEnd),
		})
	}
	for _, e := range r.Exons {
		t.Exons = append(t.Exons, Exon{
			ExonID:    e.ID,
			ExonStart: atoi(e.Start),
			ExonEnd:   atoi(e.End),
			Rank:      atoi(e.Rank),
			Strand:    atoi(e.Strand),
			Version:   atoi(e.Version),
		})
	}
	for _, u := range r.UTRs {
		t.UTRs = append(t.UTRs, UTR{
			Type:   u.Type,
			Start:  atoi(u.Start),
			End:    atoi(u.End),
			Strand: atoi(u.Strand),
		})
	}
	return t
}

// schemaVersion rises when the stored shape changes, so that an index built by
// an older binary is rebuilt rather than served.
const schemaVersion = "1"

// Build reads ensembl_biomart_transcripts.json.gz into dbPath, replacing any
// index already there.
func Build(sourcePath, dbPath string, logf func(string, ...any)) error {
	f, err := os.Open(sourcePath)
	if err != nil {
		return fmt.Errorf("open ensembl transcripts: %w", err)
	}
	defer f.Close()
	gz, err := gzip.NewReader(f)
	if err != nil {
		return fmt.Errorf("gunzip ensembl transcripts: %w", err)
	}
	defer gz.Close()

	for _, suffix := range []string{"", "-wal", "-shm"} {
		os.Remove(dbPath + suffix)
	}
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

	if _, err := db.Exec(`
		CREATE TABLE transcript (
			transcript_id TEXT PRIMARY KEY,
			gene_id       TEXT,
			protein_id    TEXT,
			doc           BLOB NOT NULL
		) WITHOUT ROWID;
		CREATE TABLE transcript_symbol (
			symbol        TEXT NOT NULL,
			transcript_id TEXT NOT NULL
		);
		CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT) WITHOUT ROWID;
	`); err != nil {
		return fmt.Errorf("create schema: %w", err)
	}

	tx, err := db.Begin()
	if err != nil {
		return err
	}
	defer tx.Rollback()

	insertTx, err := tx.Prepare(`INSERT OR REPLACE INTO transcript
		(transcript_id, gene_id, protein_id, doc) VALUES (?, ?, ?, ?)`)
	if err != nil {
		return err
	}
	defer insertTx.Close()
	insertSym, err := tx.Prepare(`INSERT INTO transcript_symbol (symbol, transcript_id) VALUES (?, ?)`)
	if err != nil {
		return err
	}
	defer insertSym.Close()

	scanner := bufio.NewScanner(gz)
	scanner.Buffer(make([]byte, 0, 1<<20), 1<<24)
	var rows int
	for scanner.Scan() {
		line := scanner.Bytes()
		if len(strings.TrimSpace(string(line))) == 0 {
			continue
		}
		var rec sourceRecord
		if err := json.Unmarshal(line, &rec); err != nil {
			return fmt.Errorf("parse record %d: %w", rows+1, err)
		}
		if rec.TranscriptStableID == "" {
			continue
		}
		doc, err := json.Marshal(rec.toWire())
		if err != nil {
			return fmt.Errorf("encode record %s: %w", rec.TranscriptStableID, err)
		}
		if _, err := insertTx.Exec(rec.TranscriptStableID, rec.GeneStableID, rec.ProteinStableID, doc); err != nil {
			return fmt.Errorf("insert %s: %w", rec.TranscriptStableID, err)
		}
		for _, sym := range rec.HgncSymbols {
			if sym == "" {
				continue
			}
			if _, err := insertSym.Exec(sym, rec.TranscriptStableID); err != nil {
				return fmt.Errorf("insert symbol %s: %w", sym, err)
			}
		}
		rows++
	}
	if err := scanner.Err(); err != nil {
		return fmt.Errorf("read ensembl transcripts: %w", err)
	}
	if rows == 0 {
		return fmt.Errorf("no transcripts read from %s", sourcePath)
	}

	for key, value := range map[string]string{
		"schema_version": schemaVersion,
		"transcripts":    strconv.Itoa(rows),
	} {
		if _, err := tx.Exec(`INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)`, key, value); err != nil {
			return err
		}
	}
	if err := tx.Commit(); err != nil {
		return fmt.Errorf("commit: %w", err)
	}

	if _, err := db.Exec(`
		CREATE INDEX idx_transcript_gene ON transcript(gene_id);
		CREATE INDEX idx_transcript_protein ON transcript(protein_id);
		CREATE INDEX idx_symbol ON transcript_symbol(symbol);
	`); err != nil {
		return fmt.Errorf("create indexes: %w", err)
	}
	if logf != nil {
		logf("built ensembl transcript index: %d transcripts", rows)
	}
	return nil
}
