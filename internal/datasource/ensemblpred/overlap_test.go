package ensemblpred

import (
	"bufio"
	"compress/gzip"
	"crypto/md5"
	"encoding/hex"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
)

// TestProteinMD5Overlap checks what percentage of GENCODE protein sequences
// have matching MD5s in the Ensembl variation database.
func TestProteinMD5Overlap(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping overlap test in short mode")
	}

	home, err := os.UserHomeDir()
	if err != nil {
		t.Skip("cannot determine home dir")
	}
	cacheDir := filepath.Join(home, ".vibe-vep", "grch38")

	// Load Ensembl translation MD5s.
	md5Path := filepath.Join(cacheDir, "translation_md5.txt.gz")
	if _, err := os.Stat(md5Path); err != nil {
		t.Skipf("translation_md5.txt.gz not found at %s", md5Path)
	}

	ensemblMD5s := make(map[string]bool)
	f, err := os.Open(md5Path)
	if err != nil {
		t.Fatal(err)
	}
	gz, err := gzip.NewReader(f)
	if err != nil {
		t.Fatal(err)
	}
	scanner := bufio.NewScanner(gz)
	for scanner.Scan() {
		line := scanner.Text()
		if tab := strings.IndexByte(line, '\t'); tab >= 0 {
			ensemblMD5s[line[tab+1:]] = true
		}
	}
	gz.Close()
	f.Close()
	t.Logf("Ensembl translation MD5s: %d", len(ensemblMD5s))

	// Load GENCODE transcripts.
	gtfPath := filepath.Join(cacheDir, "gencode.v45.annotation.gtf.gz")
	fastaPath := filepath.Join(cacheDir, "gencode.v45.pc_transcripts.fa.gz")
	if _, err := os.Stat(gtfPath); err != nil {
		t.Skipf("GENCODE GTF not found at %s", gtfPath)
	}

	c := cache.New()
	loader := cache.NewGENCODELoader(gtfPath, fastaPath)
	canonicalPath := filepath.Join(cacheDir, "ensembl_biomart_canonical_transcripts_per_hgnc.txt")
	if _, err := os.Stat(canonicalPath); err == nil {
		msk, ens, _ := cache.LoadBiomartCanonicals(canonicalPath)
		loader.SetCanonicalOverrides(msk, ens)
	}
	if err := loader.Load(c); err != nil {
		t.Fatalf("load GENCODE: %v", err)
	}

	totalProteins := 0
	matched := 0
	missedCanonical := 0
	matchedCanonical := 0

	for _, chrom := range c.Chromosomes() {
		for _, tr := range c.FindTranscriptsByChrom(chrom) {
			if tr.CDSSequence == "" || len(tr.CDSSequence) < 3 {
				continue
			}
			// Translate CDS → protein.
			protein := translateCDS(tr.CDSSequence)
			if protein == "" {
				continue
			}
			totalProteins++
			h := md5.Sum([]byte(protein))
			md5hex := hex.EncodeToString(h[:])
			if ensemblMD5s[md5hex] {
				matched++
				if tr.IsCanonicalMSK {
					matchedCanonical++
				}
			} else if tr.IsCanonicalMSK {
				missedCanonical++
			}
		}
	}

	pct := float64(matched) / float64(totalProteins) * 100
	t.Logf("GENCODE proteins: %d", totalProteins)
	t.Logf("Matched Ensembl: %d (%.1f%%)", matched, pct)
	t.Logf("Missed: %d (%.1f%%)", totalProteins-matched, 100-pct)
	t.Logf("Canonical: %d matched, %d missed", matchedCanonical, missedCanonical)

	// We expect most proteins to match. Warn if coverage is low.
	if pct < 50 {
		t.Errorf("only %.1f%% of GENCODE proteins match Ensembl 115 — version mismatch too large", pct)
	}
}

// translateCDS translates a CDS DNA sequence to a protein sequence (single-letter AA codes).
// Stops at the first stop codon (not included). Returns "" if too short.
func translateCDS(cds string) string {
	if len(cds) < 3 {
		return ""
	}
	var b strings.Builder
	for i := 0; i+2 < len(cds); i += 3 {
		aa := annotate.TranslateCodon(cds[i : i+3])
		if aa == '*' {
			break
		}
		b.WriteByte(aa)
	}
	return b.String()
}
