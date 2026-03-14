package genomicindex

import (
	"database/sql"
	"os"
	"path/filepath"
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
	_ "modernc.org/sqlite"
)

// setupTestDB creates a small test database with known variants stored in
// canonical (MAF-style) form — no VCF anchor bases, positions adjusted.
func setupTestDB(t *testing.T) string {
	t.Helper()
	dir := t.TempDir()
	dbPath := filepath.Join(dir, "test.sqlite")

	db, err := sql.Open("sqlite", dbPath)
	if err != nil {
		t.Fatal(err)
	}
	defer db.Close()

	_, err = db.Exec(`CREATE TABLE genomic_annotations (
		chrom TEXT NOT NULL,
		pos INTEGER NOT NULL,
		ref TEXT NOT NULL,
		alt TEXT NOT NULL,
		am_score REAL NOT NULL DEFAULT 0,
		am_class TEXT NOT NULL DEFAULT '',
		cv_clnsig TEXT NOT NULL DEFAULT '',
		cv_revstat TEXT NOT NULL DEFAULT '',
		cv_clndn TEXT NOT NULL DEFAULT '',
		sig_mut_status TEXT NOT NULL DEFAULT '',
		sig_count TEXT NOT NULL DEFAULT '',
		sig_freq TEXT NOT NULL DEFAULT '',
		gnomad_af TEXT NOT NULL DEFAULT '',
		gnomad_ac TEXT NOT NULL DEFAULT '',
		gnomad_an TEXT NOT NULL DEFAULT '',
		gnomad_nhomalt TEXT NOT NULL DEFAULT '',
		gnomad_version TEXT NOT NULL DEFAULT '',
		PRIMARY KEY (chrom, pos, ref, alt)
	) WITHOUT ROWID`)
	if err != nil {
		t.Fatal(err)
	}

	// SNP: AM-only variant (KRAS G12C)
	_, err = db.Exec(`INSERT INTO genomic_annotations (chrom, pos, ref, alt, am_score, am_class) VALUES ('12', 25245350, 'C', 'A', 0.9876, 'likely_pathogenic')`)
	if err != nil {
		t.Fatal(err)
	}

	// SNP: ClinVar-only variant
	_, err = db.Exec(`INSERT INTO genomic_annotations (chrom, pos, ref, alt, cv_clnsig, cv_revstat, cv_clndn) VALUES ('1', 12345, 'A', 'G', 'Pathogenic', 'reviewed_by_expert_panel', 'Some_disease')`)
	if err != nil {
		t.Fatal(err)
	}

	// SNP: Combined AM + ClinVar variant
	_, err = db.Exec(`INSERT INTO genomic_annotations (chrom, pos, ref, alt, am_score, am_class, cv_clnsig, cv_revstat) VALUES ('7', 140753336, 'A', 'T', 0.5432, 'ambiguous', 'Uncertain_significance', 'criteria_provided,_single_submitter')`)
	if err != nil {
		t.Fatal(err)
	}

	// SNP: SIGNAL variant
	_, err = db.Exec(`INSERT INTO genomic_annotations (chrom, pos, ref, alt, sig_mut_status, sig_count, sig_freq) VALUES ('17', 7577538, 'C', 'T', 'germline', '5', '0.00123')`)
	if err != nil {
		t.Fatal(err)
	}

	// Deletion: ClinVar, stored in canonical form (pos=first deleted base, no anchor)
	// Original VCF: POS=66926 REF=AG ALT=A → normalized: pos=66927 ref=G alt=""
	_, err = db.Exec(`INSERT INTO genomic_annotations (chrom, pos, ref, alt, cv_clnsig) VALUES ('1', 66927, 'G', '', 'Pathogenic')`)
	if err != nil {
		t.Fatal(err)
	}

	// Multi-base deletion: ClinVar
	// Original VCF: POS=930328 REF=GAA ALT=G → normalized: pos=930329 ref=AA alt=""
	_, err = db.Exec(`INSERT INTO genomic_annotations (chrom, pos, ref, alt, cv_clnsig) VALUES ('1', 930329, 'AA', '', 'Likely_pathogenic')`)
	if err != nil {
		t.Fatal(err)
	}

	// Insertion: ClinVar
	// Original VCF: POS=809284 REF=T ALT=TGGTCAATCA → normalized: pos=809284 ref="" alt="GGTCAATCA"
	_, err = db.Exec(`INSERT INTO genomic_annotations (chrom, pos, ref, alt, cv_clnsig) VALUES ('1', 809284, '', 'GGTCAATCA', 'Uncertain_significance')`)
	if err != nil {
		t.Fatal(err)
	}

	// 1-base insertion: ClinVar
	// Original VCF: POS=935908 REF=A ALT=AG → normalized: pos=935908 ref="" alt="G"
	_, err = db.Exec(`INSERT INTO genomic_annotations (chrom, pos, ref, alt, cv_clnsig) VALUES ('1', 935908, '', 'G', 'Benign')`)
	if err != nil {
		t.Fatal(err)
	}

	// gnomAD variant
	_, err = db.Exec(`INSERT INTO genomic_annotations (chrom, pos, ref, alt, gnomad_af, gnomad_ac, gnomad_an, gnomad_nhomalt, gnomad_version) VALUES ('2', 47403, 'G', 'A', '0.00123', '5', '4060', '0', '4.1')`)
	if err != nil {
		t.Fatal(err)
	}

	// Combined AM + gnomAD variant
	_, err = db.Exec(`INSERT INTO genomic_annotations (chrom, pos, ref, alt, am_score, am_class, gnomad_af, gnomad_ac, gnomad_an, gnomad_nhomalt, gnomad_version) VALUES ('3', 178936091, 'G', 'A', 0.8765, 'likely_pathogenic', '0.0001', '2', '20000', '0', '4.1')`)
	if err != nil {
		t.Fatal(err)
	}

	return dbPath
}

func TestOpen(t *testing.T) {
	dbPath := setupTestDB(t)
	store, err := Open(dbPath)
	if err != nil {
		t.Fatal(err)
	}
	defer store.Close()
}

func TestLookup(t *testing.T) {
	dbPath := setupTestDB(t)
	store, err := Open(dbPath)
	if err != nil {
		t.Fatal(err)
	}
	defer store.Close()

	tests := []struct {
		name  string
		chrom string
		pos   int64
		ref   string
		alt   string
		want  Result
		found bool
	}{
		{
			name: "AM hit", chrom: "12", pos: 25245350, ref: "C", alt: "A",
			want:  Result{AMScore: 0.9876, AMClass: "likely_pathogenic"},
			found: true,
		},
		{
			name: "ClinVar hit", chrom: "1", pos: 12345, ref: "A", alt: "G",
			want:  Result{CVClnSig: "Pathogenic", CVClnRevStat: "reviewed_by_expert_panel", CVClnDN: "Some_disease"},
			found: true,
		},
		{
			name: "combined hit", chrom: "7", pos: 140753336, ref: "A", alt: "T",
			want: Result{
				AMScore: 0.5432, AMClass: "ambiguous",
				CVClnSig: "Uncertain_significance", CVClnRevStat: "criteria_provided,_single_submitter",
			},
			found: true,
		},
		{
			name: "SIGNAL hit", chrom: "17", pos: 7577538, ref: "C", alt: "T",
			want:  Result{SigMutStatus: "germline", SigCount: "5", SigFreq: "0.00123"},
			found: true,
		},
		// Indels stored in canonical form, looked up with canonical keys.
		{
			name: "deletion canonical", chrom: "1", pos: 66927, ref: "G", alt: "",
			want: Result{CVClnSig: "Pathogenic"}, found: true,
		},
		{
			name: "multi-base deletion", chrom: "1", pos: 930329, ref: "AA", alt: "",
			want: Result{CVClnSig: "Likely_pathogenic"}, found: true,
		},
		{
			name: "insertion canonical", chrom: "1", pos: 809284, ref: "", alt: "GGTCAATCA",
			want: Result{CVClnSig: "Uncertain_significance"}, found: true,
		},
		{
			name: "1-base insertion", chrom: "1", pos: 935908, ref: "", alt: "G",
			want: Result{CVClnSig: "Benign"}, found: true,
		},
		{
			name: "gnomAD hit", chrom: "2", pos: 47403, ref: "G", alt: "A",
			want:  Result{GnomadAF: "0.00123", GnomadAC: "5", GnomadAN: "4060", GnomadNhomalt: "0", GnomadVersion: "4.1"},
			found: true,
		},
		{
			name: "combined AM + gnomAD", chrom: "3", pos: 178936091, ref: "G", alt: "A",
			want: Result{
				AMScore: 0.8765, AMClass: "likely_pathogenic",
				GnomadAF: "0.0001", GnomadAC: "2", GnomadAN: "20000", GnomadNhomalt: "0", GnomadVersion: "4.1",
			},
			found: true,
		},
		{
			name: "miss", chrom: "1", pos: 99999, ref: "A", alt: "G",
			found: false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, ok := store.Lookup(tt.chrom, tt.pos, tt.ref, tt.alt)
			if ok != tt.found {
				t.Fatalf("found=%v, want %v", ok, tt.found)
			}
			if !ok {
				return
			}
			if got.AMClass != tt.want.AMClass {
				t.Errorf("AMClass=%q, want %q", got.AMClass, tt.want.AMClass)
			}
			if got.CVClnSig != tt.want.CVClnSig {
				t.Errorf("CVClnSig=%q, want %q", got.CVClnSig, tt.want.CVClnSig)
			}
			if got.CVClnRevStat != tt.want.CVClnRevStat {
				t.Errorf("CVClnRevStat=%q, want %q", got.CVClnRevStat, tt.want.CVClnRevStat)
			}
			if got.CVClnDN != tt.want.CVClnDN {
				t.Errorf("CVClnDN=%q, want %q", got.CVClnDN, tt.want.CVClnDN)
			}
			if got.SigMutStatus != tt.want.SigMutStatus {
				t.Errorf("SigMutStatus=%q, want %q", got.SigMutStatus, tt.want.SigMutStatus)
			}
			if got.SigCount != tt.want.SigCount {
				t.Errorf("SigCount=%q, want %q", got.SigCount, tt.want.SigCount)
			}
			if got.SigFreq != tt.want.SigFreq {
				t.Errorf("SigFreq=%q, want %q", got.SigFreq, tt.want.SigFreq)
			}
			if got.GnomadAF != tt.want.GnomadAF {
				t.Errorf("GnomadAF=%q, want %q", got.GnomadAF, tt.want.GnomadAF)
			}
			if got.GnomadAC != tt.want.GnomadAC {
				t.Errorf("GnomadAC=%q, want %q", got.GnomadAC, tt.want.GnomadAC)
			}
			if got.GnomadAN != tt.want.GnomadAN {
				t.Errorf("GnomadAN=%q, want %q", got.GnomadAN, tt.want.GnomadAN)
			}
			if got.GnomadNhomalt != tt.want.GnomadNhomalt {
				t.Errorf("GnomadNhomalt=%q, want %q", got.GnomadNhomalt, tt.want.GnomadNhomalt)
			}
			if got.GnomadVersion != tt.want.GnomadVersion {
				t.Errorf("GnomadVersion=%q, want %q", got.GnomadVersion, tt.want.GnomadVersion)
			}
			// Float comparison with tolerance
			diff := got.AMScore - tt.want.AMScore
			if diff < -0.001 || diff > 0.001 {
				t.Errorf("AMScore=%v, want %v", got.AMScore, tt.want.AMScore)
			}
		})
	}
}

func TestNormalizeAlleles(t *testing.T) {
	tests := []struct {
		name            string
		pos             int64
		ref, alt        string
		wantPos         int64
		wantRef, wantAlt string
	}{
		// SNPs — no change.
		{"snp", 100, "A", "T", 100, "A", "T"},

		// Already canonical (MAF-style) — no-op.
		{"maf deletion", 100, "AT", "", 100, "AT", ""},
		{"maf insertion", 100, "", "CG", 100, "", "CG"},

		// VCF deletion: strip anchor base, advance pos.
		// POS=99 REF=GAT ALT=G → 100, "AT", ""
		{"vcf deletion", 99, "GAT", "G", 100, "AT", ""},

		// VCF 1-base deletion.
		// POS=66926 REF=AG ALT=A → 66927, "G", ""
		{"vcf 1bp deletion", 66926, "AG", "A", 66927, "G", ""},

		// VCF insertion: strip anchor, pos stays at anchor (MAF convention).
		// POS=100 REF=A ALT=ACG → 100, "", "CG"
		{"vcf insertion", 100, "A", "ACG", 100, "", "CG"},

		// VCF 1-base insertion.
		// POS=935908 REF=A ALT=AG → 935908, "", "G"
		{"vcf 1bp insertion", 935908, "A", "AG", 935908, "", "G"},

		// VCF large insertion.
		// POS=809284 REF=T ALT=TGGTCAATCA → 809284, "", "GGTCAATCA"
		{"vcf large insertion", 809284, "T", "TGGTCAATCA", 809284, "", "GGTCAATCA"},

		// VCF deletion with shared suffix.
		// POS=99 REF=GATC ALT=GC → trim prefix G (pos=100, ATC, C) → trim suffix C → 100, "AT", ""
		{"vcf del with suffix", 99, "GATC", "GC", 100, "AT", ""},

		// VCF insertion with shared suffix.
		// POS=99 REF=GT ALT=GCGT → trim prefix G (pos=100, T, CGT) → trim suffix T → insertion: pos=99, "", "CG"
		{"vcf ins with suffix", 99, "GT", "GCGT", 99, "", "CG"},

		// Complex (MNV/delins) — no common prefix.
		{"complex no prefix", 99, "GA", "TCC", 99, "GA", "TCC"},

		// Complex with shared prefix.
		// POS=99 REF=GCAT ALT=GTA → trim G → 100, "CAT", "TA"
		{"complex with prefix", 99, "GCAT", "GTA", 100, "CAT", "TA"},

		// MNV same length, no shared prefix.
		{"mnv", 100, "CT", "TC", 100, "CT", "TC"},

		// Single base, no change.
		{"identity ref", 100, "A", "A", 101, "", ""},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotPos, gotRef, gotAlt := NormalizeAlleles(tt.pos, tt.ref, tt.alt)
			if gotPos != tt.wantPos || gotRef != tt.wantRef || gotAlt != tt.wantAlt {
				t.Errorf("NormalizeAlleles(%d, %q, %q) = (%d, %q, %q), want (%d, %q, %q)",
					tt.pos, tt.ref, tt.alt,
					gotPos, gotRef, gotAlt,
					tt.wantPos, tt.wantRef, tt.wantAlt)
			}
		})
	}
}

// TestLookupVCFIndels verifies that VCF-style indels (with anchor base) can be
// normalized at lookup time to match canonical (MAF-style) entries in the DB.
func TestLookupVCFIndels(t *testing.T) {
	dbPath := setupTestDB(t)
	store, err := Open(dbPath)
	if err != nil {
		t.Fatal(err)
	}
	defer store.Close()

	tests := []struct {
		name  string
		chrom string
		pos   int64
		ref   string
		alt   string
		want  string // expected CVClnSig
	}{
		// VCF deletion: POS=66926 REF=AG ALT=A → canonical (66927, G, "")
		{"vcf deletion", "1", 66926, "AG", "A", "Pathogenic"},
		// VCF multi-base deletion: POS=930328 REF=GAA ALT=G → canonical (930329, AA, "")
		{"vcf multi-del", "1", 930328, "GAA", "G", "Likely_pathogenic"},
		// VCF insertion: POS=809284 REF=T ALT=TGGTCAATCA → canonical (809284, "", "GGTCAATCA")
		{"vcf insertion", "1", 809284, "T", "TGGTCAATCA", "Uncertain_significance"},
		// VCF 1bp insertion: POS=935908 REF=A ALT=AG → canonical (935908, "", "G")
		{"vcf 1bp insertion", "1", 935908, "A", "AG", "Benign"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Normalize the VCF-style variant before lookup (same as GenomicSource.Annotate does).
			nPos, nRef, nAlt := NormalizeAlleles(tt.pos, tt.ref, tt.alt)
			r, ok := store.Lookup(tt.chrom, nPos, nRef, nAlt)
			if !ok {
				t.Fatalf("expected hit for VCF variant (%s, %d, %s, %s) → normalized (%d, %s, %s)",
					tt.chrom, tt.pos, tt.ref, tt.alt, nPos, nRef, nAlt)
			}
			if r.CVClnSig != tt.want {
				t.Errorf("CVClnSig=%q, want %q", r.CVClnSig, tt.want)
			}
		})
	}
}

// TestLookupMAFIndels verifies that MAF-style indels (no anchor base) match
// canonical entries directly — normalization is a no-op for these.
func TestLookupMAFIndels(t *testing.T) {
	dbPath := setupTestDB(t)
	store, err := Open(dbPath)
	if err != nil {
		t.Fatal(err)
	}
	defer store.Close()

	tests := []struct {
		name  string
		chrom string
		pos   int64
		ref   string
		alt   string
		want  string
	}{
		// MAF deletion: Start=66927 Ref=G Alt="" (already canonical)
		{"maf deletion", "1", 66927, "G", "", "Pathogenic"},
		// MAF multi-base deletion: Start=930329 Ref=AA Alt=""
		{"maf multi-del", "1", 930329, "AA", "", "Likely_pathogenic"},
		// MAF insertion: Start=809284 Ref="" Alt="GGTCAATCA"
		{"maf insertion", "1", 809284, "", "GGTCAATCA", "Uncertain_significance"},
		// MAF 1bp insertion: Start=935908 Ref="" Alt="G"
		{"maf 1bp ins", "1", 935908, "", "G", "Benign"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Normalize should be a no-op for MAF-style.
			nPos, nRef, nAlt := NormalizeAlleles(tt.pos, tt.ref, tt.alt)
			if nPos != tt.pos || nRef != tt.ref || nAlt != tt.alt {
				t.Errorf("NormalizeAlleles changed MAF variant: (%d,%q,%q) → (%d,%q,%q)",
					tt.pos, tt.ref, tt.alt, nPos, nRef, nAlt)
			}

			r, ok := store.Lookup(tt.chrom, nPos, nRef, nAlt)
			if !ok {
				t.Fatalf("expected hit for MAF variant (%s, %d, %q, %q)", tt.chrom, tt.pos, tt.ref, tt.alt)
			}
			if r.CVClnSig != tt.want {
				t.Errorf("CVClnSig=%q, want %q", r.CVClnSig, tt.want)
			}
		})
	}
}

// TestLoadGnomad verifies that loadGnomad correctly parses a VCF file,
// normalizes indels, and inserts into the database.
func TestLoadGnomad(t *testing.T) {
	dir := t.TempDir()
	dbPath := filepath.Join(dir, "test.sqlite")

	db, err := sql.Open("sqlite", dbPath)
	if err != nil {
		t.Fatal(err)
	}
	defer db.Close()

	_, err = db.Exec(`CREATE TABLE genomic_annotations (
		chrom TEXT NOT NULL, pos INTEGER NOT NULL, ref TEXT NOT NULL, alt TEXT NOT NULL,
		am_score REAL NOT NULL DEFAULT 0, am_class TEXT NOT NULL DEFAULT '',
		cv_clnsig TEXT NOT NULL DEFAULT '', cv_revstat TEXT NOT NULL DEFAULT '', cv_clndn TEXT NOT NULL DEFAULT '',
		sig_mut_status TEXT NOT NULL DEFAULT '', sig_count TEXT NOT NULL DEFAULT '', sig_freq TEXT NOT NULL DEFAULT '',
		gnomad_af TEXT NOT NULL DEFAULT '', gnomad_ac TEXT NOT NULL DEFAULT '', gnomad_an TEXT NOT NULL DEFAULT '',
		gnomad_nhomalt TEXT NOT NULL DEFAULT '', gnomad_version TEXT NOT NULL DEFAULT '',
		PRIMARY KEY (chrom, pos, ref, alt)
	) WITHOUT ROWID`)
	if err != nil {
		t.Fatal(err)
	}

	// Write a small test VCF file (uncompressed).
	vcfContent := `##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer>
##INFO=<ID=AN,Number=1,Type=Integer>
##INFO=<ID=AF,Number=A,Type=Float>
##INFO=<ID=nhomalt,Number=A,Type=Integer>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10042	.	T	A	.	PASS	AC=4;AN=1526;AF=0.002621;nhomalt=0
chr1	10055	.	T	C	.	AC0	AC=0;AN=0;AF=0;nhomalt=0
chr2	20000	.	GA	G	.	PASS	AC=10;AN=5000;AF=0.002;nhomalt=1
chr3	30000	.	T	TCGA	.	PASS	AC=3;AN=10000;AF=0.0003;nhomalt=0
`
	vcfPath := filepath.Join(dir, "test_gnomad.vcf")
	if err := os.WriteFile(vcfPath, []byte(vcfContent), 0644); err != nil {
		t.Fatal(err)
	}

	n, err := loadGnomad(db, vcfPath, "4.1")
	if err != nil {
		t.Fatal(err)
	}

	// Should load 3 variants (1 SNP, 1 deletion, 1 insertion; the AC0 variant is filtered).
	if n != 3 {
		t.Fatalf("loaded %d variants, want 3", n)
	}

	// Verify the SNP.
	store, err := Open(dbPath)
	if err != nil {
		t.Fatal(err)
	}
	defer store.Close()

	r, ok := store.Lookup("1", 10042, "T", "A")
	if !ok {
		t.Fatal("expected SNP hit at chr1:10042")
	}
	if r.GnomadAF != "0.002621" {
		t.Errorf("GnomadAF=%q, want %q", r.GnomadAF, "0.002621")
	}
	if r.GnomadAC != "4" {
		t.Errorf("GnomadAC=%q, want %q", r.GnomadAC, "4")
	}
	if r.GnomadAN != "1526" {
		t.Errorf("GnomadAN=%q, want %q", r.GnomadAN, "1526")
	}
	if r.GnomadVersion != "4.1" {
		t.Errorf("GnomadVersion=%q, want %q", r.GnomadVersion, "4.1")
	}

	// Verify the deletion: VCF GA→G at pos 20000 → normalized (20001, "A", "").
	r, ok = store.Lookup("2", 20001, "A", "")
	if !ok {
		t.Fatal("expected deletion hit at chr2:20001 (normalized from VCF pos 20000 GA>G)")
	}
	if r.GnomadAF != "0.002" {
		t.Errorf("GnomadAF=%q, want %q", r.GnomadAF, "0.002")
	}
	if r.GnomadNhomalt != "1" {
		t.Errorf("GnomadNhomalt=%q, want %q", r.GnomadNhomalt, "1")
	}

	// Verify the insertion: VCF T→TCGA at pos 30000 → normalized (30000, "", "CGA").
	r, ok = store.Lookup("3", 30000, "", "CGA")
	if !ok {
		t.Fatal("expected insertion hit at chr3:30000 (normalized from VCF pos 30000 T>TCGA)")
	}
	if r.GnomadAF != "0.0003" {
		t.Errorf("GnomadAF=%q, want %q", r.GnomadAF, "0.0003")
	}
}

// TestGenomicSourceAnnotateGnomad verifies that GenomicSource.Annotate correctly
// distributes gnomAD data to annotations via the Extra map.
func TestGenomicSourceAnnotateGnomad(t *testing.T) {
	dbPath := setupTestDB(t)
	store, err := Open(dbPath)
	if err != nil {
		t.Fatal(err)
	}
	defer store.Close()

	src := NewSource(store, "1.0")

	// Verify columns include gnomAD.
	cols := src.Columns()
	gnomadCols := 0
	for _, c := range cols {
		if len(c.Name) > 7 && c.Name[:7] == "gnomad." {
			gnomadCols++
		}
	}
	if gnomadCols != 5 {
		t.Errorf("expected 5 gnomAD columns, got %d", gnomadCols)
	}

	// Test gnomAD-only variant (chr2:47403 G>A).
	v := &vcf.Variant{Chrom: "2", Pos: 47403, Ref: "G", Alt: "A"}
	anns := []*annotate.Annotation{
		{Consequence: "missense_variant"},
		{Consequence: "synonymous_variant"},
	}
	src.Annotate(v, anns)

	// gnomAD should be on all annotations.
	for i, ann := range anns {
		if got := ann.GetExtraKey("gnomad.af"); got != "0.00123" {
			t.Errorf("ann[%d] gnomad.af=%q, want %q", i, got, "0.00123")
		}
		if got := ann.GetExtraKey("gnomad.ac"); got != "5" {
			t.Errorf("ann[%d] gnomad.ac=%q, want %q", i, got, "5")
		}
		if got := ann.GetExtraKey("gnomad.an"); got != "4060" {
			t.Errorf("ann[%d] gnomad.an=%q, want %q", i, got, "4060")
		}
		if got := ann.GetExtraKey("gnomad.nhomalt"); got != "0" {
			t.Errorf("ann[%d] gnomad.nhomalt=%q, want %q", i, got, "0")
		}
		if got := ann.GetExtraKey("gnomad.version"); got != "4.1" {
			t.Errorf("ann[%d] gnomad.version=%q, want %q", i, got, "4.1")
		}
	}

	// Test combined AM + gnomAD variant (chr3:178936091 G>A).
	// AlphaMissense should only apply to missense, but gnomAD applies to all.
	v2 := &vcf.Variant{Chrom: "3", Pos: 178936091, Ref: "G", Alt: "A"}
	anns2 := []*annotate.Annotation{
		{Consequence: "missense_variant"},
		{Consequence: "synonymous_variant"},
	}
	src.Annotate(v2, anns2)

	// Missense: should have both AM and gnomAD.
	if got := anns2[0].GetExtraKey("alphamissense.score"); got == "" {
		t.Error("missense annotation should have alphamissense.score")
	}
	if got := anns2[0].GetExtraKey("gnomad.af"); got != "0.0001" {
		t.Errorf("missense gnomad.af=%q, want %q", got, "0.0001")
	}

	// Synonymous: should have gnomAD but not AM.
	if got := anns2[1].GetExtraKey("alphamissense.score"); got != "" {
		t.Errorf("synonymous should not have alphamissense.score, got %q", got)
	}
	if got := anns2[1].GetExtraKey("gnomad.af"); got != "0.0001" {
		t.Errorf("synonymous gnomad.af=%q, want %q", got, "0.0001")
	}

	// Test miss: no annotations set.
	v3 := &vcf.Variant{Chrom: "99", Pos: 1, Ref: "A", Alt: "G"}
	anns3 := []*annotate.Annotation{{Consequence: "missense_variant"}}
	src.Annotate(v3, anns3)
	if got := anns3[0].GetExtraKey("gnomad.af"); got != "" {
		t.Errorf("miss should have empty gnomad.af, got %q", got)
	}
}

func TestReady(t *testing.T) {
	dir := t.TempDir()
	dbPath := filepath.Join(dir, "test.sqlite")

	sources := BuildSources{AlphaMissenseTSV: "/nonexistent"}

	// DB doesn't exist yet.
	if Ready(dbPath, sources) {
		t.Error("Ready should return false when DB doesn't exist")
	}

	// Empty file — not a valid SQLite database.
	f, err := os.Create(dbPath)
	if err != nil {
		t.Fatal(err)
	}
	f.Close()
	if Ready(dbPath, sources) {
		t.Error("Ready should return false for empty file")
	}

	// Corrupted file (non-SQLite content).
	if err := os.WriteFile(dbPath, []byte("not a sqlite database"), 0644); err != nil {
		t.Fatal(err)
	}
	if Ready(dbPath, sources) {
		t.Error("Ready should return false for corrupted file")
	}

	// Valid SQLite but wrong schema (missing genomic_annotations table).
	os.Remove(dbPath)
	wrongDB, err := sql.Open("sqlite", dbPath)
	if err != nil {
		t.Fatal(err)
	}
	wrongDB.Exec("CREATE TABLE wrong_table (id INTEGER)")
	wrongDB.Close()
	if Ready(dbPath, sources) {
		t.Error("Ready should return false for wrong schema")
	}

	// Valid database — use setupTestDB.
	os.Remove(dbPath)
	goodPath := setupTestDB(t)
	goodSources := BuildSources{AlphaMissenseTSV: "/nonexistent"}
	if !Ready(goodPath, goodSources) {
		t.Error("Ready should return true for valid database")
	}
}
