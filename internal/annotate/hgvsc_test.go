package annotate

import (
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

func TestFormatHGVSc_KRAS_SNV(t *testing.T) {
	// KRAS G12C: c.34G>T on reverse strand
	// Genomic pos 25245351, ref=C (genomic), alt=A (genomic)
	// On coding strand: ref=G, alt=T → c.34G>T
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245351,
		Ref:   "C",
		Alt:   "A",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.34G>T", hgvsc)
}

func TestFormatHGVSc_KRAS_Intronic(t *testing.T) {
	v := &vcf.Variant{
		Chrom: "12",
		Pos:   25245272, // 2bp before exon 2 start (intronic, on the 3' side for reverse strand)
		Ref:   "A",
		Alt:   "G",
	}

	transcript := createKRASTranscript()
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("HGVSc for intronic position: %s (consequence: %s)", hgvsc, result.Consequence)
	assert.NotEmpty(t, hgvsc)
}

func TestFormatHGVSc_ForwardStrand_SNV(t *testing.T) {
	transcript := createForwardTranscript()

	v := &vcf.Variant{
		Chrom: "1",
		Pos:   1005,
		Ref:   "A",
		Alt:   "G",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("Forward strand SNV HGVSc: %s", hgvsc)
	assert.NotEmpty(t, hgvsc)
	if len(hgvsc) > 2 {
		assert.True(t, strings.HasPrefix(hgvsc, "c."), "expected HGVSc to start with 'c.', got %q", hgvsc)
	}
}

func TestFormatHGVSc_Deletion(t *testing.T) {
	transcript := createForwardTranscript()

	v := &vcf.Variant{
		Chrom: "1",
		Pos:   1005,
		Ref:   "AG",
		Alt:   "A",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("Deletion HGVSc: %s", hgvsc)
	assert.NotEmpty(t, hgvsc)
	assert.True(t, strings.HasSuffix(hgvsc, "del"), "expected HGVSc to end with 'del', got %q", hgvsc)
}

func TestFormatHGVSc_Insertion(t *testing.T) {
	transcript := createForwardTranscript()

	v := &vcf.Variant{
		Chrom: "1",
		Pos:   1005,
		Ref:   "A",
		Alt:   "AGG",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("Insertion HGVSc: %s", hgvsc)
	assert.NotEmpty(t, hgvsc)
}

func TestFormatHGVSc_UTR(t *testing.T) {
	transcript := createForwardTranscript()

	v := &vcf.Variant{
		Chrom: "1",
		Pos:   998, // In exon 1 but before CDSStart=1000
		Ref:   "A",
		Alt:   "G",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	t.Logf("5'UTR HGVSc: %s (consequence: %s)", hgvsc, result.Consequence)
	assert.NotEmpty(t, hgvsc)
	if len(hgvsc) > 2 {
		assert.True(t, strings.HasPrefix(hgvsc, "c.-"), "expected HGVSc to start with 'c.-', got %q", hgvsc)
	}
}

func TestFormatHGVSc_UpstreamDownstream(t *testing.T) {
	transcript := createForwardTranscript()

	v := &vcf.Variant{
		Chrom: "1",
		Pos:   500,
		Ref:   "A",
		Alt:   "G",
	}

	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Empty(t, hgvsc)
}

func TestGenomicToHGVScPos_CDS(t *testing.T) {
	transcript := createForwardTranscript()

	pos := genomicToHGVScPos(1005, transcript)
	assert.Equal(t, "6", pos)
}

func TestGenomicToHGVScPos_FivePrimeUTR(t *testing.T) {
	transcript := createForwardTranscript()

	pos := genomicToHGVScPos(998, transcript)
	assert.Equal(t, "-2", pos)
}

func TestFormatHGVSc_DupPrecedingBase(t *testing.T) {
	// Test: insertion that duplicates the preceding (anchor) base.
	// Forward strand transcript, CDS = "ATGATCGATCG..."
	// Insert at pos 1002 (CDS pos 3 = 'G'), ref=G, alt=GG
	// Inserted 'G' matches CDS pos 3 → c.3dup
	transcript := createDupTestTranscript()

	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "G", Alt: "GG"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.3dup", hgvsc)
}

func TestFormatHGVSc_DupFollowingBase(t *testing.T) {
	// Test: insertion that duplicates the following base (not the anchor).
	// Forward strand transcript, CDS = "ATGATCGATCG..."
	// CDS pos 6='C', pos 7='G'
	// Insert at pos 1005 (CDS pos 6 = 'C'), ref=C, alt=CG
	// Inserted 'G' doesn't match CDS pos 6 ('C'), but matches CDS pos 7 ('G') → c.7dup
	transcript := createDupTestTranscript()

	v := &vcf.Variant{Chrom: "1", Pos: 1005, Ref: "C", Alt: "CG"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.7dup", hgvsc)
}

func TestFormatHGVSc_DupMultiBase(t *testing.T) {
	// Test: multi-base duplication.
	// CDS = "ATGATCGATCG..."
	// CDS pos 1-3 = "ATG"
	// Insert at pos 1002 (CDS 3), ref=G, alt=GATG
	// Inserted "ATG" matches CDS pos 1-3, but 3' shifting moves it rightward
	// through matching bases (A→T→C partial matches via rotation) → c.3_5dup
	transcript := createDupTestTranscript()

	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "G", Alt: "GATG"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.3_5dup", hgvsc)
}

func TestFormatHGVSc_InsertionNotDup(t *testing.T) {
	// Test: insertion that does NOT match adjacent bases → plain insertion.
	// CDS = "ATGATCGATCG..."
	// Insert at pos 1002 (CDS 3 = 'G'), ref=G, alt=GCC
	// Inserted 'CC' doesn't match "TG" (preceding) or "AT" (following) → insertion
	transcript := createDupTestTranscript()

	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "G", Alt: "GCC"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Contains(t, hgvsc, "ins")
	assert.NotContains(t, hgvsc, "dup")
}

func TestFormatHGVSc_DupReverseStrand(t *testing.T) {
	// Test: duplication on reverse strand.
	// KRAS reverse strand. CDS pos 34 = 'G', CDS pos 35 = 'G' (GGT → Gly12)
	// Insert 'G' at CDS pos 34: with 3' shift, shifts right to pos 35
	// (both positions are G), then stops at pos 36='T' → c.35dup
	transcript := createKRASTranscript()

	// VCF: pos=25245351, ref=C, alt=CC → on coding strand: ref=G, alt=GG → ins G
	v := &vcf.Variant{Chrom: "12", Pos: 25245351, Ref: "C", Alt: "CC"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.35dup", hgvsc)
}

func TestShiftDeletionThreePrime(t *testing.T) {
	tests := []struct {
		name     string
		cdsSeq   string
		delStart int
		delEnd   int
		wantS    int
		wantE    int
	}{
		{
			name:     "poly-A shift",
			cdsSeq:   "ATGAAAACCC",
			delStart: 3, delEnd: 5, // delete "AAA" at positions 3-5
			wantS: 4, wantE: 6, // shifts right by 1 (next base is also A)
		},
		{
			name:     "no shift - next base differs",
			cdsSeq:   "ATGCCCAAA",
			delStart: 3, delEnd: 5, // delete "CCC" at positions 3-5
			wantS: 3, wantE: 5, // C≠A, no shift
		},
		{
			name:     "single base shift",
			cdsSeq:   "ATGAACCC",
			delStart: 3, delEnd: 3, // delete single "A" at position 3
			wantS: 4, wantE: 4, // shifts right by 1 (next base is also A)
		},
		{
			name:     "shift to end of sequence",
			cdsSeq:   "ATGAAA",
			delStart: 3, delEnd: 4, // delete "AA" at positions 3-4
			wantS: 4, wantE: 5, // shifts right by 1 (position 5 is A), then stops at end
		},
		{
			name:     "at end of sequence - no shift",
			cdsSeq:   "ATGCCC",
			delStart: 3, delEnd: 5, // delete last 3 bases
			wantS: 3, wantE: 5, // can't shift past end
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotS, gotE := shiftDeletionThreePrime(tt.delStart, tt.delEnd, tt.cdsSeq)
			assert.Equal(t, tt.wantS, gotS, "delStart")
			assert.Equal(t, tt.wantE, gotE, "delEnd")
		})
	}
}

func TestShiftInsertionThreePrime(t *testing.T) {
	tests := []struct {
		name       string
		cdsSeq     string
		insertSeq  string
		anchorIdx  int
		wantSeq    string
		wantIdx    int
	}{
		{
			name:      "shift through poly-A",
			cdsSeq:    "ATGAAACCC",
			insertSeq: "A", anchorIdx: 3,
			wantSeq: "A", wantIdx: 5, // shifts through AAA run, stops at C
		},
		{
			name:      "no shift - next base differs",
			cdsSeq:    "ATGAAACCC",
			insertSeq: "G", anchorIdx: 3,
			wantSeq: "G", wantIdx: 3, // A≠G, no shift
		},
		{
			name:      "multi-base rotation",
			cdsSeq:    "ATGATGATGCCC",
			insertSeq: "ATG", anchorIdx: 2,
			wantSeq: "ATG", wantIdx: 8, // shifts through all three ATG repeats
		},
		{
			name:      "shift single base to end",
			cdsSeq:    "ATGAAA",
			insertSeq: "A", anchorIdx: 3,
			wantSeq: "A", wantIdx: 5, // shifts to last position
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotSeq, gotIdx := shiftInsertionThreePrime(tt.insertSeq, tt.anchorIdx, tt.cdsSeq)
			assert.Equal(t, tt.wantSeq, gotSeq, "shifted sequence")
			assert.Equal(t, tt.wantIdx, gotIdx, "shifted anchor index")
		})
	}
}

func TestFormatHGVSc_DeletionShiftPolyA(t *testing.T) {
	// Forward strand, deletion in poly-A tract should shift right.
	// CDS: "ATGAAAACCC..." with A's at positions 4-7 (1-based)
	// Delete first A (genomic pos 1003, ref=AA, alt=A) → VCF deletes position 1004 (CDS pos 5)
	// Without shift: c.5del, with shift: c.7del (shifted to last A in run)
	cds := "ATGAAAACCCATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGA"
	transcript := &cache.Transcript{
		ID: "ENST_DEL_SHIFT", GeneName: "DELSHIFT", Chrom: "1",
		Start: 990, End: 1210, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	v := &vcf.Variant{Chrom: "1", Pos: 1003, Ref: "AA", Alt: "A"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.7del", hgvsc)
}

func TestFormatHGVSc_DeletionShiftMultiBase(t *testing.T) {
	// ATG repeat: positions 1-9 are ATGATGATG
	// Delete positions 4-6 (second ATG) → should shift to 7-9 (last ATG in repeat)
	cds := "ATGATGATGCCCATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGA"
	transcript := &cache.Transcript{
		ID: "ENST_DEL_MULTI", GeneName: "DELMULTI", Chrom: "1",
		Start: 990, End: 1210, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	// VCF: pos=1002 (CDS 3=G), ref=GATG, alt=G → deletes positions 1003-1005 (CDS 4-6: ATG)
	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "GATG", Alt: "G"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.7_9del", hgvsc)
}

func TestFormatHGVSc_InsertionBecomesDupAfterShift(t *testing.T) {
	// CDS: "ATGCCCAAAGGG..."
	// Insert "A" at pos 1005 (anchor CDS pos 6='C'), ref=C, alt=CA
	// At original position: inserted A doesn't match preceding C → not dup
	// After shifting: A matches CDS[7]='A', shifts to anchor 7, then CDS[8]='A' → anchor 8,
	// then CDS[9]='G' → stops. Anchor=8 (0-based). Check dup: CDS[8]='A' == inserted 'A' → c.9dup
	cds := "ATGCCCAAAGGGGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
	transcript := &cache.Transcript{
		ID: "ENST_INS_DUP", GeneName: "INSDUP", Chrom: "1",
		Start: 990, End: 1210, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	v := &vcf.Variant{Chrom: "1", Pos: 1005, Ref: "C", Alt: "CA"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.9dup", hgvsc)
}

func TestFormatHGVSc_DupPositionShift(t *testing.T) {
	// CDS: "ATGAAACCC..."
	// Insert "A" at pos 1003 (anchor CDS pos 4='A'), ref=A, alt=AA
	// At original position: dup of pos 4. But HGVS requires 3' shift.
	// After shifting: anchor shifts from 4 to 5 (CDS[5]='A'), then to 5 (CDS[6]='C'≠'A') → anchor=5
	// Dup check: CDS[5]='A' == inserted 'A' → c.6dup
	cds := "ATGAAACCCATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGA"
	transcript := &cache.Transcript{
		ID: "ENST_DUP_SHIFT", GeneName: "DUPSHIFT", Chrom: "1",
		Start: 990, End: 1210, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	v := &vcf.Variant{Chrom: "1", Pos: 1003, Ref: "A", Alt: "AA"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.6dup", hgvsc)
}

func TestFormatHGVSc_DeletionNoShift(t *testing.T) {
	// CDS: "ATGCCCAAA..." — delete "CCC" at positions 4-6
	// Next base after deletion is 'A' ≠ 'C' → no shift
	cds := "ATGCCCAAAGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
	transcript := &cache.Transcript{
		ID: "ENST_DEL_NOSHIFT", GeneName: "DELNOSHIFT", Chrom: "1",
		Start: 990, End: 1210, Strand: 1, Biotype: "protein_coding",
		CDSStart: 1000, CDSEnd: 1101,
		Exons:       []cache.Exon{{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0}},
		CDSSequence: cds,
	}

	// VCF: pos=1002 (CDS 3=G), ref=GCCC, alt=G → deletes CDS 4-6
	v := &vcf.Variant{Chrom: "1", Pos: 1002, Ref: "GCCC", Alt: "G"}
	result := PredictConsequence(v, transcript)
	hgvsc := FormatHGVSc(v, transcript, result)

	assert.Equal(t, "c.4_6del", hgvsc)
}

// createDupTestTranscript creates a forward-strand transcript with a known CDS for dup testing.
func createDupTestTranscript() *cache.Transcript {
	// Single exon forward strand transcript for simplicity.
	// Exon: 990-1200, CDS: 1000-1101 (102bp = 34 codons)
	// CDS: pos 1=1000, pos 2=1001, ..., pos 102=1101
	// CDSSequence[0..101] = "ATGATCGATCG..." (repeating pattern)
	cds := "ATGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
	return &cache.Transcript{
		ID:       "ENST00000DUP001",
		GeneName: "DUPTEST",
		Chrom:    "1",
		Start:    990,
		End:      1210,
		Strand:   1,
		Biotype:  "protein_coding",
		CDSStart: 1000,
		CDSEnd:   1101,
		Exons: []cache.Exon{
			{Number: 1, Start: 990, End: 1210, CDSStart: 1000, CDSEnd: 1101, Frame: 0},
		},
		CDSSequence: cds,
	}
}

// createForwardTranscript creates a simple forward-strand transcript for testing.
func createForwardTranscript() *cache.Transcript {
	// Simple 2-exon forward strand transcript
	// Exon 1: 990-1020 (CDS starts at 1000)
	// Intron: 1021-1099
	// Exon 2: 1100-1200 (CDS ends at 1180)
	//
	// CDS: 1000-1020 (21bp from exon 1) + 1100-1180 (81bp from exon 2) = 102bp = 34 codons
	return &cache.Transcript{
		ID:          "ENST00000000001",
		GeneID:      "ENSG00000000001",
		GeneName:    "TEST",
		Chrom:       "1",
		Start:       990,
		End:         1200,
		Strand:      1,
		Biotype:     "protein_coding",
		IsCanonical: true,
		CDSStart:    1000,
		CDSEnd:      1180,
		Exons: []cache.Exon{
			{Number: 1, Start: 990, End: 1020, CDSStart: 1000, CDSEnd: 1020, Frame: 0},
			{Number: 2, Start: 1100, End: 1200, CDSStart: 1100, CDSEnd: 1180, Frame: 0},
		},
		CDSSequence: "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG",
	}
}

func BenchmarkFormatHGVSc_SNV(b *testing.B) {
	transcript := createKRASTranscript()
	v := &vcf.Variant{Chrom: "12", Pos: 25245351, Ref: "C", Alt: "A"}
	result := PredictConsequence(v, transcript)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		FormatHGVSc(v, transcript, result)
	}
}

func BenchmarkFormatHGVSc_Deletion(b *testing.B) {
	transcript := createKRASTranscript()
	v := &vcf.Variant{Chrom: "12", Pos: 25245350, Ref: "GG", Alt: "G"}
	result := PredictConsequence(v, transcript)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		FormatHGVSc(v, transcript, result)
	}
}

func BenchmarkFormatHGVSc_Insertion(b *testing.B) {
	transcript := createKRASTranscript()
	v := &vcf.Variant{Chrom: "12", Pos: 25245350, Ref: "C", Alt: "CCC"}
	result := PredictConsequence(v, transcript)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		FormatHGVSc(v, transcript, result)
	}
}

// TestAllocRegression_FormatHGVSc verifies allocation count for HGVSc formatting.
func TestAllocRegression_FormatHGVSc(t *testing.T) {
	transcript := createKRASTranscript()
	snv := &vcf.Variant{Chrom: "12", Pos: 25245351, Ref: "C", Alt: "A"}
	result := PredictConsequence(snv, transcript)

	allocs := testing.AllocsPerRun(100, func() {
		FormatHGVSc(snv, transcript, result)
	})

	// SNV HGVSc: ReverseComplement (1 alloc each for ref+alt) + position string + concat.
	// Budget: ~5 allocs.
	const maxAllocs = 8
	if int(allocs) > maxAllocs {
		t.Errorf("FormatHGVSc allocation regression: got %.0f, want <= %d", allocs, maxAllocs)
	}
	t.Logf("FormatHGVSc(SNV) allocs: %.0f (budget: %d)", allocs, maxAllocs)
}
