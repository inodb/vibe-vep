package output

import (
	"bytes"
	"encoding/json"
	"strings"
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

func krasVariantAndAnnotation() (*vcf.Variant, *annotate.Annotation) {
	v := &vcf.Variant{Chrom: "12", Pos: 25398284, Ref: "C", Alt: "T"}
	ann := &annotate.Annotation{
		TranscriptID:    "ENST00000311936.3",
		GeneName:        "KRAS",
		GeneID:          "ENSG00000133703",
		Consequence:     "missense_variant",
		Impact:          "MODERATE",
		Allele:          "T",
		Biotype:         "protein_coding",
		ProteinPosition: 12,
		AminoAcidChange: "G12D",
		CodonChange:     "ggt/gAt",
		HGVSc:           "c.35G>A",
		HGVSp:           "p.Gly12Asp",
		RefSeqIDs:       []string{"NM_004985.3", "NM_033360.4"},
	}
	return v, ann
}

func vepLineFromWriter(t *testing.T, v *vcf.Variant, ann *annotate.Annotation) string {
	t.Helper()
	var buf bytes.Buffer
	w := NewJSONLWriter(&buf, "ensembl-vep-jsonl", "GRCh37")
	w.SetInput("12,25398284,25398284,C,T")
	if err := w.Write(v, ann); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}
	return buf.String()
}

// genome-nexus's RefSeqResolver reads refseq_transcript_ids and reports the
// first element, so both presence and order matter.
func TestVEPEmitsRefSeqTranscriptIDs(t *testing.T) {
	v, ann := krasVariantAndAnnotation()

	var result VEPVariantAnnotation
	if err := json.Unmarshal([]byte(vepLineFromWriter(t, v, ann)), &result); err != nil {
		t.Fatalf("invalid JSON: %v", err)
	}
	if len(result.TranscriptConsequences) != 1 {
		t.Fatalf("got %d transcript consequences, want 1", len(result.TranscriptConsequences))
	}
	got := result.TranscriptConsequences[0].RefSeqTranscriptIDs
	want := []string{"NM_004985.3", "NM_033360.4"}
	if len(got) != len(want) {
		t.Fatalf("refseq_transcript_ids=%v, want %v", got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("refseq_transcript_ids[%d]=%q, want %q", i, got[i], want[i])
		}
	}
}

// A nil slice must omit the key entirely rather than emit null/[], so the
// output stays shaped like real VEP JSON.
func TestVEPOmitsRefSeqWhenAbsent(t *testing.T) {
	v, ann := krasVariantAndAnnotation()
	ann.RefSeqIDs = nil

	line := vepLineFromWriter(t, v, ann)
	if strings.Contains(line, "refseq_transcript_ids") {
		t.Errorf("expected refseq_transcript_ids to be omitted, got: %s", line)
	}
}

// marshalVEP (annotate stream) and MarshalVEPAnnotation (serve) are
// hand-maintained duplicates. This pins them together so RefSeq cannot be
// added to one and forgotten in the other.
func TestVEPEmittersAgreeOnRefSeq(t *testing.T) {
	v, ann := krasVariantAndAnnotation()

	standalone, err := MarshalVEPAnnotation("12,25398284,25398284,C,T", v, []*annotate.Annotation{ann}, "GRCh37")
	if err != nil {
		t.Fatal(err)
	}

	var fromWriter, fromStandalone VEPVariantAnnotation
	if err := json.Unmarshal([]byte(vepLineFromWriter(t, v, ann)), &fromWriter); err != nil {
		t.Fatal(err)
	}
	if err := json.Unmarshal(standalone, &fromStandalone); err != nil {
		t.Fatal(err)
	}

	a := fromWriter.TranscriptConsequences[0].RefSeqTranscriptIDs
	b := fromStandalone.TranscriptConsequences[0].RefSeqTranscriptIDs
	if len(a) == 0 {
		t.Fatal("JSONLWriter emitted no refseq_transcript_ids")
	}
	if len(a) != len(b) {
		t.Fatalf("emitters disagree: writer=%v standalone=%v", a, b)
	}
	for i := range a {
		if a[i] != b[i] {
			t.Errorf("emitters disagree at [%d]: writer=%q standalone=%q", i, a[i], b[i])
		}
	}
}

// The genome-nexus response shape carries the full list on the consequence and
// the first accession on the summary, mirroring RefSeqResolver.
func TestGNResponseCarriesRefSeq(t *testing.T) {
	v, ann := krasVariantAndAnnotation()

	data, err := MarshalGNAnnotation("12,25398284,25398284,C,T", v, []*annotate.Annotation{ann}, "GRCh37")
	if err != nil {
		t.Fatal(err)
	}

	var gn GNAnnotation
	if err := json.Unmarshal(data, &gn); err != nil {
		t.Fatalf("invalid JSON: %v\n%s", err, data)
	}
	if len(gn.TranscriptConsequences) == 0 {
		t.Fatal("no transcript consequences")
	}
	if got := gn.TranscriptConsequences[0].RefseqTranscriptIds; len(got) != 2 || got[0] != "NM_004985.3" {
		t.Errorf("refseq_transcript_ids=%v, want [NM_004985.3 NM_033360.4]", got)
	}
}
