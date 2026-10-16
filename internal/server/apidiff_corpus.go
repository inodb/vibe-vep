//go:build apidiff

package server

import (
	"encoding/json"
	"fmt"
	"strings"
)

// The corpus the differential harness drives both servers with. It is checked
// in rather than generated so that a run is reproducible and a difference can
// be chased back to a named case.
//
// Variants are chosen for coverage of the shapes that have actually gone wrong
// before, not for statistical representativeness: the MAF is 70% missense, so
// sampling it uniformly would barely exercise the cases that break.

// apiCase is one request made to both servers.
type apiCase struct {
	// Name identifies the case in failure output.
	Name string
	// Method and Path are the request; Body is sent for POST.
	Method string
	Path   string
	Body   string
	// SortKey names a field to sort top-level arrays by before comparing,
	// where the order is not meaningful.
	SortKey string
}

// hardVariants are the cases earlier work found to be decisive. Each one cost
// real time to diagnose, so each one is worth a test.
var hardVariants = []struct {
	name, genomic, why string
}{
	{"BRAF V600E", "7,140453136,140453136,A,T", "the canonical missense"},
	{"CDKN2A dual canonical", "9,21971173,21971173,A,G", "two MSKCC canonical transcripts; order decides"},
	{"CDKN2A frameshift", "9,21971091,21971092,-,C", "same, on an insertion"},
	{"RTEL1 readthrough", "20,62319056,62319056,G,A", "also in RTEL1-TNFRSF6B; the cancer gene filter decides"},
	{"APLNR NMD canonical", "11,57003943,57003943,A,G", "the MSKCC canonical is an NMD transcript"},
	{"NSD2 retired symbol", "4,1980527,1980527,C,T", "GENCODE calls it WHSC1"},
	{"KIT reVUE", "4,55593576,55593606,CCACAGAAACCCATGTATGAAGTACAGTGGA,-", "published effect overrules the prediction"},
	{"ERBB2 reVUE", "17,37876039,37876039,G,A", "same, and absent from a stale catalogue"},
	{"TP53 5'UTR splice", "17,7590693,7590693,A,G", "c.-29+2T>C anchors at residue 1"},
	{"PBRM1 frameshift", "3,52675445,52675445,T,-", "protein position is the unshifted one"},
}

// hardTranscripts exercise the lookups that do not involve a variant.
var hardTranscripts = []string{
	"ENST00000288602", // BRAF
	"ENST00000269305", // TP53
	"ENST00000311936", // KRAS
	"ENST00000477973", // RYBP, the one cds_start_NF transcript selected
	"ENST00000288602.6",
}

// hardSymbols include a retired symbol, a readthrough, and one no gene claims.
var hardSymbols = []string{"BRAF", "TP53", "NSD2", "WHSC1", "RTEL1", "RTEL1-TNFRSF6B", "CDKN2A", "NOT_A_GENE"}

// apiCases expands the corpus into one request per endpoint per input.
func apiCases() []apiCase {
	var cases []apiCase
	add := func(name, method, path, body, sortKey string) {
		cases = append(cases, apiCase{Name: name, Method: method, Path: path, Body: body, SortKey: sortKey})
	}

	// Annotation, the endpoints that carry the most surface.
	for _, v := range hardVariants {
		add("annotation genomic GET "+v.name, "GET",
			"/annotation/genomic/"+v.genomic+"?fields=annotation_summary&isoformOverrideSource=mskcc", "", "")
		add("annotation genomic GET all fields "+v.name, "GET",
			"/annotation/genomic/"+v.genomic+
				"?fields=annotation_summary,clinvar,hotspots,ptms,signal&isoformOverrideSource=mskcc", "", "")
	}
	add("annotation genomic POST batch", "POST",
		"/annotation/genomic?fields=annotation_summary&isoformOverrideSource=mskcc",
		genomicLocationsJSON(), "")

	// The isoform override sources, which pick different transcripts.
	for _, source := range []string{"mskcc", "uniprot", "ensembl", "genome_nexus", ""} {
		label := source
		if label == "" {
			label = "unset"
		}
		add("canonical transcript CDKN2A isoform "+label, "GET",
			"/ensembl/canonical-transcript/hgnc/CDKN2A?isoformOverrideSource="+source, "", "")
	}

	for _, sym := range hardSymbols {
		add("canonical transcript "+sym, "GET", "/ensembl/canonical-transcript/hgnc/"+sym, "", "")
		add("canonical gene "+sym, "GET", "/ensembl/canonical-gene/hgnc/"+sym, "", "")
	}
	add("canonical transcript batch", "POST", "/ensembl/canonical-transcript/hgnc",
		jsonList(hardSymbols), "transcriptId")
	add("canonical gene batch", "POST", "/ensembl/canonical-gene/hgnc",
		jsonList(hardSymbols), "hugoSymbol")

	for _, tx := range hardTranscripts {
		add("transcript "+tx, "GET", "/ensembl/transcript/"+tx, "", "")
		add("ptms "+tx, "GET", "/ptm/experimental?ensemblTranscriptId="+tx, "", "")
		add("hotspots transcript "+tx, "GET", "/cancer_hotspots/transcript/"+tx, "", "residue")
		add("hotspots transcript v2 "+tx, "GET", "/cancer_hotspots/transcript/"+tx+"?version=v2", "", "residue")
	}
	add("transcript filter by ids", "POST", "/ensembl/transcript",
		`{"transcriptIds":`+jsonList(hardTranscripts)+`}`, "transcriptId")
	add("transcript filter by symbols", "POST", "/ensembl/transcript",
		`{"hugoSymbols":["BRAF"]}`, "transcriptId")
	add("transcript filter empty", "POST", "/ensembl/transcript", `{}`, "transcriptId")
	add("transcript list no params", "GET", "/ensembl/transcript", "", "transcriptId")

	for _, v := range hardVariants {
		add("hotspots genomic "+v.name, "GET", "/cancer_hotspots/genomic/"+v.genomic, "", "residue")
	}
	add("hotspots genomic batch", "POST", "/cancer_hotspots/genomic", genomicLocationsJSON(), "")

	add("pfam domain PF00069", "GET", "/pfam/domain/PF00069", "", "")
	add("pfam domain batch", "POST", "/pfam/domain", `["PF00069","PF07714"]`, "pfamAccession")
	add("version", "GET", "/version", "", "")

	return cases
}

// genomicLocationsJSON renders the hard variants as the body the genomic
// endpoints take.
func genomicLocationsJSON() string {
	var b strings.Builder
	b.WriteByte('[')
	for i, v := range hardVariants {
		parts := strings.Split(v.genomic, ",")
		if len(parts) != 5 {
			continue
		}
		if i > 0 {
			b.WriteByte(',')
		}
		fmt.Fprintf(&b, `{"chromosome":%q,"start":%s,"end":%s,"referenceAllele":%q,"variantAllele":%q}`,
			parts[0], parts[1], parts[2], parts[3], parts[4])
	}
	b.WriteByte(']')
	return b.String()
}

// jsonList renders a list of identifiers as a JSON array.
func jsonList(items []string) string {
	out, _ := json.Marshal(items)
	return string(out)
}
