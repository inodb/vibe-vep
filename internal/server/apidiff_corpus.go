//go:build apidiff

package server

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
