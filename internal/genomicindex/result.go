// Package genomicindex provides a unified SQLite-backed index for genomic
// annotation sources (AlphaMissense, ClinVar, SIGNAL). One database per assembly
// with a WITHOUT ROWID clustered primary key on (chrom, pos, ref, alt) gives
// ~1-5μs point lookups with near-zero Go heap via mmap.
//
// # Coordinate normalization
//
// All variants are stored with normalized coordinates:
//   - Chromosome: without "chr" prefix (e.g. "12", not "chr12")
//   - Position and alleles: canonical MAF-style (no anchor base for indels)
//
// VCF-style indels (with anchor base) are normalized via [NormalizeAlleles]:
// left/right-trim shared bases, adjust position. This ensures the same
// physical variant produces the same key whether it came from VCF or MAF input.
//
// AlphaMissense is SNV-only so format is unambiguous. ClinVar VCF indels are
// normalized during build. SIGNAL uses MAF format natively. gnomAD VCF indels
// are normalized during build. Lookups also normalize, so both VCF and MAF
// input match correctly.
package genomicindex

// Result holds the combined annotation data for a single genomic position.
type Result struct {
	AMScore      float32
	AMClass      string
	CVClnSig     string
	CVClnRevStat string
	CVClnDN      string
	SigMutStatus string
	SigCount     string
	SigFreq      string
	GnomadAF     string
	GnomadAC     string
	GnomadAN     string
	GnomadNhomalt string
	GnomadVersion string
}

// BuildSources holds paths to the source data files for building the index.
type BuildSources struct {
	AlphaMissenseTSV string // gzipped TSV (e.g. AlphaMissense_hg38.tsv.gz)
	ClinVarVCF       string // gzipped VCF (e.g. clinvar.vcf.gz)
	SignalTSV        string // plain TSV (e.g. signaldb_all_variants_frequencies.txt)
	GnomadVCF        string // gzipped VCF (e.g. gnomad.genomes.v4.1.sites.vcf.bgz)
	GnomadVersion    string // version string (e.g. "4.1" or "2.1.1")
}
