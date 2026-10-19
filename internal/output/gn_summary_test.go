package output

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestFormatHGVSg(t *testing.T) {
	tests := []struct {
		name                 string
		chrom                string
		start, end           int64
		ref, alt             string
		want, wantWithRefSeq string
	}{
		{
			name:  "substitution carries both alleles either way",
			chrom: "2", start: 216809708, end: 216809708, ref: "C", alt: "T",
			want: "2:g.216809708C>T",
		},
		{
			name:  "single base deletion",
			chrom: "13", start: 32914438, end: 32914438, ref: "T", alt: "-",
			want: "13:g.32914438del",
		},
		{
			name:  "multi base deletion names the range",
			chrom: "1", start: 206811015, end: 206811016, ref: "AC", alt: "-",
			want: "1:g.206811015_206811016del",
		},
		{
			name:  "insertion falls between two bases and removes nothing",
			chrom: "17", start: 36002277, end: 36002278, ref: "-", alt: "A",
			want: "17:g.36002277_36002278insA",
		},
		{
			name:  "multi base delins splits into del and ins",
			chrom: "2", start: 216809708, end: 216809709, ref: "CA", alt: "T",
			want: "2:g.216809708_216809709delinsT",
		},
		{
			name:  "single base delins",
			chrom: "17", start: 7579363, end: 7579363, ref: "A", alt: "TTT",
			want: "17:g.7579363delinsTTT",
		},
		{
			name:  "a position before the first base cannot be named",
			chrom: "1", start: 0, end: 0, ref: "A", alt: "T",
			want: "",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			assert.Equal(t, tt.want, formatHGVSg(tt.chrom, tt.start, tt.end, tt.ref, tt.alt))
		})
	}
}

func TestGNVariantClassification(t *testing.T) {
	tests := []struct {
		name        string
		term        string
		variantType string
		inframe     bool
		want        string
	}{
		{"a frameshift deletion", "frameshift_variant", "DEL", false, "Frame_Shift_Del"},
		{"a frameshift insertion", "frameshift_variant", "INS", false, "Frame_Shift_Ins"},
		{
			name: "a frameshift that preserves the frame keeps the bare name",
			term: "frameshift_variant", variantType: "DEL", inframe: true, want: "Frame_Shift",
		},
		{
			name: "the alleles decide what a protein altering variant is",
			term: "protein_altering_variant", variantType: "DEL", inframe: true, want: "In_Frame_Del",
		},
		{
			name: "and what a coding sequence variant is",
			term: "coding_sequence_variant", variantType: "INS", inframe: true, want: "In_Frame_Ins",
		},
		{"a missense substitution", "missense_variant", "SNP", true, "Missense_Mutation"},
		{
			name: "a polypyrimidine tract variant on its own is an intron variant",
			term: "splice_polypyrimidine_tract_variant", variantType: "SNP", inframe: true, want: "Intron",
		},
		{
			name: "an unknown term without alleles has nothing to go on",
			term: "something_new", variantType: "", want: "Targeted_Region",
		},
		{
			name: "a splice term is classified even without alleles",
			term: "splice_acceptor_variant", variantType: "", want: "Splice_Site",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			assert.Equal(t, tt.want, gnVariantClassification(tt.term, tt.variantType, tt.inframe))
		})
	}
}

func TestHgvspShortFromHgvsc(t *testing.T) {
	tests := []struct {
		name          string
		tc            GNTranscriptConsequence
		proteinLength int
		variantType   string
		inframe       bool
		want          string
	}{
		{
			name: "a splice site is named for the residue it sits beside",
			tc: GNTranscriptConsequence{
				HGVSc:            "ENST00000269305.4:c.151-2A>G",
				ConsequenceTerms: []string{"splice_acceptor_variant"},
			},
			want: "p.X51_splice",
		},
		{
			name: "a splice site before the protein starts is named for its first residue",
			tc: GNTranscriptConsequence{
				HGVSc:            "ENST00000269305.4:c.-29+2T>C",
				ConsequenceTerms: []string{"splice_donor_variant"},
			},
			want: "p.X1_splice",
		},
		{
			name: "a splice site past the stop is named for the last residue",
			tc: GNTranscriptConsequence{
				HGVSc:            "ENST00000269305.4:c.*8+2T>C",
				ConsequenceTerms: []string{"splice_donor_variant"},
			},
			proteinLength: 393,
			want:          "p.X393_splice",
		},
		{
			name: "a range starting before the protein is named from its coding end",
			tc: GNTranscriptConsequence{
				HGVSc:            "ENST00000349496.5:c.-309_60+142del",
				ConsequenceTerms: []string{"splice_donor_variant"},
			},
			want: "p.X20_splice",
		},
		{
			name: "a splice region variant is not a splice site",
			tc: GNTranscriptConsequence{
				HGVSc:            "ENST00000342988.3:c.1140-25_1140-2del",
				ConsequenceTerms: []string{"splice_region_variant"},
			},
			variantType: "DEL",
			want:        "p.*380*",
		},
		{
			name: "a frameshift with no amino acids reaches a stop",
			tc: GNTranscriptConsequence{
				HGVSc:            "ENST00000000001.1:c.100del",
				ConsequenceTerms: []string{"frameshift_variant"},
			},
			variantType: "DEL",
			want:        "p.*34fs*",
		},
		{
			name: "the 3'UTR has no protein position to report",
			tc: GNTranscriptConsequence{
				HGVSc:            "ENST00000000001.1:c.*42A>G",
				ConsequenceTerms: []string{"3_prime_UTR_variant"},
			},
			want: "",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := hgvspShortFromHgvsc(&tt.tc, tt.proteinLength, tt.variantType, tt.inframe)
			assert.Equal(t, tt.want, got)
		})
	}
}

func TestShortenHgvsp(t *testing.T) {
	assert.Equal(t, "p.V600E", shortenHgvsp("p.Val600Glu"))
	assert.Equal(t, "p.Q1067*", shortenHgvsp("p.Gln1067Ter"))
	assert.Equal(t, "p.=", normalizeHgvsp("ENSP00000269305.4:(p.%3D)"))
	assert.Equal(t, "p.Val600Glu", normalizeHgvsp("ENSP00000288602.6:p.Val600Glu"))
}

func TestGNResolveProteinPosition(t *testing.T) {
	// A duplication and a splice marker carry a position the VEP-level
	// protein_start does not, so the name is read back for those.
	dup := GNTranscriptConsequence{ProteinStart: 10, ProteinEnd: 10}
	assert.Equal(t, &GNIntegerRange{Start: 44, End: 45}, gnResolveProteinPosition(&dup, "p.A44_P45dup"))
	assert.Equal(t, &GNIntegerRange{Start: 51, End: 51}, gnResolveProteinPosition(&dup, "p.X51_splice"))

	// Everything else uses the reported range.
	assert.Equal(t, &GNIntegerRange{Start: 10, End: 10}, gnResolveProteinPosition(&dup, "p.V600E"))

	none := GNTranscriptConsequence{}
	assert.Nil(t, gnResolveProteinPosition(&none, ""))
}

func TestResolveCanonicalNarrowsToTheRightGene(t *testing.T) {
	overrides := map[string]bool{"ENST00000360203": true, "ENST00000482936": true}
	cancerGenes := map[string]bool{"RTEL1": true}

	readthrough := &annotate.Annotation{
		TranscriptID: "ENST00000482936", GeneName: "RTEL1-TNFRSF6B", GeneID: "ENSG00000206190",
		Biotype: "nonsense_mediated_decay", Consequence: "missense_variant,NMD_transcript_variant",
	}
	gene := &annotate.Annotation{
		TranscriptID: "ENST00000360203", GeneName: "RTEL1", GeneID: "ENSG00000258366",
		Biotype: "protein_coding", Consequence: "missense_variant",
	}

	// The readthrough is listed first, so only the narrowing can put the row
	// on RTEL1.
	anns := []*annotate.Annotation{readthrough, gene}
	opt := GNMarshalOptions{OverrideTranscripts: overrides, CancerGenes: cancerGenes}
	got, marked := resolveCanonical(anns, opt, "missense_variant")
	require.NotNil(t, got)
	assert.Equal(t, "ENST00000360203", got.TranscriptID)
	assert.True(t, marked["ENST00000360203"])

	// With no override table at all the annotations keep their own flags.
	gene.IsCanonicalEnsembl = true
	got, _ = resolveCanonical(anns, GNMarshalOptions{}, "missense_variant")
	require.NotNil(t, got)
	assert.Equal(t, "ENST00000360203", got.TranscriptID)
}

func TestAlphaMissenseFromAnnotation(t *testing.T) {
	// A scored missense substitution carries both halves of the call.
	ann := &annotate.Annotation{}
	ann.SetExtraKey("alphamissense.score", "0.9871")
	ann.SetExtraKey("alphamissense.class", "likely_pathogenic")

	am := alphaMissenseFrom(ann)
	require.NotNil(t, am)
	assert.Equal(t, "likely_pathogenic", am.Pathogenicity)
	require.NotNil(t, am.Score)
	assert.InDelta(t, 0.9871, *am.Score, 1e-9)

	// A variant AlphaMissense does not score reports nothing at all, rather
	// than an empty object.
	assert.Nil(t, alphaMissenseFrom(&annotate.Annotation{}))

	// A class with no parseable score still reports the class.
	partial := &annotate.Annotation{}
	partial.SetExtraKey("alphamissense.class", "ambiguous")
	am = alphaMissenseFrom(partial)
	require.NotNil(t, am)
	assert.Equal(t, "ambiguous", am.Pathogenicity)
	assert.Nil(t, am.Score)
}
