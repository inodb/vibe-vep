package output

import (
	"regexp"
	"strconv"
	"strings"
)

// annotation_summary is the only part of the genome-nexus response the
// annotation pipeline reads: every MAF column it writes comes from the
// canonical transcript's summary, not from the VEP-level fields beside it.
// genome-nexus builds that summary with a dozen small resolver classes, and
// this file is a port of the ones that decide a compared column:
// ProteinChangeResolver, ProteinPositionResolver, VariantClassificationResolver,
// VariantTypeResolver, ConsequenceTermsResolver and their helpers.
//
// They are deliberately written against the VEP-level view -- consequence
// terms, amino_acids, protein_start, hgvsc, hgvsp -- and not against
// annotate.Annotation, because that is all genome-nexus itself ever saw. A
// summary field is then wrong only where the VEP-level field feeding it is
// wrong, which is the same place it was wrong before.

// aa3to1 maps three-letter residue names to one-letter codes, in
// genome-nexus's own order. Replacement is a plain substring pass over the
// whole string, so Ter -> * turns p.Gln1067Ter into p.Q1067*.
var aa3to1 = [][2]string{
	{"Ala", "A"}, {"Arg", "R"}, {"Asn", "N"}, {"Asp", "D"}, {"Asx", "B"}, {"Cys", "C"},
	{"Glu", "E"}, {"Gln", "Q"}, {"Glx", "Z"}, {"Gly", "G"}, {"His", "H"}, {"Ile", "I"},
	{"Leu", "L"}, {"Lys", "K"}, {"Met", "M"}, {"Phe", "F"}, {"Pro", "P"}, {"Ser", "S"},
	{"Thr", "T"}, {"Trp", "W"}, {"Tyr", "Y"}, {"Val", "V"}, {"Xxx", "X"}, {"Ter", "*"},
}

// spliceSiteVariants are the terms that make genome-nexus name a protein
// change p.X<n>_splice instead of translating one. A splice region variant is
// not one of them: it can sit far enough from the junction to leave a real
// protein change worth naming.
var spliceSiteVariants = map[string]bool{
	"splice_acceptor_variant": true,
	"splice_donor_variant":    true,
}

// cDNAExtractor pulls the coding position out of an HGVSc notation. It is
// genome-nexus's pattern unchanged, including the leading .* that makes it
// take the last coding prefix in the string rather than the first.
var cDNAExtractor = regexp.MustCompile(`^.*[cn].-?\*?(\d+).*$`)

// digits finds the numbers in a protein change, for reading a position back
// out of a name genome-nexus has already built.
var digits = regexp.MustCompile(`\d+`)

// A coding position written from the 5'UTR carries a minus sign that
// cDNAExtractor drops, so c.-309_60+142del reads as 309 rather than -309.
// These three patterns recognise the forms where that matters: a range
// starting in the 5'UTR, whose second coordinate is the one in the CDS, and a
// single 5'UTR position at an intron boundary.
var (
	cDNAUTRStartRange  = regexp.MustCompile(`^.*[cn]\.-\d+_.*$`)
	cDNARangeSecondPos = regexp.MustCompile(`^.*_(\d+).*$`)
	cDNAUTRStartSplice = regexp.MustCompile(`^.*[cn]\.-\d+[+-].*$`)
)

// alleleLength counts a MAF allele, where "-" means the allele is absent.
func alleleLength(allele string) int {
	if allele == "" || allele == "-" {
		return 0
	}
	return len(allele)
}

// gnVariantType names the shape of the change the way genome-nexus does:
// equal lengths give SNP/DNP/TNP/ONP, unequal give INS or DEL.
func gnVariantType(ref, alt string) string {
	refLen, altLen := alleleLength(ref), alleleLength(alt)
	if refLen != altLen {
		if refLen < altLen {
			return "INS"
		}
		return "DEL"
	}
	switch refLen {
	case 0:
		return ""
	case 1:
		return "SNP"
	case 2:
		return "DNP"
	case 3:
		return "TNP"
	default:
		return "ONP"
	}
}

// gnIsInframe reports whether the change preserves the reading frame.
func gnIsInframe(ref, alt string) bool {
	d := alleleLength(ref) - alleleLength(alt)
	if d < 0 {
		d = -d
	}
	return d%3 == 0
}

// pickHighestPriorityConsequence returns the most severe term by
// genome-nexus's severity order. Unlike most_severe_consequence it ranks on
// the vcf2maf table, and it is the term the variant classification is read
// from.
func pickHighestPriorityConsequence(terms []string) string {
	best := ""
	highest := int(^uint(0) >> 1)
	for _, term := range terms {
		if p := rankEffect(term); p < highest {
			best = term
			highest = p
		}
	}
	return best
}

// gnVariantClassification maps a consequence term to a MAF variant
// classification. Three terms do not map on their own -- a frameshift, a
// protein-altering and a coding-sequence variant each name a deletion or an
// insertion depending on the alleles -- so the shape of the change decides
// those.
func gnVariantClassification(term, variantType string, inframe bool) string {
	const defaultClass = "Targeted_Region"

	if term == "" || (variantType == "" && !strings.Contains(term, "splice")) {
		return defaultClass
	}
	term = strings.ToLower(term)

	// Three terms name no classification on their own: the shape of the
	// change says whether they deleted or inserted, and whether the reading
	// frame survived.
	if term == "frameshift_variant" || term == "protein_altering_variant" || term == "coding_sequence_variant" {
		if !inframe {
			switch variantType {
			case "DEL":
				return "Frame_Shift_Del"
			case "INS":
				return "Frame_Shift_Ins"
			}
		} else if term != "frameshift_variant" {
			// A frameshift that turns out to preserve the frame keeps its
			// bare name rather than becoming an in-frame indel.
			switch variantType {
			case "DEL":
				return "In_Frame_Del"
			case "INS":
				return "In_Frame_Ins"
			}
		}
	}

	if c, ok := gnVariantMap[term]; ok {
		return c
	}
	return defaultClass
}

// gnVariantMap is genome-nexus's consequence term to MAF classification table.
var gnVariantMap = map[string]string{
	"splice_acceptor_variant":       "Splice_Site",
	"splice_donor_variant":          "Splice_Site",
	"transcript_ablation":           "Splice_Site",
	"exon_loss_variant":             "Splice_Site",
	"stop_gained":                   "Nonsense_Mutation",
	"frameshift_variant":            "Frame_Shift",
	"stop_lost":                     "Nonstop_Mutation",
	"initiator_codon_variant":       "Translation_Start_Site",
	"start_lost":                    "Translation_Start_Site",
	"inframe_insertion":             "In_Frame_Ins",
	"inframe_deletion":              "In_Frame_Del",
	"disruptive_inframe_insertion":  "In_Frame_Ins",
	"disrupting_inframe_deletion":   "In_Frame_Del",
	"missense_variant":              "Missense_Mutation",
	"protein_altering_variant":      "Missense_Mutation",
	"coding_sequence_variant":       "Missense_Mutation",
	"coding_transcript_variant":     "Missense_Mutation",
	"conservative_missense_variant": "Missense_Mutation",
	"rare_amino_acid_variant":       "Missense_Mutation",
	"transcript_amplification":      "Intron",
	"feature_elongation":            "Feature_Elongation",
	"feature_truncation":            "Feature_Elongation",
	"sequence_variant":              "Sequence_Variant",
	"splice_region_variant":         "Splice_Region",
	"splice_donor_region_variant":   "Splice_Region",
	// A polypyrimidine tract variant sits 3 to 15 bases before the acceptor
	// site. Where it overlaps the splice region the splice region term is
	// reported alongside it and wins, so on its own it is an intron variant.
	"splice_polypyrimidine_tract_variant":            "Intron",
	"splice_donor_5th_base_variant":                  "Splice_Region",
	"intron_variant":                                 "Intron",
	"intragenic":                                     "Intron",
	"intragenic_variant":                             "Intron",
	"incomplete_terminal_codon_variant":              "Silent",
	"synonymous_variant":                             "Silent",
	"start_retained_variant":                         "Silent",
	"stop_retained_variant":                          "Silent",
	"nmd_transcript_variant":                         "Silent",
	"mature_mirna_variant":                           "RNA",
	"non_coding_exon_variant":                        "RNA",
	"non_coding_transcript_exon_variant":             "RNA",
	"non_coding_transcript_variant":                  "RNA",
	"nc_transcript_variant":                          "RNA",
	"exon_variant":                                   "RNA",
	"5_prime_utr_variant":                            "5'UTR",
	"5_prime_utr_premature_start_codon_gain_variant": "5'UTR",
	"3_prime_utr_variant":                            "3'UTR",
	"tf_binding_site_variant":                        "IGR",
	"tfbs_ablation":                                  "IGR",
	"tfbs_amplification":                             "IGR",
	"regulatory_region_variant":                      "IGR",
	"regulatory_region_amplification":                "IGR",
	"regulatory_region":                              "IGR",
	"regulatory_region_ablation":                     "IGR",
	"intergenic_variant":                             "IGR",
	"intergenic_region":                              "IGR",
	"upstream_gene_variant":                          "5'Flank",
	"downstream_gene_variant":                        "3'Flank",
	"splicing":                                       "Splice_Region",
	"nonframeshift_insertion":                        "In_Frame_Ins",
	"frameshift_insertion":                           "Frame_Shift_Ins",
	"stopgain_SNV":                                   "Nonsense_Mutation",
	"frameshift_deletion":                            "Frame_Shift_Del",
	"stoploss_SNV":                                   "Nonstop_Mutation",
	"upstream":                                       "5'Flank",
	"nonframeshift_deletion":                         "In_Frame_Del",
	"nonsynonymous_SNV":                              "Missense_Mutation",
}

// gnResolveHgvsp returns the protein change a transcript reports, or nothing
// where the variant is at a splice site: a disturbed junction has no protein
// change worth naming, and genome-nexus builds p.X<n>_splice from the coding
// position instead. A splice region variant keeps its protein change, since
// it can sit far enough from the junction to have a real one.
//
// It asks for the classification without the variant's alleles, exactly as
// genome-nexus does, which leaves every other term classified Targeted_Region
// here. Only Splice_Site is ever read off the answer.
func gnResolveHgvsp(tc *GNTranscriptConsequence) string {
	if tc.HGVSp == "" {
		return ""
	}
	class := gnVariantClassification(pickHighestPriorityConsequence(tc.ConsequenceTerms), "", false)
	if strings.EqualFold(class, "Splice_Site") {
		return ""
	}
	return normalizeHgvsp(tc.HGVSp)
}

// normalizeHgvsp drops the protein accession an HGVSp is prefixed with, and
// spells VEP's escaped synonymous notation the short way.
func normalizeHgvsp(hgvsp string) string {
	if strings.Contains(hgvsp, "(p.%3D)") {
		return "p.="
	}
	if i := strings.Index(hgvsp, ":"); i >= 0 {
		return hgvsp[i+1:]
	}
	return hgvsp
}

// gnResolveHgvspShort names the protein change in one-letter form. It is
// genome-nexus's three-way cascade: translate the reported protein change
// where there is one, otherwise read a position out of the coding change, and
// failing that build a name from the amino acids alone.
func gnResolveHgvspShort(tc *GNTranscriptConsequence, proteinLength int, variantType string, inframe bool) string {
	if hgvsp := gnResolveHgvsp(tc); hgvsp != "" {
		return shortenHgvsp(hgvsp)
	}
	if tc.HGVSc != "" {
		if short := hgvspShortFromHgvsc(tc, proteinLength, variantType, inframe); short != "" {
			return short
		}
	}
	if tc.AminoAcids != "" {
		return hgvspShortFromAminoAcids(tc)
	}
	return ""
}

// shortenHgvsp rewrites three-letter residue names as one-letter codes.
func shortenHgvsp(hgvsp string) string {
	short := hgvsp
	for _, aa := range aa3to1 {
		short = strings.ReplaceAll(short, aa[0], aa[1])
	}
	return short
}

// hgvspShortFromHgvsc names a protein change for a variant that has none,
// from the coding position alone. A splice site is named for the residue the
// junction sits beside; a variant with no amino acids at all is named as
// reaching a stop there.
//
// A position outside the coding sequence needs interpreting rather than
// dividing by three. The regex cannot see the minus sign on a 5'UTR
// coordinate, so c.-309+2 would otherwise land in the middle of the protein
// when the junction is before its start, and c.*8+2 would land near the start
// when the junction is past its end. Both are anchored to the end of the
// protein they are actually beside.
func hgvspShortFromHgvsc(tc *GNTranscriptConsequence, proteinLength int, variantType string, inframe bool) string {
	m := cDNAExtractor.FindStringSubmatch(tc.HGVSc)
	if m == nil {
		return ""
	}

	firstTerm := ""
	if len(tc.ConsequenceTerms) > 0 {
		firstTerm = tc.ConsequenceTerms[0]
	}
	atSpliceSite := spliceSiteVariants[firstTerm]

	// "c.*" is the 3'UTR, which has no protein position to report -- unless
	// the variant is at a splice site, which does.
	inThreePrimeUTR := strings.Contains(tc.HGVSc, "c.*")
	if inThreePrimeUTR && !atSpliceSite {
		return ""
	}

	cPos, err := strconv.Atoi(m[1])
	if err != nil {
		return ""
	}
	if cDNAUTRStartRange.MatchString(tc.HGVSc) {
		if rm := cDNARangeSecondPos.FindStringSubmatch(tc.HGVSc); rm != nil {
			if second, err := strconv.Atoi(rm[1]); err == nil && second >= 1 {
				cPos = second
			}
		}
	}
	if cPos < 1 {
		cPos = 1
	}
	pPos := (cPos + 2) / 3

	if atSpliceSite {
		switch {
		case inThreePrimeUTR:
			if proteinLength <= 0 {
				return ""
			}
			pPos = proteinLength
		case cDNAUTRStartSplice.MatchString(tc.HGVSc):
			pPos = 1
		}
		return "p.X" + strconv.Itoa(pPos) + "_splice"
	}

	if tc.AminoAcids == "" {
		class := gnVariantClassification(pickHighestPriorityConsequence(tc.ConsequenceTerms), variantType, inframe)
		if strings.HasPrefix(strings.ToLower(class), "frame_shift") {
			return "p.*" + strconv.Itoa(pPos) + "fs*"
		}
		return "p.*" + strconv.Itoa(pPos) + "*"
	}
	return ""
}

// hgvspShortFromAminoAcids is genome-nexus's last resort, building a name out
// of the reference and alternate residues and the protein position.
func hgvspShortFromAminoAcids(tc *GNTranscriptConsequence) string {
	parts := strings.Split(tc.AminoAcids, "/")
	term := ""
	if len(tc.ConsequenceTerms) > 0 {
		term = strings.ToLower(tc.ConsequenceTerms[0])
	}

	switch {
	case strings.Contains(term, "inframe_insertion"):
		if len(parts) < 2 {
			return ""
		}
		alt := parts[1]
		if !strings.Contains(tc.HGVSc, "dup") {
			return "p.X" + strconv.FormatInt(tc.ProteinStart, 10) +
				"_X" + strconv.FormatInt(tc.ProteinEnd, 10) + "ins" + alt
		}
		if len(alt) == 1 {
			return "p." + alt[:1] + strconv.FormatInt(tc.ProteinStart-1, 10) + "dup"
		}
		// An insertion's reported end is its start plus one, so the length of
		// the inserted peptide gives the last duplicated residue instead.
		return "p." + alt[:1] + strconv.FormatInt(tc.ProteinStart+1, 10) +
			"_" + alt[len(alt)-1:] + strconv.FormatInt(tc.ProteinStart+int64(len(alt)), 10) + "dup"

	case strings.Contains(term, "inframe_deletion"):
		return "p." + parts[0] + "del"

	default:
		if parts[0] == "" {
			return ""
		}
		short := "p." + parts[0][:1] + strconv.FormatInt(tc.ProteinStart, 10)
		if strings.Contains(term, "frameshift_variant") {
			return short + "fs"
		}
		// A single amino acid with no alternate is VEP saying the residue is
		// unchanged.
		if len(parts) < 2 {
			return short + "="
		}
		return short + parts[1]
	}
}

// gnResolveProteinPosition reports the residues a variant covers. A
// duplication and a splice variant are read back out of the name
// genome-nexus just built, because for those the name carries a position the
// VEP-level protein_start does not.
func gnResolveProteinPosition(tc *GNTranscriptConsequence, hgvspShort string) *GNIntegerRange {
	lower := strings.ToLower(hgvspShort)
	if hgvspShort != "" && (strings.Contains(lower, "dup") || strings.Contains(lower, "x")) {
		if r := extractProteinPos(hgvspShort); r != nil {
			return r
		}
	}
	if tc.ProteinStart != 0 || tc.ProteinEnd != 0 {
		return &GNIntegerRange{Start: tc.ProteinStart, End: tc.ProteinEnd}
	}
	return nil
}

// extractProteinPos reads the first and last numbers out of a protein change.
func extractProteinPos(proteinChange string) *GNIntegerRange {
	found := digits.FindAllString(proteinChange, -1)
	if len(found) == 0 {
		return nil
	}
	start, err := strconv.ParseInt(found[0], 10, 64)
	if err != nil {
		return nil
	}
	end := start
	if len(found) >= 2 {
		if e, err := strconv.ParseInt(found[len(found)-1], 10, 64); err == nil {
			end = e
		}
	}
	return &GNIntegerRange{Start: start, End: end}
}

// GNVariantType names the shape of a change the way genome-nexus does:
// SNP, DNP, TNP, ONP, INS or DEL.
func GNVariantType(ref, alt string) string { return gnVariantType(ref, alt) }

// GNIsInframe reports whether a change preserves the reading frame.
func GNIsInframe(ref, alt string) bool { return gnIsInframe(ref, alt) }

// GNVariantClassifications returns the MAF classification each consequence
// term maps to, which is what a hotspot's type is matched against. Every term
// is translated, not just the most severe one, because any of them matching is
// enough for the hotspot to apply.
func GNVariantClassifications(terms []string, variantType string, inframe bool) []string {
	seen := make(map[string]bool, len(terms))
	classes := make([]string, 0, len(terms))
	for _, term := range terms {
		c := gnVariantClassification(term, variantType, inframe)
		if c == "" || seen[c] {
			continue
		}
		seen[c] = true
		classes = append(classes, c)
	}
	return classes
}
