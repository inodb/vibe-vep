package annotate

import "testing"

// Cases taken from a VEP111 MSK-IMPACT MAF: the Codons string VEP reports
// alongside its Consequence, for frameshift-length indels.
func TestAddStopGainedForFrameshift(t *testing.T) {
	tests := []struct {
		name        string
		codons      string
		consequence string
		want        string
	}{
		// The inserted bases spell a stop outright.
		{"insertion spelling TAA", "tac/taAAc", ConsequenceFrameshiftVariant,
			ConsequenceStopGained + "," + ConsequenceFrameshiftVariant},
		{"insertion spelling TGA", "tac/tGACTac", ConsequenceFrameshiftVariant,
			ConsequenceStopGained + "," + ConsequenceFrameshiftVariant},
		{"deletion joining into TAG", "taTAgt/tagt", ConsequenceFrameshiftVariant,
			ConsequenceStopGained + "," + ConsequenceFrameshiftVariant},
		// The stop is further along the inserted run, past the first codon,
		// and the protein still reads p.T20Sfs*6.
		{"stop in a later inserted codon", "gct/gcCTCTGGCCCTGGCTGCTGATGGCt",
			ConsequenceFrameshiftVariant,
			ConsequenceStopGained + "," + ConsequenceFrameshiftVariant},
		// p.S70*, but the affected codons do not spell a stop: the shifted
		// frame finds one downstream, which is not what VEP tests.
		{"downstream stop is not stop_gained", "caGTCA/ca", ConsequenceFrameshiftVariant,
			ConsequenceFrameshiftVariant},
		{"ordinary frameshift", "cCt/ct", ConsequenceFrameshiftVariant,
			ConsequenceFrameshiftVariant},
		// TGA appears out of frame (position 4), so it is not a stop codon.
		{"out-of-frame stop ignored", "gaa/gaTGAa", ConsequenceFrameshiftVariant,
			ConsequenceFrameshiftVariant},
		// Co-terms survive the upgrade and keep Ensembl's severity order.
		{"splice co-term kept", "tac/taAAc",
			ConsequenceFrameshiftVariant + "," + ConsequenceSpliceRegion,
			ConsequenceStopGained + "," + ConsequenceFrameshiftVariant + "," + ConsequenceSpliceRegion},
		// The reference already holds the stop, so it is lost, not gained.
		{"ref stop is not a gain", "tAAg/tg", ConsequenceFrameshiftStopLost,
			ConsequenceFrameshiftStopLost},
		// Nothing to test without codons, and non-frameshifts are left alone.
		{"no codons", "", ConsequenceFrameshiftVariant, ConsequenceFrameshiftVariant},
		{"not a frameshift", "tac/taAAc", ConsequenceMissenseVariant, ConsequenceMissenseVariant},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			r := &ConsequenceResult{Consequence: tt.consequence, CodonChange: tt.codons}
			addStopGainedFromCodons(r)
			if r.Consequence != tt.want {
				t.Errorf("Codons %q, %q -> %q, want %q", tt.codons, tt.consequence, r.Consequence, tt.want)
			}
		})
	}
}

// The same test applies to in-frame indels and delins, where VEP keeps the
// length term alongside: stop_gained,inframe_insertion for SMAD2 gaa/gaGTAa.
func TestAddStopGainedFromCodonsOnInFrameIndels(t *testing.T) {
	tests := []struct {
		name        string
		codons      string
		consequence string
		want        string
	}{
		{"SMAD2 gaa/gaGTAa", "gaa/gaGTAa", ConsequenceInframeInsertion,
			ConsequenceStopGained + "," + ConsequenceInframeInsertion},
		{"ERCC4 tTGAag/tag", "tTGAag/tag", ConsequenceInframeDeletion,
			ConsequenceStopGained + "," + ConsequenceInframeDeletion},
		{"GNAS cgCTCCAACGag/cgTTag", "cgCTCCAACGag/cgTTag", ConsequenceProteinAltering,
			ConsequenceStopGained + "," + ConsequenceProteinAltering},
		// The reference already holds the stop, so it is lost, not gained.
		{"ref already stops", "tTAAg/tg", ConsequenceInframeDeletion, ConsequenceInframeDeletion},
		// No stop in the altered codons.
		{"no stop", "gaa/gaCCCa", ConsequenceInframeInsertion, ConsequenceInframeInsertion},
		// Substitutions are left to the single-base path.
		{"substitution untouched", "tGG/tAA", ConsequenceMissenseVariant, ConsequenceMissenseVariant},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			r := &ConsequenceResult{Consequence: tt.consequence, CodonChange: tt.codons}
			addStopGainedFromCodons(r)
			if r.Consequence != tt.want {
				t.Errorf("Codons %q, %q -> %q, want %q", tt.codons, tt.consequence, r.Consequence, tt.want)
			}
		})
	}
}

// Whether VEP keeps the length term alongside stop_gained. Cases from a
// VEP111 MSK-IMPACT MAF, with the consequence VEP reports.
func TestStopGainedDropsTheLengthTerm(t *testing.T) {
	tests := []struct {
		name        string
		codons      string
		consequence string
		want        string
	}{
		// The reference codon is replaced by a stop: nothing survives there,
		// so there is no replacement left to describe.
		{"delins whose first codon becomes a stop", "CCGGga/Tga", ConsequenceProteinAltering,
			ConsequenceStopGained},
		{"insertion whose first codon becomes a stop", "tct/tGATGGTCCCAGAGTTct",
			ConsequenceInframeInsertion, ConsequenceStopGained},
		// The first altered codon still spells a residue.
		{"delins whose first codon still codes", "cgCTCCAACGag/cgTTag", ConsequenceProteinAltering,
			ConsequenceStopGained + "," + ConsequenceProteinAltering},
		{"insertion after a surviving residue", "gaa/gaGTAa", ConsequenceInframeInsertion,
			ConsequenceStopGained + "," + ConsequenceInframeInsertion},
		// Whole codons inserted at a boundary replace no reference codon.
		{"boundary insertion", "-/TAACAA", ConsequenceInframeInsertion,
			ConsequenceStopGained + "," + ConsequenceInframeInsertion},
		// A pure deletion keeps the term even when the junction spells a stop.
		{"pure deletion", "tTGAag/tag", ConsequenceInframeDeletion,
			ConsequenceStopGained + "," + ConsequenceInframeDeletion},
		// Frameshifts always keep it.
		{"frameshift", "tct/tAGTct", ConsequenceFrameshiftVariant,
			ConsequenceStopGained + "," + ConsequenceFrameshiftVariant},
		// Splice co-terms are not length terms and survive either way.
		{"splice term survives the drop", "CCGGga/Tga",
			ConsequenceProteinAltering + "," + ConsequenceSpliceRegion,
			ConsequenceStopGained + "," + ConsequenceSpliceRegion},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			r := &ConsequenceResult{Consequence: tt.consequence, CodonChange: tt.codons}
			addStopGainedFromCodons(r)
			if r.Consequence != tt.want {
				t.Errorf("Codons %q, %q -> %q, want %q", tt.codons, tt.consequence, r.Consequence, tt.want)
			}
		})
	}
}
