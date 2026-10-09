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
			addStopGainedForFrameshift(r)
			if r.Consequence != tt.want {
				t.Errorf("Codons %q, %q -> %q, want %q", tt.codons, tt.consequence, r.Consequence, tt.want)
			}
		})
	}
}
