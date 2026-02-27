package annotate

import "testing"

func TestGetImpact(t *testing.T) {
	tests := []struct {
		consequence string
		want        string
	}{
		{"missense_variant", ImpactModerate},
		{"stop_gained", ImpactHigh},
		{"synonymous_variant", ImpactLow},
		{"intron_variant", ImpactModifier},
		{"frameshift_variant,splice_region_variant", ImpactHigh},
		{"splice_region_variant,intron_variant", ImpactLow},
		{"missense_variant,splice_region_variant", ImpactModerate},
	}
	for _, tt := range tests {
		t.Run(tt.consequence, func(t *testing.T) {
			if got := GetImpact(tt.consequence); got != tt.want {
				t.Errorf("GetImpact(%q) = %q, want %q", tt.consequence, got, tt.want)
			}
		})
	}
}

// TestAllocRegression_GetImpact verifies zero allocations for GetImpact.
func TestAllocRegression_GetImpact(t *testing.T) {
	// Single consequence (fast path)
	allocs := testing.AllocsPerRun(100, func() {
		GetImpact("missense_variant")
	})
	if allocs > 0 {
		t.Errorf("GetImpact(single) allocs: %.0f, want 0", allocs)
	}

	// Compound consequence (comma-separated)
	allocs = testing.AllocsPerRun(100, func() {
		GetImpact("frameshift_variant,splice_region_variant")
	})
	if allocs > 0 {
		t.Errorf("GetImpact(compound) allocs: %.0f, want 0", allocs)
	}
}
