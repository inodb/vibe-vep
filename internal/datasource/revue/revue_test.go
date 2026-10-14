package revue

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

const kitCatalogue = `[
  {
    "hugoGeneSymbol": "KIT",
    "revisedProteinEffects": [
      {
        "transcriptId": "ENST00000288135",
        "genomicLocation": "4,55593576,55593606,CCACAGAAACCCATGTATGAAGTACAGTGGA,-",
        "revisedProteinEffect": "p.K550_K558del",
        "revisedVariantClassification": "Splice_Exon_Shortening_In_Frame",
        "revisedStandardVariantClassification": "In_Frame_Del",
        "confirmed": true
      },
      {
        "transcriptId": "ENST00000288135",
        "genomicLocation": "4,55593577,55593577,A,-",
        "revisedProteinEffect": "p.K550Nfs*3",
        "revisedVariantClassification": "Splice_Exon_Shortening_Frame_Shift",
        "confirmed": false
      }
    ]
  }
]`

func TestParseAndLookup(t *testing.T) {
	store, err := Parse([]byte(kitCatalogue))
	require.NoError(t, err)
	assert.Equal(t, 2, store.Count())

	entry, ok := store.Lookup("ENST00000288135", "4,55593576,55593606,CCACAGAAACCCATGTATGAAGTACAGTGGA,-")
	require.True(t, ok)
	assert.Equal(t, "p.K550_K558del", entry.ProteinEffect)
	assert.Equal(t, "In_Frame_Del", entry.VariantClassification)
	assert.True(t, entry.Confirmed)

	// An entry with no standard classification falls back to translating
	// reVUE's own name for the mechanism.
	entry, ok = store.Lookup("ENST00000288135", "4,55593577,55593577,A,-")
	require.True(t, ok)
	assert.Equal(t, "Frame_Shift_Del", entry.VariantClassification)
	assert.False(t, entry.Confirmed)

	// Another transcript at the same location is a different variant.
	_, ok = store.Lookup("ENST00000412167", "4,55593576,55593606,CCACAGAAACCCATGTATGAAGTACAGTGGA,-")
	assert.False(t, ok)

	// A nil store answers for every variant without being checked for first.
	var missing *Store
	_, ok = missing.Lookup("ENST00000288135", "4,55593576,55593606,CCACAGAAACCCATGTATGAAGTACAGTGGA,-")
	assert.False(t, ok)
	assert.Equal(t, 0, missing.Count())
}

func TestStandardVariantClassification(t *testing.T) {
	// reVUE renamed Out_Of_Frame to Frame_Shift; both spellings name the same
	// mechanism and have to reach the same classification.
	assert.Equal(t, "Frame_Shift_Del", StandardVariantClassification("Splice_Exon_Skip_Out_Of_Frame"))
	assert.Equal(t, "Frame_Shift_Del", StandardVariantClassification("Splice_Exon_Skip_Frame_Shift"))
	assert.Equal(t, "Frame_Shift_Ins", StandardVariantClassification("Splice_Exon_Extension_Frame_Shift"))
	assert.Equal(t, "In_Frame_Del", StandardVariantClassification("Splice_Exon_Skip_In_Frame"))
	assert.Equal(t, "In_Frame_Ins", StandardVariantClassification("Splice_Intron_Retention_In_Frame"))

	// A mechanism with no MAF spelling leaves the prediction's own
	// classification in place.
	assert.Equal(t, "", StandardVariantClassification("Something_New"))
}
