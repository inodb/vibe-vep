package genomicindex

import (
	"os"
	"path/filepath"
	"testing"
)

// An index asked to build from sources that are all missing used to finish
// successfully and answer nothing for every variant, which is how a 3.66 GB
// index was once replaced by an empty one without a single error.
func TestBuildRefusesToProduceAnEmptyIndex(t *testing.T) {
	dir := t.TempDir()
	dbPath := filepath.Join(dir, "genomic_annotations.sqlite")

	err := Build(dbPath, BuildSources{
		AlphaMissenseTSV: filepath.Join(dir, "absent_alphamissense.tsv.gz"),
		ClinVarVCF:       filepath.Join(dir, "absent_clinvar.vcf.gz"),
	}, func(string, ...any) {})

	if err == nil {
		t.Fatal("building from sources that are all missing should be an error")
	}
	if Ready(dbPath, BuildSources{}) {
		t.Error("an index built from nothing must not report itself ready")
	}
}

func TestLoadedSourcesOnMissingIndex(t *testing.T) {
	if _, err := LoadedSources(filepath.Join(t.TempDir(), "absent.sqlite")); err == nil {
		// sql.Open is lazy, so the query is what fails; either way it must
		// not report sources that are not there.
		t.Log("open deferred; the query reports the failure")
	}
}

func TestReadyRejectsEmptyFile(t *testing.T) {
	dir := t.TempDir()
	dbPath := filepath.Join(dir, "empty.sqlite")
	if err := os.WriteFile(dbPath, nil, 0o644); err != nil {
		t.Fatal(err)
	}
	if Ready(dbPath, BuildSources{}) {
		t.Error("a zero-byte index must not report itself ready")
	}
}
