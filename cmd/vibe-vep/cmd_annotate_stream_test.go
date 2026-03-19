package main

import (
	"testing"
)

func TestStreamSubcommandHelp(t *testing.T) {
	out, _, err := executeCommand("annotate", "stream", "--help")
	if err != nil {
		t.Fatalf("annotate stream --help failed: %v", err)
	}
	for _, want := range []string{
		"--input-format",
		"--output-format",
		"--assembly",
		"genome-nexus-genomic-location-jsonl",
		"ensembl-vep-jsonl",
		"vibe-vep-jsonl",
	} {
		if !containsString(out, want) {
			t.Errorf("annotate stream --help should mention %q", want)
		}
	}
}

func TestStreamInvalidInputFormat(t *testing.T) {
	_, _, err := executeCommand("annotate", "stream", "--input-format", "bad")
	if err == nil {
		t.Fatal("expected error for bad input format")
	}
}

func TestStreamInvalidOutputFormat(t *testing.T) {
	_, _, err := executeCommand("annotate", "stream", "--output-format", "bad")
	if err == nil {
		t.Fatal("expected error for bad output format")
	}
}
