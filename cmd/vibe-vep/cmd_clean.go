package main

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/spf13/cobra"
)

func newCleanCmd() *cobra.Command {
	var (
		assembly string
		all      bool
	)

	cmd := &cobra.Command{
		Use:   "clean",
		Short: "Remove raw downloaded data files",
		Long: `Remove raw/ subdirectories containing downloaded source data.
Prepared indexes (transcripts.gob, *.sqlite) are kept intact.
Run this after 'vibe-vep prepare' to reclaim disk space.

The checksums.sha256 file in each assembly directory is preserved
so you can verify which raw files were used to build the indexes.`,
		Example: `  vibe-vep clean --assembly GRCh38
  vibe-vep clean --all`,
		Args: cobra.NoArgs,
		RunE: func(cmd *cobra.Command, args []string) error {
			if all {
				return cleanAll()
			}
			normalized, err := normalizeAssembly(assembly)
			if err != nil {
				return err
			}
			return cleanAssembly(normalized)
		},
	}

	cmd.Flags().StringVar(&assembly, "assembly", "GRCh38", "Genome assembly to clean")
	cmd.Flags().BoolVar(&all, "all", false, "Clean raw data for all assemblies")

	return cmd
}

func cleanAssembly(assembly string) error {
	dir := DefaultGENCODEPath(assembly)
	if dir == "" {
		return fmt.Errorf("cannot determine data directory")
	}

	rawDir := filepath.Join(dir, "raw")
	info, err := os.Stat(rawDir)
	if err != nil {
		if os.IsNotExist(err) {
			fmt.Printf("No raw data to clean for %s\n", assembly)
			return nil
		}
		return err
	}
	if !info.IsDir() {
		return nil
	}

	// Show what we're deleting.
	var totalSize int64
	filepath.Walk(rawDir, func(_ string, fi os.FileInfo, _ error) error {
		if fi != nil && !fi.IsDir() {
			totalSize += fi.Size()
		}
		return nil
	})

	fmt.Printf("Removing %s (%s)...\n", rawDir, formatSize(totalSize))
	if err := os.RemoveAll(rawDir); err != nil {
		return fmt.Errorf("removing %s: %w", rawDir, err)
	}
	fmt.Printf("Done. Checksums preserved in %s\n", filepath.Join(dir, "checksums.sha256"))
	return nil
}

func cleanAll() error {
	home, err := os.UserHomeDir()
	if err != nil {
		return fmt.Errorf("cannot determine home directory: %w", err)
	}
	baseDir := filepath.Join(home, ".vibe-vep")

	entries, err := os.ReadDir(baseDir)
	if err != nil {
		return fmt.Errorf("reading %s: %w", baseDir, err)
	}

	for _, e := range entries {
		if !e.IsDir() {
			continue
		}
		asm := strings.ToUpper(e.Name())
		if asm == "GRCH37" || asm == "GRCH38" {
			if err := cleanAssembly(asm); err != nil {
				return err
			}
		}
	}
	return nil
}
