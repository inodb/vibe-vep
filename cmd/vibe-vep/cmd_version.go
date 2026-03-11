package main

import (
	"fmt"
	"os"
	"strings"
	"text/tabwriter"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

func newVersionCmd(verbose *bool) *cobra.Command {
	var mafColumns bool

	cmd := &cobra.Command{
		Use:   "version",
		Short: "Show version and data source information",
		Long:  "Show vibe-vep version, loaded data sources, and optionally the MAF column mapping.",
		PreRunE: func(cmd *cobra.Command, args []string) error {
			return viper.BindPFlags(cmd.Flags())
		},
		RunE: func(cmd *cobra.Command, args []string) error {
			fmt.Printf("%-15s%s\n", "vibe-vep", version)

			// Show GENCODE info
			assembly, _ := normalizeAssembly(viper.GetString("assembly"))
			if assembly == "" {
				assembly = "GRCh38"
			}
			if _, _, _, found := FindGENCODEFiles(assembly); found {
				fmt.Printf("%-15s%s (%s)\n", "GENCODE", GencodeVersionForAssembly(assembly), assembly)
			}

			// Show configured annotation sources (lightweight, no data loading).
			cacheDir := DefaultGENCODEPath(assembly)
			type sourceInfo struct {
				name, match, assembly, version, status string
				columns                                []annotate.ColumnDef
			}
			var infos []sourceInfo

			// fileModDate returns the file modification date as YYYY-MM-DD, or "" on error.
			fileModDate := func(path string) string {
				fi, err := os.Stat(path)
				if err != nil {
					return ""
				}
				return fi.ModTime().Format("2006-01-02")
			}

			// OncoKB
			if cglPath := viper.GetString("oncokb.cancer-gene-list"); cglPath != "" {
				status := "configured"
				ver := ""
				if fi, err := os.Stat(cglPath); err != nil {
					status = "file not found"
				} else {
					ver = fi.ModTime().Format("2006-01-02")
				}
				infos = append(infos, sourceInfo{"oncokb", string(annotate.MatchGene), "any", ver, status,
					[]annotate.ColumnDef{{Name: "gene_type", Description: "Gene classification (ONCOGENE/TSG)"}}})
			}

			// Genomic index (AlphaMissense + ClinVar + SIGNAL)
			needGenomic := viper.GetBool("annotations.alphamissense") || viper.GetBool("annotations.clinvar") ||
				(viper.GetBool("annotations.signal") && assembly == "grch37")
			if needGenomic {
				dbPath := genomicIndexPath(cacheDir)
				status := "ready"
				ver := fileModDate(dbPath)
				if ver == "" {
					status = "not prepared (run: vibe-vep prepare)"
				}

				if viper.GetBool("annotations.alphamissense") {
					infos = append(infos, sourceInfo{"alphamissense", string(annotate.MatchGenomic), "GRCh38", ver, status,
						[]annotate.ColumnDef{
							{Name: "score", Description: "Pathogenicity score (0-1)"},
							{Name: "class", Description: "likely_benign/ambiguous/likely_pathogenic"},
						}})
				}
				if viper.GetBool("annotations.clinvar") {
					infos = append(infos, sourceInfo{"clinvar", string(annotate.MatchGenomic), "GRCh38", ver, status,
						[]annotate.ColumnDef{
							{Name: "clnsig", Description: "Clinical significance (e.g. Pathogenic, Benign)"},
							{Name: "clnrevstat", Description: "Review status"},
							{Name: "clndn", Description: "Disease name(s)"},
						}})
				}
				if viper.GetBool("annotations.signal") {
					if assembly != "grch37" {
						infos = append(infos, sourceInfo{"signal", string(annotate.MatchGenomic), "GRCh37", "", "skipped (GRCh37 only)", nil})
					} else {
						infos = append(infos, sourceInfo{"signal", string(annotate.MatchGenomic), "GRCh37", ver, status,
							[]annotate.ColumnDef{
								{Name: "mutation_status", Description: "Germline mutation status"},
								{Name: "count_carriers", Description: "Number of carriers in SIGNAL cohort"},
								{Name: "frequency", Description: "Overall allele frequency in SIGNAL cohort"},
							}})
					}
				}
			}

			// Hotspots
			if hsPath := viper.GetString("annotations.hotspots"); hsPath != "" {
				status := "configured"
				ver := fileModDate(hsPath)
				if ver == "" {
					status = "file not found"
				}
				infos = append(infos, sourceInfo{"hotspots", string(annotate.MatchProteinPosition), "any", ver, status,
					[]annotate.ColumnDef{
						{Name: "hotspot", Description: "Y if position is a known cancer hotspot"},
						{Name: "type", Description: "Hotspot type: single residue, in-frame indel, 3d, splice"},
						{Name: "qvalue", Description: "Statistical significance (q-value)"},
					}})
			}

			if len(infos) > 0 {
				fmt.Println()
				fmt.Println("Annotation Sources:")
				tw := tabwriter.NewWriter(os.Stdout, 0, 4, 2, ' ', 0)
				fmt.Fprintln(tw, "  NAME\tMATCH\tASSEMBLY\tVERSION\tSTATUS")
				for _, info := range infos {
					fmt.Fprintf(tw, "  %s\t%s\t%s\t%s\t%s\n", info.name, info.match, info.assembly, info.version, info.status)
				}
				tw.Flush()
			}

			if mafColumns {
				// Parse excluded columns from config
				var excludedCols map[string]bool
				if excl := viper.GetString("exclude-columns"); excl != "" {
					excludedCols = make(map[string]bool)
					for _, s := range strings.Split(excl, ",") {
						s = strings.TrimSpace(s)
						if s != "" {
							excludedCols[s] = true
						}
					}
				}

				fmt.Println()
				fmt.Println("MAF Output Columns:")
				if len(excludedCols) > 0 {
					var names []string
					for k := range excludedCols {
						names = append(names, k)
					}
					fmt.Printf("  Excluded: %s\n", strings.Join(names, ", "))
				}
				fmt.Println()
				tw := tabwriter.NewWriter(os.Stdout, 0, 4, 2, ' ', 0)
				fmt.Fprintln(tw, "COLUMN\tSOURCE\tDESCRIPTION")
				for _, col := range annotate.CoreColumns {
					status := ""
					if excludedCols[col.Name] {
						status = " (excluded)"
					}
					fmt.Fprintf(tw, "vibe.%s\tvibe-vep\t%s%s\n", col.Name, col.Description, status)
				}
				allEffectsStatus := ""
				if excludedCols["all_effects"] {
					allEffectsStatus = " (excluded)"
				}
				fmt.Fprintf(tw, "vibe.all_effects\tvibe-vep\tAll transcript consequences%s\n", allEffectsStatus)
				for _, info := range infos {
					for _, col := range info.columns {
						fmt.Fprintf(tw, "vibe.%s.%s\t%s\t%s\n", info.name, col.Name, info.name, col.Description)
					}
				}
				tw.Flush()
			}

			return nil
		},
	}

	cmd.Flags().BoolVar(&mafColumns, "maf-columns", false, "Show MAF output column mapping")

	return cmd
}
