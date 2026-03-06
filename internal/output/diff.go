package output

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"sort"
	"strings"
)

// DiffWriter compares two sets of variant annotations column by column,
// producing inline TSV output (one row per variant) with L_/R_ prefixed columns.
type DiffWriter struct {
	w           *bufio.Writer
	columns     []string                     // columns to compare (display names)
	leftNames   []string                     // left file column names (parallel to columns)
	rightNames  []string                     // right file column names (parallel to columns)
	matchCount  map[string]int               // column → exact match count
	diffCount   map[string]int               // column → diff count
	transitions map[string]map[[2]string]int // column → [left,right] → count
	total       int
	leftOnly    int
	rightOnly   int
	showAll     bool
	maxDiffs    int // 0 = unlimited
	diffsShown  int

	// Categorizer support
	categorizer *Categorizer
	catCounts   map[string]map[Category]int // column → category → count
	rowCallback func(key string, left, right map[string]string, cats map[string]Category)
}

// NewDiffWriter creates a DiffWriter that compares the given columns.
// leftNames and rightNames are the actual column names in each file (parallel arrays).
// displayNames are used for output headers.
func NewDiffWriter(w io.Writer, displayNames, leftNames, rightNames []string, showAll bool, maxDiffs int) *DiffWriter {
	matchCount := make(map[string]int, len(displayNames))
	diffCount := make(map[string]int, len(displayNames))
	transitions := make(map[string]map[[2]string]int, len(displayNames))
	for _, col := range displayNames {
		transitions[col] = make(map[[2]string]int)
	}
	return &DiffWriter{
		w:           bufio.NewWriter(w),
		columns:     displayNames,
		leftNames:   leftNames,
		rightNames:  rightNames,
		matchCount:  matchCount,
		diffCount:   diffCount,
		transitions: transitions,
		showAll:     showAll,
		maxDiffs:    maxDiffs,
	}
}

// SetCategorizer enables semantic categorization for column comparisons.
func (d *DiffWriter) SetCategorizer(c *Categorizer) {
	d.categorizer = c
	d.catCounts = make(map[string]map[Category]int, len(d.columns))
	for _, col := range d.columns {
		d.catCounts[col] = make(map[Category]int)
	}
}

// SetRowCallback sets a function called for each shared variant when categorizer is set.
func (d *DiffWriter) SetRowCallback(fn func(key string, left, right map[string]string, cats map[string]Category)) {
	d.rowCallback = fn
}

// CategoryCounts returns the per-column category counts (only populated when categorizer is set).
func (d *DiffWriter) CategoryCounts() map[string]map[Category]int {
	return d.catCounts
}

// WriteHeader writes the TSV header with L_/R_ prefixed column names.
// When categorizer is set, adds Cat_<col> columns after each L_/R_ pair.
func (d *DiffWriter) WriteHeader() error {
	parts := []string{"Variant"}
	for _, col := range d.columns {
		parts = append(parts, "L_"+col, "R_"+col)
		if d.categorizer != nil {
			parts = append(parts, "Cat_"+col)
		}
	}
	_, err := fmt.Fprintln(d.w, strings.Join(parts, "\t"))
	return err
}

// WriteDiff compares column values for a matched variant and writes the row.
// left and right map column names (from their respective files) to values.
func (d *DiffWriter) WriteDiff(variantKey string, left, right map[string]string) error {
	d.total++

	if d.categorizer != nil {
		return d.writeDiffCategorized(variantKey, left, right)
	}
	return d.writeDiffSimple(variantKey, left, right)
}

// writeDiffSimple uses simple string equality comparison.
func (d *DiffWriter) writeDiffSimple(variantKey string, left, right map[string]string) error {
	hasDiff := false
	for i, col := range d.columns {
		lv := left[d.leftNames[i]]
		rv := right[d.rightNames[i]]
		if lv == rv {
			d.matchCount[col]++
		} else {
			d.diffCount[col]++
			d.transitions[col][[2]string{lv, rv}]++
			hasDiff = true
		}
	}

	if !hasDiff && !d.showAll {
		return nil
	}

	if !hasDiff {
		return d.writeRow(variantKey, left, right, nil)
	}

	// Has diff — check maxDiffs limit
	if d.maxDiffs > 0 && d.diffsShown >= d.maxDiffs {
		return nil
	}
	d.diffsShown++
	return d.writeRow(variantKey, left, right, nil)
}

// writeDiffCategorized uses the Categorizer for semantic comparison.
func (d *DiffWriter) writeDiffCategorized(variantKey string, left, right map[string]string) error {
	cats := d.categorizer.CategorizeRow(d.columns, left, right, d.leftNames, d.rightNames)

	// Update category counts and transition tracking
	hasDiff := false
	for i, col := range d.columns {
		cat := cats[col]
		d.catCounts[col][cat]++

		lv := left[d.leftNames[i]]
		rv := right[d.rightNames[i]]
		if cat == CatMatch || cat == CatBothEmpty {
			d.matchCount[col]++
		} else {
			d.diffCount[col]++
			d.transitions[col][[2]string{lv, rv}]++
			hasDiff = true
		}
	}

	if d.rowCallback != nil {
		d.rowCallback(variantKey, left, right, cats)
	}

	if !hasDiff && !d.showAll {
		return nil
	}

	if !hasDiff {
		return d.writeRow(variantKey, left, right, cats)
	}

	if d.maxDiffs > 0 && d.diffsShown >= d.maxDiffs {
		return nil
	}
	d.diffsShown++
	return d.writeRow(variantKey, left, right, cats)
}

func (d *DiffWriter) writeRow(variantKey string, left, right map[string]string, cats map[string]Category) error {
	parts := []string{variantKey}
	for i, col := range d.columns {
		parts = append(parts, left[d.leftNames[i]], right[d.rightNames[i]])
		if d.categorizer != nil && cats != nil {
			parts = append(parts, string(cats[col]))
		}
	}
	_, err := fmt.Fprintln(d.w, strings.Join(parts, "\t"))
	return err
}

// WriteLeftOnly records a variant present only in the left file.
func (d *DiffWriter) WriteLeftOnly(variantKey string) {
	d.leftOnly++
}

// WriteRightOnly records a variant present only in the right file.
func (d *DiffWriter) WriteRightOnly(variantKey string) {
	d.rightOnly++
}

// Flush flushes the buffered writer.
func (d *DiffWriter) Flush() error {
	return d.w.Flush()
}

// WriteSummary writes a concordance summary report to w.
func (d *DiffWriter) WriteSummary(w io.Writer, leftFile, rightFile string, leftCount, rightCount int) {
	shared := d.total
	fmt.Fprintf(w, "\nComparison Summary\n")
	fmt.Fprintf(w, "==================\n")
	fmt.Fprintf(w, "Left:  %s (%s variants)\n", leftFile, formatCount(leftCount))
	fmt.Fprintf(w, "Right: %s (%s variants)\n", rightFile, formatCount(rightCount))

	fmt.Fprintf(w, "\nVariant Overlap:\n")
	if leftCount+rightCount > 0 {
		totalUnique := shared + d.leftOnly + d.rightOnly
		if totalUnique > 0 {
			fmt.Fprintf(w, "  Shared:      %s (%.1f%%)\n", formatCount(shared), 100*float64(shared)/float64(totalUnique))
		} else {
			fmt.Fprintf(w, "  Shared:      0\n")
		}
	} else {
		fmt.Fprintf(w, "  Shared:      0\n")
	}
	fmt.Fprintf(w, "  Left only:   %s\n", formatCount(d.leftOnly))
	fmt.Fprintf(w, "  Right only:  %s\n", formatCount(d.rightOnly))

	if shared == 0 {
		return
	}

	fmt.Fprintf(w, "\nColumn Concordance (%s shared variants):\n", formatCount(shared))
	for _, col := range d.columns {
		matches := d.matchCount[col]
		diffs := d.diffCount[col]
		pct := 100 * float64(matches) / float64(shared)
		fmt.Fprintf(w, "  %-30s %5.1f%%  (%s differ)\n", col, pct, formatCount(diffs))
	}

	// Category breakdown (when categorizer is set)
	if d.categorizer != nil && d.catCounts != nil {
		fmt.Fprintf(w, "\nCategory Breakdown (%s shared variants):\n", formatCount(shared))
		for _, col := range d.columns {
			cats := d.catCounts[col]
			if len(cats) == 0 {
				continue
			}
			fmt.Fprintf(w, "\n  %s:\n", col)

			type catCount struct {
				cat   Category
				count int
			}
			var sorted []catCount
			for cat, count := range cats {
				sorted = append(sorted, catCount{cat, count})
			}
			sort.Slice(sorted, func(i, j int) bool {
				if sorted[i].count != sorted[j].count {
					return sorted[i].count > sorted[j].count
				}
				return sorted[i].cat < sorted[j].cat
			})

			for _, cc := range sorted {
				fmt.Fprintf(w, "    %-25s %s\n", cc.cat, formatCount(cc.count))
			}
		}
	}

	// Transition matrices for low-cardinality columns
	for _, col := range d.columns {
		trans := d.transitions[col]
		if len(trans) == 0 {
			continue
		}
		// Count distinct values
		vals := make(map[string]bool)
		for pair := range trans {
			vals[pair[0]] = true
			vals[pair[1]] = true
		}
		if len(vals) > 500 {
			continue // skip high-cardinality columns
		}

		diffs := d.diffCount[col]
		fmt.Fprintf(w, "\n%s Transitions (%s changes):\n", col, formatCount(diffs))

		// Sort transitions by count descending
		type transition struct {
			from, to string
			count    int
		}
		var sorted []transition
		for pair, count := range trans {
			sorted = append(sorted, transition{pair[0], pair[1], count})
		}
		sort.Slice(sorted, func(i, j int) bool {
			if sorted[i].count != sorted[j].count {
				return sorted[i].count > sorted[j].count
			}
			return sorted[i].from < sorted[j].from
		})

		// Show top 20
		limit := 20
		if len(sorted) < limit {
			limit = len(sorted)
		}
		for _, t := range sorted[:limit] {
			from := t.from
			if from == "" {
				from = "(empty)"
			}
			to := t.to
			if to == "" {
				to = "(empty)"
			}
			fmt.Fprintf(w, "  %-40s → %-30s %s\n", from, to, formatCount(t.count))
		}
		if len(sorted) > limit {
			fmt.Fprintf(w, "  ... and %d more\n", len(sorted)-limit)
		}
	}
}

// Stats returns diff statistics.
func (d *DiffWriter) Stats() (total, leftOnly, rightOnly, diffsShown int) {
	return d.total, d.leftOnly, d.rightOnly, d.diffsShown
}

// formatCount formats an integer with comma separators.
func formatCount(n int) string {
	if n < 0 {
		return "-" + formatCount(-n)
	}
	s := fmt.Sprintf("%d", n)
	if len(s) <= 3 {
		return s
	}
	var b strings.Builder
	remainder := len(s) % 3
	if remainder > 0 {
		b.WriteString(s[:remainder])
	}
	for i := remainder; i < len(s); i += 3 {
		if b.Len() > 0 {
			b.WriteByte(',')
		}
		b.WriteString(s[i : i+3])
	}
	return b.String()
}

// NormalizeVariantKey creates a normalized variant key from chrom, pos, ref, alt.
// Strips "chr" prefix and replaces "-" with "" for consistent matching.
func NormalizeVariantKey(chrom string, pos, ref, alt string) string {
	chrom = strings.TrimPrefix(chrom, "chr")
	ref = strings.ReplaceAll(ref, "-", "")
	alt = strings.ReplaceAll(alt, "-", "")
	return chrom + ":" + pos + " " + ref + ">" + alt
}

// ResolveColumns determines which columns to compare between two file headers.
// columns filters to specific column names (empty = all shared).
// colMap provides explicit left→right name mappings.
// Returns parallel arrays of display names, left names, and right names.
func ResolveColumns(leftHeader, rightHeader []string, columns []string, colMap map[string]string) (displayNames, leftNames, rightNames []string, warnings []string) {
	rightSet := make(map[string]bool, len(rightHeader))
	for _, h := range rightHeader {
		rightSet[h] = true
	}
	leftSet := make(map[string]bool, len(leftHeader))
	for _, h := range leftHeader {
		leftSet[h] = true
	}

	// Build auto-matched columns (intersection, preserving left file order)
	var autoMatched []string
	for _, h := range leftHeader {
		if rightSet[h] {
			autoMatched = append(autoMatched, h)
		}
	}

	if len(columns) > 0 {
		// Filter auto-matched to only requested columns
		requested := make(map[string]bool, len(columns))
		for _, c := range columns {
			requested[c] = true
		}
		var filtered []string
		for _, h := range autoMatched {
			if requested[h] {
				filtered = append(filtered, h)
				delete(requested, h)
			}
		}
		autoMatched = filtered

		// Warn about columns that don't exist in either file
		for c := range requested {
			// Check if it's in colMap (might be a mapped column)
			if _, ok := colMap[c]; !ok {
				if !leftSet[c] && !rightSet[c] {
					warnings = append(warnings, fmt.Sprintf("column %q not found in either file", c))
				}
			}
		}
	}

	// Add auto-matched
	for _, h := range autoMatched {
		displayNames = append(displayNames, h)
		leftNames = append(leftNames, h)
		rightNames = append(rightNames, h)
	}

	// Add explicit mappings from colMap
	for left, right := range colMap {
		if leftSet[left] && rightSet[right] {
			displayNames = append(displayNames, left)
			leftNames = append(leftNames, left)
			rightNames = append(rightNames, right)
		} else {
			if !leftSet[left] {
				warnings = append(warnings, fmt.Sprintf("mapped column %q not found in left file", left))
			}
			if !rightSet[right] {
				warnings = append(warnings, fmt.Sprintf("mapped column %q not found in right file", right))
			}
		}
	}

	return
}

// ParseColumnMap parses a --map flag value like "Protein_Change=HGVSp_Short,Effect=Variant_Classification"
// into a map of left_name → right_name.
func ParseColumnMap(mapStr string) (map[string]string, error) {
	if mapStr == "" {
		return nil, nil
	}
	result := make(map[string]string)
	for _, pair := range strings.Split(mapStr, ",") {
		pair = strings.TrimSpace(pair)
		if pair == "" {
			continue
		}
		parts := strings.SplitN(pair, "=", 2)
		if len(parts) != 2 || parts[0] == "" || parts[1] == "" {
			return nil, fmt.Errorf("invalid column mapping %q: expected format left_name=right_name", pair)
		}
		result[strings.TrimSpace(parts[0])] = strings.TrimSpace(parts[1])
	}
	return result, nil
}

// MAF variant key columns — excluded from column comparison by default.
var mafKeyColumns = map[string]bool{
	"Chromosome":       true,
	"Start_Position":   true,
	"End_Position":     true,
	"Reference_Allele": true,
	"Tumor_Seq_Allele2": true,
}

// ReadMAFFile reads a MAF file into a map keyed by normalized variant key.
// Returns the header columns, variant map (key → column values), ordered keys, and total count.
func ReadMAFFile(path string) (header []string, variants map[string][]map[string]string, orderedKeys []string, err error) {
	f, err := openFile(path)
	if err != nil {
		return nil, nil, nil, err
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	// Find header line (skip comments)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") || line == "" {
			continue
		}
		header = strings.Split(line, "\t")
		break
	}
	if header == nil {
		return nil, nil, nil, fmt.Errorf("no header line found in %s", path)
	}

	// Find key column indices
	colIdx := make(map[string]int, len(header))
	for i, h := range header {
		colIdx[h] = i
	}

	chromIdx, ok1 := colIdx["Chromosome"]
	posIdx, ok2 := colIdx["Start_Position"]
	refIdx, ok3 := colIdx["Reference_Allele"]
	altIdx, ok4 := colIdx["Tumor_Seq_Allele2"]
	if !ok1 || !ok2 || !ok3 || !ok4 {
		return nil, nil, nil, fmt.Errorf("missing required key columns (Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) in %s", path)
	}

	variants = make(map[string][]map[string]string)
	keySeen := make(map[string]bool)

	for scanner.Scan() {
		line := scanner.Text()
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) <= chromIdx || len(fields) <= posIdx || len(fields) <= refIdx || len(fields) <= altIdx {
			continue
		}

		key := NormalizeVariantKey(fields[chromIdx], fields[posIdx], fields[refIdx], fields[altIdx])

		row := make(map[string]string, len(header))
		for i, h := range header {
			if i < len(fields) {
				row[h] = fields[i]
			}
		}

		variants[key] = append(variants[key], row)
		if !keySeen[key] {
			orderedKeys = append(orderedKeys, key)
			keySeen[key] = true
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, nil, nil, fmt.Errorf("reading %s: %w", path, err)
	}

	return header, variants, orderedKeys, nil
}

// openFile opens a file, supporting gzipped files.
func openFile(path string) (io.ReadCloser, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}

	// Check for gzip magic bytes
	buf := make([]byte, 2)
	n, err := f.Read(buf)
	if err != nil || n < 2 {
		f.Seek(0, 0)
		return f, nil
	}
	f.Seek(0, 0)

	if buf[0] == 0x1f && buf[1] == 0x8b {
		gz, err := gzip.NewReader(f)
		if err != nil {
			f.Close()
			return nil, err
		}
		return &gzipReadCloser{gz: gz, f: f}, nil
	}

	return f, nil
}

type gzipReadCloser struct {
	gz *gzip.Reader
	f  *os.File
}

func (g *gzipReadCloser) Read(p []byte) (int, error) { return g.gz.Read(p) }
func (g *gzipReadCloser) Close() error {
	g.gz.Close()
	return g.f.Close()
}

// ReadVCFFile reads a VCF file into a map keyed by normalized variant key.
// Returns header column names (from INFO fields), variant map, ordered keys, and total count.
func ReadVCFFile(path string) (header []string, variants map[string][]map[string]string, orderedKeys []string, err error) {
	f, err := openFile(path)
	if err != nil {
		return nil, nil, nil, err
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	// Collect INFO field keys from ##INFO headers and find #CHROM line
	var infoKeys []string
	infoKeySet := make(map[string]bool)
	var vcfCols []string

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "##INFO=<") {
			// Extract ID from ##INFO=<ID=...,
			if idx := strings.Index(line, "ID="); idx >= 0 {
				rest := line[idx+3:]
				if end := strings.IndexByte(rest, ','); end >= 0 {
					key := rest[:end]
					if !infoKeySet[key] {
						infoKeys = append(infoKeys, key)
						infoKeySet[key] = true
					}
				}
			}
			continue
		}
		if strings.HasPrefix(line, "##") {
			continue
		}
		if strings.HasPrefix(line, "#CHROM") {
			vcfCols = strings.Split(line[1:], "\t") // strip leading #
			break
		}
	}

	if vcfCols == nil {
		return nil, nil, nil, fmt.Errorf("no #CHROM header line found in %s", path)
	}

	// Find column indices
	colIdx := make(map[string]int, len(vcfCols))
	for i, h := range vcfCols {
		colIdx[h] = i
	}

	chromIdx, ok1 := colIdx["CHROM"]
	posIdx, ok2 := colIdx["POS"]
	refIdx, ok3 := colIdx["REF"]
	altIdx, ok4 := colIdx["ALT"]
	infoIdx, hasInfo := colIdx["INFO"]
	if !ok1 || !ok2 || !ok3 || !ok4 {
		return nil, nil, nil, fmt.Errorf("missing required VCF columns (CHROM, POS, REF, ALT) in %s", path)
	}

	// Header for comparison = INFO keys + standard VCF columns (QUAL, FILTER)
	header = append([]string{"QUAL", "FILTER"}, infoKeys...)

	variants = make(map[string][]map[string]string)
	keySeen := make(map[string]bool)

	for scanner.Scan() {
		line := scanner.Text()
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) <= altIdx {
			continue
		}

		key := NormalizeVariantKey(fields[chromIdx], fields[posIdx], fields[refIdx], fields[altIdx])

		row := make(map[string]string, len(header))

		// Add QUAL and FILTER
		if qualIdx, ok := colIdx["QUAL"]; ok && qualIdx < len(fields) {
			row["QUAL"] = fields[qualIdx]
		}
		if filterIdx, ok := colIdx["FILTER"]; ok && filterIdx < len(fields) {
			row["FILTER"] = fields[filterIdx]
		}

		// Parse INFO field
		if hasInfo && infoIdx < len(fields) {
			infoStr := fields[infoIdx]
			if infoStr != "." {
				for _, kv := range strings.Split(infoStr, ";") {
					parts := strings.SplitN(kv, "=", 2)
					if len(parts) == 2 {
						row[parts[0]] = parts[1]
					} else {
						row[parts[0]] = "true" // flag fields
					}
				}
			}
		}

		variants[key] = append(variants[key], row)
		if !keySeen[key] {
			orderedKeys = append(orderedKeys, key)
			keySeen[key] = true
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, nil, nil, fmt.Errorf("reading %s: %w", path, err)
	}

	return header, variants, orderedKeys, nil
}

// CompareFiles runs a two-file comparison, writing diffs to stdout and summary to summaryW.
// When categorizer is non-nil, semantic categorization is applied to column comparisons.
func CompareFiles(
	leftHeader, rightHeader []string,
	leftVariants, rightVariants map[string][]map[string]string,
	leftKeys, rightKeys []string,
	leftFile, rightFile string,
	columns []string,
	colMap map[string]string,
	showAll bool,
	maxDiffs int,
	categorizer *Categorizer,
	out io.Writer,
	summaryW io.Writer,
) error {
	// Filter out key columns from headers for auto-matching
	filteredLeft := filterKeyColumns(leftHeader)
	filteredRight := filterKeyColumns(rightHeader)

	displayNames, leftNames, rightNames, warnings := ResolveColumns(filteredLeft, filteredRight, columns, colMap)

	for _, w := range warnings {
		fmt.Fprintf(summaryW, "Warning: %s\n", w)
	}

	if len(displayNames) == 0 {
		return fmt.Errorf("no columns to compare between the two files")
	}

	dw := NewDiffWriter(out, displayNames, leftNames, rightNames, showAll, maxDiffs)
	if categorizer != nil {
		dw.SetCategorizer(categorizer)
	}
	if err := dw.WriteHeader(); err != nil {
		return err
	}

	// Count total variants per file
	leftCount := 0
	for _, rows := range leftVariants {
		leftCount += len(rows)
	}
	rightCount := 0
	for _, rows := range rightVariants {
		rightCount += len(rows)
	}

	// Track which left keys have been matched
	matchedLeft := make(map[string]bool)

	// Process right file variants, looking up in left
	for _, key := range rightKeys {
		rightRows := rightVariants[key]
		leftRows := leftVariants[key]

		for i, rightRow := range rightRows {
			if i < len(leftRows) {
				matchedLeft[key] = true
				if err := dw.WriteDiff(key, leftRows[i], rightRow); err != nil {
					return err
				}
			} else {
				dw.WriteRightOnly(key)
			}
		}
	}

	// Report left-only variants
	for _, key := range leftKeys {
		leftRows := leftVariants[key]
		rightRows := rightVariants[key]
		if len(rightRows) == 0 {
			for range leftRows {
				dw.WriteLeftOnly(key)
			}
		} else if len(leftRows) > len(rightRows) {
			// Extra rows on the left side
			for i := len(rightRows); i < len(leftRows); i++ {
				dw.WriteLeftOnly(key)
			}
		}
	}

	if err := dw.Flush(); err != nil {
		return err
	}

	dw.WriteSummary(summaryW, leftFile, rightFile, leftCount, rightCount)

	return nil
}

// filterKeyColumns removes MAF key columns from a header list.
func filterKeyColumns(header []string) []string {
	var filtered []string
	for _, h := range header {
		if !mafKeyColumns[h] {
			filtered = append(filtered, h)
		}
	}
	return filtered
}
