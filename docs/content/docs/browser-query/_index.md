---
title: Browser Query
weight: 7
---

# Browser-Based Variant Lookup

Query pre-computed variant annotations directly in the browser using [DuckDB-WASM](https://duckdb.org/docs/api/wasm/overview.html). No backend server required — DuckDB makes HTTP range requests against a sorted Parquet file, fetching only the relevant row groups (~20-50 KB per query).

## How It Works

```
vibe-vep export parquet → annotations.parquet → S3/static hosting → DuckDB-WASM → browser
```

1. **Export**: `vibe-vep export parquet` annotates variants and writes a sorted Parquet file
2. **Host**: Upload to S3, GCS, or any static host with CORS and range request support
3. **Query**: DuckDB-WASM in the browser reads only the needed row groups via HTTP range requests

## Try It

The demo below auto-loads a bundled TCGA CHOL example (3,764 variants, 152 KB). You can also point it at your own Parquet file.

{{< parquet-query >}}

## Generating a Parquet File

```bash
# Export a MAF file (canonical transcript, one per variant)
vibe-vep export parquet --canonical --pick -o annotations.parquet input.maf

# Export a VCF file
vibe-vep export parquet -o annotations.parquet input.vcf

# Custom row group size (default: 20,000)
vibe-vep export parquet --row-group-size 10000 -o annotations.parquet input.maf
```

## Hosting Requirements

The Parquet file must be served with:

- **CORS headers**: `Access-Control-Allow-Origin: *` (or your domain)
- **Range request support**: `Accept-Ranges: bytes`
- **Content-Length header**: Required for range requests

Most S3-compatible storage and CDNs support this out of the box. For local testing:

```bash
# npx serve supports range requests
npx serve -l 8080 -C .
```

Note: Python's `http.server` does **not** support range requests.

## JavaScript API

```javascript
import { VibeVEPClient } from './parquet-query.js';

const client = new VibeVEPClient('https://example.com/annotations.parquet');
await client.init();

// Point query by variant
const results = await client.lookupVariant('12', 25245350, 'C', 'A');

// Gene query
const krasVariants = await client.lookupGene('KRAS');

// Free-form SQL
const custom = await client.query(`
  SELECT gene_name, consequence, count(*)
  FROM annotations
  GROUP BY gene_name, consequence
  ORDER BY count(*) DESC
  LIMIT 20
`);

await client.close();
```
