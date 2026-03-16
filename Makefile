.PHONY: build test lint clean install download-testdata download-tcga download-grch37 download-datahub-gdc download-datahub-all docs docs-build wasm wasm-exec parquet-export

# Binary name
BINARY=vibe-vep

# Build flags
LDFLAGS=-ldflags "-s -w"

build:
	go build $(LDFLAGS) -o $(BINARY) ./cmd/vibe-vep

test:
	go test -v -race -cover ./...

lint:
	golangci-lint run ./...

clean:
	rm -f $(BINARY)
	go clean -testcache

install:
	go install $(LDFLAGS) ./cmd/vibe-vep

# Run tests with coverage report
coverage:
	go test -coverprofile=coverage.out ./...
	go tool cover -html=coverage.out -o coverage.html

# Format code
fmt:
	gofmt -s -w .
	goimports -w .

# Download all test data for validation
download-testdata: download-tcga download-grch37

# Download TCGA test data for GRCh38 validation (~1.6GB)
download-tcga:
	./scripts/download_tcga.sh

# Download GRCh37 test data for validation
download-grch37:
	./scripts/download_grch37.sh

# Download datahub GDC studies for GRCh38 validation (~5GB)
download-datahub-gdc:
	./scripts/download_datahub_gdc.sh

# Download all datahub studies for full validation (~50GB)
download-datahub-all:
	./scripts/download_datahub_all.sh

# Download dependencies
deps:
	go mod download
	go mod tidy

# Export annotations to Parquet (override INPUT and OUTPUT as needed)
INPUT ?= testdata/tcga/chol_tcga_gdc_data_mutations.txt
OUTPUT ?= annotations.parquet
parquet-export: build
	./$(BINARY) export parquet --canonical --pick -o $(OUTPUT) $(INPUT)

# Build WASM binary for interactive tutorial
wasm:
	GOOS=js GOARCH=wasm go build $(LDFLAGS) -o docs/static/wasm/vibe-vep.wasm ./cmd/vibe-vep-wasm

# Copy wasm_exec.js from Go installation
wasm-exec:
	cp "$$(go env GOROOT)/lib/wasm/wasm_exec.js" docs/static/js/wasm_exec.js

# Local docs preview (with drafts)
docs:
	cd docs && hugo server -D

# Build docs for production (includes WASM)
docs-build: wasm
	cd docs && hugo --minify
