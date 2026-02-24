.PHONY: build test lint clean install download-testdata

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

# Download TCGA test data for validation (~1.6GB)
download-testdata:
	./scripts/download_tcga.sh

# Download dependencies
deps:
	go mod download
	go mod tidy
