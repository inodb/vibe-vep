<!--
Sync Impact Report
==================
Version change: 0.0.0 → 1.0.0
Bump rationale: Initial constitution adoption (MAJOR)

Modified principles: N/A (initial version)
Added sections:
  - Core Principles (5): Test-First, Simplicity, Performance, Compatibility, Documentation
  - Technical Standards
  - Development Workflow
  - Governance

Templates requiring updates:
  - .specify/templates/plan-template.md: ✅ No updates needed (generic Constitution Check)
  - .specify/templates/spec-template.md: ✅ No updates needed
  - .specify/templates/tasks-template.md: ✅ No updates needed
  - .specify/templates/checklist-template.md: ✅ No updates needed

Follow-up TODOs: None
-->

# vibe-vep Constitution

## Core Principles

### I. Test-First

All functionality MUST be developed using test-driven development (TDD):

- Tests MUST be written before implementation code
- Tests MUST fail before implementation begins (Red phase)
- Implementation MUST make tests pass with minimal code (Green phase)
- Code MUST be refactored only after tests pass (Refactor phase)
- Genomic correctness is critical: variant effect predictions MUST be validated against
  known VEP outputs

**Rationale**: Variant effect prediction errors can lead to incorrect clinical
interpretations. TDD ensures correctness from the start and enables confident refactoring.

### II. Simplicity

vibe-vep MUST remain simple and user-friendly:

- Single static binary distribution with no runtime dependencies
- CLI interface MUST be intuitive with sensible defaults
- Configuration MUST support both flags and config files
- YAGNI: Do not add features until explicitly needed
- Prefer standard library over external dependencies when reasonable
- Start with human/cancer focus; other species are future scope

**Rationale**: The original VEP has significant complexity from Perl dependencies and
installation. A simple Go binary removes these barriers.

### III. Performance

High performance is a primary design goal:

- MUST efficiently reuse existing VEP data files (genome references, cached indexes)
- Smart caching MUST avoid recomputing previously seen variant→protein changes
- Cache hit/miss rates MUST be observable
- Benchmark against original VEP for key operations
- Memory usage MUST be predictable and bounded for large input files
- Parallel processing SHOULD be used where safe and beneficial

**Rationale**: Genomic analysis often involves millions of variants. Performance directly
impacts usability for cancer research workflows.

### IV. Compatibility

vibe-vep MUST maintain compatibility with the VEP ecosystem:

- MUST read existing VEP cache and data files without modification
- Output formats MUST match or be convertible to standard VEP outputs
- MUST support VCF input format as primary input
- Test cases from original VEP SHOULD be adapted for validation
- Breaking changes to data format compatibility require MAJOR version bump

**Rationale**: Users have existing workflows and data investments. Compatibility enables
adoption without data migration.

### V. Documentation

User-facing documentation MUST be comprehensive:

- README MUST include quick-start examples
- All CLI flags and options MUST be documented with `--help`
- Data file requirements and setup MUST be clearly explained
- Error messages MUST be actionable and suggest remediation
- Code comments SHOULD explain "why" not "what"

**Rationale**: User-friendliness is a key differentiator. Good documentation reduces
support burden and enables self-service.

## Technical Standards

### Language and Tooling

- **Language**: Go (latest stable version)
- **Build**: Standard `go build` producing static binary
- **Testing**: `go test` with table-driven tests preferred
- **Linting**: `golangci-lint` with project configuration
- **Formatting**: `gofmt` (enforced)

### Code Organization

- `cmd/`: CLI entry points
- `internal/`: Private packages (core logic)
- `pkg/`: Public packages (if any reusable libraries emerge)
- `testdata/`: Test fixtures including VEP-compatible sample data

### Error Handling

- Errors MUST be wrapped with context using `fmt.Errorf("context: %w", err)`
- User-facing errors MUST be human-readable
- Internal errors SHOULD include file/line context for debugging
- Exit codes MUST follow standard conventions (0=success, 1=error, 2=usage)

## Development Workflow

### Code Review Requirements

- All changes MUST be reviewed before merge
- Tests MUST pass in CI before merge
- Linting MUST pass before merge
- Performance-sensitive changes SHOULD include benchmarks

### Commit Conventions

- Use conventional commits: `feat:`, `fix:`, `docs:`, `test:`, `refactor:`, `perf:`
- Reference issues where applicable
- Keep commits atomic and focused

### Release Process

- Semantic versioning: MAJOR.MINOR.PATCH
- MAJOR: Breaking changes to CLI interface or data format compatibility
- MINOR: New features, backward compatible
- PATCH: Bug fixes, performance improvements, documentation

## Governance

This constitution supersedes all other development practices for vibe-vep.

### Amendment Procedure

1. Propose amendment via pull request to this file
2. Document rationale for change
3. Update version according to semantic versioning:
   - MAJOR: Principle removal or fundamental redefinition
   - MINOR: New principle or significant expansion
   - PATCH: Clarifications, wording improvements
4. Update `LAST_AMENDED_DATE` to amendment date
5. Propagate changes to dependent templates if affected

### Compliance

- All pull requests MUST verify compliance with Core Principles
- Constitution Check in implementation plans MUST reference these principles
- Violations MUST be documented with justification in Complexity Tracking

**Version**: 1.0.0 | **Ratified**: 2026-01-24 | **Last Amended**: 2026-01-24
