# PolkaPox: Changelog

All notable changes since the initial release are listed below. Whenever possible, links to relevant commits or issues are provided. A full commit history relative to the previous release is also provided. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - 2026-08-28
[Full history](https://github.com/CDCgov/polkapox/pull/80)
Refactoring of tests and writing of tests to ensure 100% coverage

### `Added`
- Added nf-test unit tests for all local modules: `aggregate_tsvs`, `bwa_denovo`, `create_samplesheet`, `graph_reconstruct`, `ivar_consensus_polish_cleanup`, `publish_contigs`, `samplesheet_check`, `samtools_coverage`, `sra_to_samplesheet`, `summarize_qc`, `summarize_tsv`, and `variant_convert` (issue #72).
- Added nf-test for the `filter_reads` subworkflow using publicly available S3 test data.
- Added a dedicated `.github/workflows/nf-test.yml` CI workflow to run local module and subworkflow tests on every push and pull request.
- Added local test data fixtures (FASTQs, BAMs, reference FASTA, samplesheets, MultiQC stats) co-located with each module test.

### `Fixed`
- Migrated `nextflow_schema.json` from JSON Schema draft-07 to draft/2020-12 (`$defs`, updated `$ref` paths) to satisfy `nf-core pipelines lint` 3.5.2.
- Fixed `NfcoreSchema.groovy` to handle both `$defs` (2020-12) and `definitions` (draft-07) schema keys, and to normalise the schema for the Everit library (draft-07 only) before validation. Resolves the `#: could not determine version` CI error.
- Fixed ivar version to avoid frequent segfault errors.
- Fixed `blastn` invocation by removing `makeblastdb` to avoid occasional exception errors.
- Corrected hidden error in listing final-vs-draft assembly path.
- Removed unnecessary assembly FASTA collect step in `summarize_qc`.
- Fixed Nextflow version constraint to be compatible with AWS environments (lower bound `>=23.04.0`).
- Fixed GitHub CI configuration so that test-profile resource settings are correctly respected at runtime.

### `Dependencies`
- Renamed all local module scripts from `<module>.nf` to `main.nf` to follow nf-core conventions.
- Removed a large set of unused nf-core module tests (`tests/modules/nf-core/`) and replaced them with co-located tests under each local module directory.
- Removed tests from `tests/nf-test-data/` since we moved all tests to be local.
- Updated `tests/nextflow.config` to set `errorStrategy = 'terminate'` so test failures are not silently swallowed.

### `Deprecated`
- Removed unused local subworkflows and nf-core module entries that were no longer referenced by the pipeline.
- Removed `subworkflows/local/ivar_consensus_polish` and moved to `ivar_consensus_polish_cleanup` as local module.
- Removed `subworkflows/local/samtools_flagstat_denovo`, thin wrapper no longer needed after refactoring.
- Removed `subworkflows/local/vcftools_diff`


## v1.0dev - 2026-05-21
[Full history]()  
Initial release of PolkaPox for testing, created with the [nf-core](https://nf-co.re/) template.

### `Added`
- Added thin layer subworkflows for de novo BWA MEM, samtools flagstat, iVar consensus polishing, and vcftools diff to preserve legacy PolkaPox behavior while using nf-core modules.

### `Fixed`
- Updated workflow and config syntax for compatibility with current Nextflow versions.
- Fixed de novo BWA MEM wrapper inputs and outputs to match the installed nf-core `bwa/index` and `bwa/mem` module contracts.
- Fixed de novo samtools flagstat and iVar polish version-channel outputs.
- Fixed iVar polish cleanup so the cleaned `*.final.fa` is emitted as a new output file rather than being treated as a staged input.
- Updated `conf/modules.config` selectors for wrapped local subworkflows and removed stale selectors.
- Fixed optional output declarations in affected nf-core modules to avoid `optional()` null-object errors.

### `Dependencies`
- Realigned local workflow wrappers with installed nf-core module interfaces.
- Reduced direct custom edits to nf-core modules by moving PolkaPox-specific behavior into local subworkflows.

### `Deprecated`

## v1.0dev - [date]
[Full history]()  
Initial release of PolkaPox for testing, created with the [nf-core](https://nf-co.re/) template.

### `Added`
- Standard disclaimers and notices from [CDCgov](https://github.com/CDCgov/template) template.

### `Fixed`

### `Dependencies`

### `Deprecated`
