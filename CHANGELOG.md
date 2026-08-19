# PolkaPox: Changelog

All notable changes since the initial release are listed below. Whenever possible, links to relevant commits or issues are provided. A full commit history relative to the previous release is also provided. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
- Limit CI nextflow version of version 25

### `Dependencies`

### `Deprecated`
