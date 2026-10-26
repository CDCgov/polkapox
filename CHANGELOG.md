# PolkaPox: Changelog

All notable changes since the initial release are listed below. Whenever possible, links to relevant commits or issues are provided. A full commit history relative to the previous release is also provided. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v3.1.2 - 2026-10-23
[Full history](https://github.com/CDCgov/polkapox/pull/79)
Upgrade of all nf-core modules to their latest versions, migration of software version tracking to Nextflow topic channels, and FASTQ validation of the inputs

### `Added`
- Added the nf-core `fq/lint` module to `INPUT_CHECK` to validate every samplesheet entry (paired and single-end). Reads are released downstream only after validation passes, empty FASTQs raise an error, and reports are published to `<outdir>/fq_lint/`. Documented in `docs/usage.md`.
- Added nf-core `krakentools/extractkrakenreads` (run with `--fastq-output --include-children`) to extract orthopox reads after Kraken2, replacing `seqtk/subseq`.
- Added nf-core `samtools/index` and `samtools/faidx` modules, and the `custom/sratoolsncbisettings` module.
- Added `meta.yml` files for all local modules.
- Added nf-tests for the `fq/lint` empty-FASTQ handling in `input_check`, and nf-core module tests for the updated modules.
- Added a workflow version block (`getWorkflowVersion` / `workflowVersionToYAML`) in `workflows/polkapox.nf` for the collated software versions report.
- Added `subworkflows/local/ivar_consensus_polish`, a thin wrapper around nf-core `ivar/consensus` that keeps the legacy de novo polishing interface.

### `Fixed`
- Updated `bin/summarize_qc.py` to read the current MultiQC general-stats headers and convert the new display units (millions of reads, Kbp/Mbp) back to the units of the older summary.
- Updated `conf/modules.config` to the new nf-core argument layout for `IVAR_CONSENSUS`, `IVAR_VARIANTS` and `SAMTOOLS_DEPTH` (swapped `ext.args`/`ext.args2`), and renamed the `BANDAGE` selector to `BANDAGE_IMAGE`.
- Changed the de novo polished consensus output prefix from `.polished` to `.final`, and set the `.ivar` and `.depth` prefixes that `summarize_qc.py` expects.
- Fixed `multiqc_config` handling to use a checked file path (`checkIfExists`) rather than an optional channel.
- Fixed `SRA_TOOLS` samplesheet generation to use the new combined `reads` output for single and paired-end data.
- Fixed `subworkflows/local/ref_based` test: `BWA_INDEX` setup now receives a `[meta, fasta]` tuple.
- Fixed the `denovo` and `filter_reads` tests: removed the `versions` assertion (no longer emitted), and `filter_reads` now snapshots `orthopox_reads` instead of the non-existent `seqtk_reads`.
- Refreshed the `denovo`, `filter_reads` and other snapshots for the new QUAST TSV output, `.final.mpileup` naming, KrakenTools output and fastp read names.
- Fixed `bin/summarize_qc.py` for the upgraded modules: Kraken percentages and `total_raw_reads` now count records in the KrakenTools-extracted FASTQ (previously the removed `*.opxreads.txt`), `corrected_Ns` reads `*.final.mpileup`, and `opxv_reads_*` point to `krakentools/`. Added a `summarize_qc` nf-test covering these inputs.
- Restored the iVar `-m ${params.DP_cutoff} -t ${params.af_cutoff}` thresholds in `conf/modules.config` (dropped when moving to the nf-core `ivar/variants` module), so `total_snps` is again filtered at depth >= 5 and allele frequency >= 0.01 instead of iVar defaults.
- Updated `docs/output.md` to the current published outputs: new `fq_lint/`, `bandage/` and `krakentools/` sections; Kraken2/fastp/Unicycler/QUAST/graph reconstruction file names; `collated_versions.yml` instead of `software_versions.yml`; and the full 39-column `sample_summary.tsv` list (adds `mean_depth_denovo`, `final_assembly`, `opxv_reads_1/2`, corrected column order).
- Updated `docs/usage.md` (KrakenTools read filtering, iVar `--DP_cutoff` / `--af_cutoff` / `--altDP_cutoff` thresholds, and `--genome` -> `--fasta` in the run commands) and `README.md` (fq lint and KrakenTools steps, MUMmer `dnadiff`, output directory tree, minimum Nextflow version). Replaced seqtk with KrakenTools and added fq in `CITATIONS.md`.
- Added a local `DNADIFF` module (MUMmer 3.23 `dnadiff`, with nf-test) and used it in the `DENOVO` subworkflow instead of nf-core `mummer`. The updated nf-core module only writes a `.coords` file, so `corrected_snps` / `corrected_indels` in `sample_summary.tsv` had become `NA`. The `*.report` is still published to `mummer/`.
- Set `BWA_INDEX` memory to 6 GB in `tests/nextflow.config`. The module requests 7 B per FASTA byte (~1.4 MB for the test genome), which is below Docker's 6 MB minimum and failed `--profile docker` runs.
- Disabled the fq `P001` paired-name validator in `conf/test.config`: the S3 test reads carry a trailing `kraken:taxid|...` token, so the `1:N:0` / `2:N:0` mate comments no longer matched and `FQ_LINT` failed under `-profile test`.

### `Dependencies`
- Updated all nf-core modules to their latest versions: `bandage/image`, `bwa/index`, `bwa/mem`, `fastp`, `ivar/consensus`, `ivar/variants`, `kraken2/kraken2`, `multiqc`, `mummer`, `quast`, `samtools/depth`, `samtools/flagstat`, `sratools/prefetch`, `sratools/fasterqdump` and `vcftools`.
- Moved software version reporting to Nextflow topic channels (`topic: versions`) in all local modules and subworkflows, and removed the `ch_versions` plumbing from `workflows/polkapox.nf`.
- Local `unicycler` remains pinned to v0.4.8 and `summarize_qc` intentionally does not publish to the versions topic, to avoid a circular wait with MultiQC.
- `BWA_INDEX` is now called directly from `workflows/polkapox.nf` with `[ [id:'reference'], fasta ]`.
- Raised the minimum Nextflow version to `25.04.0` (`manifest.nextflowVersion` and README), which the `topic: versions` outputs of the updated modules require, and dropped the `23.04.0` job from the CI matrix. The supported range is `>=25.04.0, <=25.10.4`.

### `Deprecated`
- Workflow no longer suports NExtflow 23.04
- Removed `CUSTOM_DUMPSOFTWAREVERSIONS` (module, template and config selectors) in favour of topic-channel version collation.
- Removed `subworkflows/local/prepare_genome` and the nf-core `seqtk/subseq` module.
- Removed `modules/nf-core/sratools/fasterqdump/nextflow.config`.


## v3.1.1 - 2026-09-28
[Full history](https://github.com/CDCgov/polkapox/pull/62)
Refactoring of tests and writing of tests to ensure 100% coverage

### `Added`
- Added nf-test unit tests for all local modules: `aggregate_tsvs`, `bwa_denovo`, `create_samplesheet`, `graph_reconstruct`, `ivar_consensus_polish_cleanup`, `publish_contigs`, `samplesheet_check`, `samtools_coverage`, `sra_to_samplesheet`, `summarize_qc`, `summarize_tsv`, and `variant_convert` (issue #72).
- Added nf-test for the `filter_reads` subworkflow using publicly available S3 test data.
- Added nf-tests for the `denovo`, `input_check`, `prepare_genome`, and `ref_based` local subworkflows and wired them into the nf-test CI workflow.
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
- Capped CPUs and memory for all test processes (2 CPUs, 6 GB) via `tests/nextflow.config` so nf-test CI runs fit within GitHub Actions runners.
- Excluded `versions.yml` from module and subworkflow test snapshots, since its embedded tool-version strings differ between environments and caused spurious snapshot mismatches.
- Excluded non-deterministic outputs (`graph_reconstruct` log, MUMmer `.report`) from snapshots and replaced them with existence/size assertions.
- Made `AssemblyGraph_gfaPy.py` ITR selection deterministic (`sorted(set(...))`) to stabilise `graph_reconstruct` snapshots.
- Skipped `QUAST` in the test profile (`ext.when = false`) because `quast.py` is unavailable in the nf-test environment.
- Removed the stale `CUSTOM_DUMPSOFTWAREVERSIONS` selector from `conf/modules.config` to silence the "no process matching config selector" warning.

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
