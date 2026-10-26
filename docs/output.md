# PolkaPox: Output

## Introduction

This document describes the output produced by the pipeline. Sub-directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data into the following outputs:

- [fq lint](#fq-lint) - Validation of the input FASTQ files
- [Kraken2](#kraken2) - Fastq sequence classification 
- [KrakenTools](#krakentools) - Keep only sequences classified as Orthopox Virus
- [FastP](#fastp) - Raw read QC and adapter + quality trimming
- [BWA mem](#bwa) - Read alignment to a reference genome
- [samtools](#samtools) - Read alignment depth + metrics
- [iVar](#ivar) - Consensus sequence generation and variant calling from reference-based assembly
- [Variant summaries](#variant-summaries) - Summary table of variant calls for select coordinates if `coords` parameter is specified
- [Unicycler](#unicycler) - De novo assembly from Fastp-trimmed reads
- [Bandage](#bandage) - Assembly graph image
- [QUAST](#quast) - Assembly quality
- [Graph_recon](#graph-reconstruction) - Assembly graph resolution
- [MUMmer](#mummer-dnadiff) - Quantify assembly corrections
- [Final assembly](#final-assembly) - Polish de novo assembly
- [summarize_qc.py](#qc-summary) - Aggregate summary metrics in a tsv file per sample
- [MultiQC](#multiqc) - Aggregate report describing QC related to read quality
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### fq lint

<details markdown="1">
<summary>Output files</summary>

- `fq_lint/`
  - `*.fq_lint.txt`: `fq lint` log for each samplesheet entry (paired-end files are validated together).

</details>

The pipeline stops if a FASTQ is invalid, a pair has mismatched read names, or a file contains no reads. See [usage](usage.md#fastq-validation).

[fq documentation](https://github.com/stjude-rust-labs/fq)

### Kraken2

<details markdown="1">
<summary>Output files</summary>

- `kraken2/`
  - `*.classified_{1,2}.fastq.gz`: Classified read pairs (single-end: `*.classified.fastq.gz`).
  - `*.unclassified_{1,2}.fastq.gz`: Unclassified read pairs (single-end: `*.unclassified.fastq.gz`).
  - `*.kraken2.classifiedreads.txt`: Report showing classification output for each read.
  - `*.kraken2.report.txt`: Taxonomic classification report.

</details>

[Kraken2 documentation](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)

### KrakenTools

<details markdown="1">
<summary>Output files</summary>

- `krakentools/`
  - `*.extracted_kraken2_read_{1,2}.fastq.gz`: Fastqs containing only sequences classified as Orthopox Virus with Kraken2 (single-end: `*.extracted_kraken2_read.fastq.gz`).

</details>

[KrakenTools documentation](https://github.com/jenniferlu717/KrakenTools)

### FastP

<details markdown="1">
<summary>Output files</summary>

- `fastp/`
  - `*_R1.fastp.fastq.gz`, `*_R2.fastp.fastq.gz`: Trimmed read output for paired-end input (single-end: `*.fastp.fastq.gz`).
  - `*.fastp.html`: QC report containing quality metrics.
  - `*.fastp.json`: QC report in json file.
  - `*.fastp.log`: fastp log.

</details>

[Fastp documentation](https://github.com/OpenGene/fastp)

### BWA mem

<details markdown="1">
<summary>Output files</summary>

- `bwa/`
  - `*.bam`: Sorted bam alignment to reference genome.
  - `bwa/*`: Reference genome index files

</details>

[BWA mem documentation](http://bio-bwa.sourceforge.net/bwa.shtml)

### samtools

<details markdown="1">
<summary>Output files</summary>

- `samtools/`
  - `*.bam.bai`: Index of the reference alignment in `bwa/`.
  - `*.depth.tsv`: Coverage depth vs reference genome.
  - `*.flagstat`: Alignment stats from `samtools flagstat`

</details>

[samtools documentation](http://www.htslib.org/)

### iVar

<details markdown="1">
<summary>Output files</summary>

- `ivar/`
  - `*.consensus.fa`: Consensus generated from BWA MEM alignment to reference genome.
  - `*.final.fa`: Intermediate de novo assembly following read mapping and correction but before the `Consensus_` header prefix and terminal Ns are removed (the cleaned file is in `final_assembly/`).
  - `*.bwa.mpileup`: Mpileup output from BWA MEM alignment to reference genome.
  - `*.bwa.qual.txt`, `*.final.qual.txt`: Per-position quality output from `ivar consensus` for the reference-based and de novo consensus.
  - `*.ivar.tsv`: Default ivar variant output with variants as tsv table for each sample.
  - `*.ivar.vcf`: VCF converted from ivar variants tsv for each sample.

</details>

[iVar documentation](https://andersen-lab.github.io/ivar/html/)

### Variant summaries

<details markdown="1">
<summary>Output files</summary>

- `variant_summaries/`
  - `all_samples.vcf.summary.txt`: Concatenation of the sample-level summary files for regions of interest specified by the `coords` parameter.
  - `*_ivar_summary.txt`: Sample-level summary files for regions of interest specified by the `coords` parameter.

</details>

### Unicycler

<details markdown="1">
<summary>Output files</summary>

- `unicycler/`
  - `*.bridges_applied.gfa`: Assembly graph output.
  - `*.scaffolds.fa.gz`: Scaffold-level de novo genome assembly output.
  - `*.unicycler.log`: Log of unicycler process.

</details>

Unicycler is pinned to v0.4.8 (local module) and runs without long-read input.

[Unicycler documentation](https://github.com/rrwick/Unicycler)

### Bandage

<details markdown="1">
<summary>Output files</summary>

- `bandage/`
  - `*.png`, `*.svg`: Image of the Unicycler assembly graph.

</details>

[Bandage documentation](https://rrwick.github.io/Bandage/)

### QUAST

<details markdown="1">
<summary>Output files</summary>

- `quast/`
  - `<sample>/`: QUAST output directory for each sample (contains `report.tsv`, which MultiQC parses).
  - `<sample>.tsv`: Link to the QUAST `report.tsv` for each sample.

</details>

[QUAST documentation](https://quast.sourceforge.net/)

### Graph reconstruction

<details markdown="1">
<summary>Output files</summary>

- `graph_recon/`
  - `*.assembly_asm.fasta`: Assembled genome.
  - `*.assembly.summary`: Summary QC metrics for graph reconstruction. 
  - `*.assembly.log`: Log of assembly graph reconstruction process.
  - `*contigs.fasta`: Reformatted fasta file from the input gfa file.
- `graph_recon_mapping/`
  - `*.bam`, `*.bam.bai`: Sorted bam alignment of reads to reconstructed genome and its index.
  - `*.denovo.flagstat`: Alignment stats from `samtools flagstat`
  - `*.denovo.coverage`: Per-contig coverage from `samtools coverage`

</details>

[documentation](/bin/)

### MUMmer dnadiff

<details markdown="1">
<summary>Output files</summary>

- `mummer/`
  - `*.report`: Summary of assembly corrections detected with `dnadiff` (reconstructed graph assembly vs. the polished assembly). `TotalSNPs` and `TotalIndels` are reported as `corrected_snps` and `corrected_indels` in `sample_summary.tsv`.

</details>

[MUMmer documentation](https://github.com/mummer4/mummer)

### Final assembly

<details markdown="1">
<summary>Output files</summary>

- `final_assembly/`
  - `*.final.fa`: Final assembly generated by de novo subworkflow
  - `*.draft.fa`: Multi-contig draft assembly written only if Graph_Recon fails to reconstruct the Unicycler graph

</details>

### QC Summary

<details markdown="1">
<summary>Output files</summary>

- `sample_summary.tsv`: A tsv file with fields summarizing QC metrics for various steps in the pipeline. Note that columns included varies with subworkflow executed. All columns will be present with `--workflow full` and `--filter true`, but other options will only include relevant outputs.
  - Column summary:
    1. `sample` - Sample name
    2. `reference_genome` - Reference genome defined with `--fasta`
    3. `total_raw_reads` - Total count of raw input reads
    4. `opx_read_count_kraken` - Count of reads classified as orthopox with kraken
    5. `opx_percent_kraken` - Percentage of total reads classified as orthopox
    6. `human_percent_kraken` - Percentage of reads classified as human
    7. `unclass_percent_kraken` - Percentage of reads not classified
    8. `kraken_db` - kraken database defined with `--kraken_db`
    9. `kraken_tax_ids` - List of taxids used for classification with `--kraken2_tax_ids`
    10. `filtered_read_count_fastp` - Count of classified reads passing fastp trimming and filtering
    11. `percent_reads_passed_fastp` - Percentage of classified reads passing fastp
    12. `percent_adapter_fastp` - Percent of classified reads with adapter contamination
    13. `gc_content_postfilter_fastp` - Percent QC of filtered reads
    14. `q30_rate_postfilter_fastp` - Q30 of filtered reads
    15. `percent_duplication_fastp` - Percent duplication in filtered reads
    16. `reads_mapped_bwa` - Count of filtered reads mapping to the reference genome
    17. `percent_mapped_bwa` - Percent of filtered reads mapping to the reference genome
    18. `average_depth_bwa` - Average coverage depth of filtered reads mapped to the reference genome
    19. `count_20xdepth_bwa` - Count of reference positions with >20x coverage (breadth 20x)
    20. `n_contigs_unicycler` - Number of contigs produced by de novo assembly with Unicycler
    21. `assembly_length_unicycler` - Total length of assembled contigs
    22. `n50_unicycler` - N50 of Unicycler assembly
    23. `mapped_reads_denovo` - Count of filtered reads mapped to the final assembly
    24. `percent_mapped_denovo` - Percent of filtered reads mapped to the final assembly
    25. `mean_depth_denovo` - Mean coverage depth of filtered reads mapped to the final assembly
    26. `orientation_copy_number` - Final assembly graph path
    27. `sequence_length` - Total length of final assembly
    28. `itr_length` - ITR length inferred from the assembly graph
    29. `gfa_status` - Status of assembly graph resolution (PASS|FAIL)
    30. `gfa_notes` - Detailed message of graph resolution status
    31. `corrected_snps` - Count of SNPs corrected when polishing the final assembly (`TotalSNPs` of the `dnadiff` report in `mummer/`)
    32. `corrected_indels` - Count of indels corrected when polishing the final assembly (`TotalIndels` of the `dnadiff` report in `mummer/`)
    33. `corrected_Ns` - Count of Ns added to the final assembly due to low coverage
    34. `total_snps` - Count of variants called by iVar relative to the reference genome (depth >= `--DP_cutoff`, allele frequency >= `--af_cutoff`)
    35. `COORD1_SNPs` - Count of filtered SNPs within defined coordinate range (column named after `--locus1_name`, e.g. `F13L_SNPs`)
    36. `COORD2_SNPs` - Count of filtered SNPs within defined coordinate range (column named after `--locus2_name`, e.g. `E9L_SNPs`)
    37. `final_assembly` - Path to the final assembly (`final_assembly/*.final.fa`, else `*.draft.fa`, else the reference-based `ivar/*.consensus.fa`)
    38. `opxv_reads_1` - Path to the KrakenTools-extracted orthopox reads, read 1
    39. `opxv_reads_2` - Path to the KrakenTools-extracted orthopox reads, read 2

Kraken columns are computed from `kraken2/*.kraken2.classifiedreads.txt` and the number of reads in `krakentools/*.extracted_kraken2_read_1.fastq.gz`. Values are `NA` when the corresponding file is missing.
 
</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `samplesheet.valid.csv`: Reformatted samplesheet files used as input to the pipeline.
  - `collated_versions.yml`: Software versions of every process, collected from the pipeline's `versions` topic channel, plus the workflow and Nextflow versions.
  - `pipeline_report_*`: Reports generated with `--email` / `--email_on_fail` parameters at runtime.
  - `execution_trace_*`: Tracing files that contain information about each process executed in the pipeline.
  - `pipeline_dag_`: Graphical representations of the directed acyclic graph corresponding to the workflow structure defined in the pipeline.
  - `execution_timeline_`: Files that contain information about the execution timeline of tasks completed by the pipeline.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
