# PolkaPox: Output

## Introduction

This document describes the output produced by the pipeline. Sub-directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data into the following outputs:

- [Kraken2](#kraken2) - Fastq sequence classification 
- [seqTK](#seqtk) - Keep only sequences classified as Monkeypox Virus
- [FastP](#fastp) - Raw read QC and adapter + quality trimming
- [BWA mem](#bwa) - Read alignment to a reference genome
- [samtools](#samtools) - Read alignment depth + metrics
- [iVar](#ivar) - Consensus sequence generation and variant calling from reference-based assembly
- [Variant summaries](#variant-summaries) - Summary table of variant calls for select coordinates if `coords` parameter is specified
- [Unicycler](#unicycler) - De novo assembly from Fastp-trimmed reads
- [QUAST](#quast) - Assembly quality
- [Graph_recon](#graph) - Assembly graph resolution
- [MUMmer](#mummer) - Quantify assembly corrections
- [Final assembly](#final-assembly) - Polish de novo assembly
- [summarize_qc.py](#qc-summary) - Aggregate summary metrics in a tsv file per sample
- [MultiQC](#multiqc) - Aggregate report describing QC related to read quality
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Kraken2

<details markdown="1">
<summary>Output files</summary>

- `kraken2/`
  - `*.classified.fastq.gz`: Classified read pairs.
  - `*.classifiedreads.txt`: Report showing classification output for each read.
  - `*.report.txt`: Taxonomic classification report.

</details>

[Kraken2 documentation](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)

### seqTK

<details markdown="1">
<summary>Output files</summary>

- `seqtk/`
  - `*.fq.gz`: Fastqs containing only sequences classified as Monkeypox Virus with Kraken2.

</details>

[seqTK documentation](https://github.com/lh3/seqtk)

### FastP

<details markdown="1">
<summary>Output files</summary>

- `fastp/`
  - `*.fastp.fastq.gz`: Trimmed read output. If paired input, a trimmed file will exist for R1 and R2.
  - `*.fastp.html`: QC report containing quality metrics.
  - `*.fastp.json`: QC report in json file.
  - `*.fastp.txt`: QC report in text file.

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
  - `*.depth.tsv`: Coverage depth vs reference genome.
  - `*.flagstat`: Alignment stats from `samtools flagstat`

</details>

[samtools documentation](http://www.htslib.org/)

### iVar

<details markdown="1">
<summary>Output files</summary>

- `ivar/`
  - `*.bwa.fa`: Consensus generated from BWA MEM alignment to reference genome.
  - `*.bwa.mpileup*`: Mpileup output from BWA MEM alignment to reference genome.
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
  - `*.assembly.gfa.gz`: Assembly graph output.
  - `*.scaffolds.fa.gz`: Scaffold-level de novo genome assembly output.
  - `*.unicycler.log`: Log of unicycler process.

</details>

[Unicycler documentation](https://github.com/rrwick/Unicycler)

### QUAST

<details markdown="1">
<summary>Output files</summary>

- `quast/`
  - `quast/*`: QUAST output files.
  - `report.tsv`: Summary of assembly metrics for all samples.

</details>

[QUAST documentation](https://quast.sourceforge.net/)

### Graph reconstruction

<details markdown="1">
<summary>Output files</summary>

- `graph_recon/`
  - `*.assembly_asm.fasta`: Assembled genome.
  - `*.assembly_longest.fasta`: Longest input contig sequence from Unicycler
  - `*.assembly.summary`: Summary QC metrics for graph reconstruction. 
  - `*.assembly.log`: Log of assembly graph reconstruction process.
  - `*contigs.fasta`: Reformatted fasta file from the input gfa file.
- `graph_recon_mapping/`
  - `*.bam`: Sorted bam alignment of reads to reconstructed genome.
  - `*flagstat`: Alignment stats from `samtools flagstat`
  - `bwa/*`: Reconstructed genome index files

</details>

[documentation](/bin/)

### MUMmer dnadiff

<details markdown="1">
<summary>Output files</summary>

- `mummer/`
  - `*.report`: Summary of assembly corrections detected with `dnadiff`.

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
    25. `orientation_copy_number` - Final assembly graph path
    26. `sequence_length` - Total length of final assembly
    27. `itr_length` - ITR length inferred from the assembly graph
    28. `gfa_status` - Status of assembly graph resolution (PASS|FAIL)
    29. `gfa_notes` - Detailed message of graph resolution status
    30. `total_snps` - Count of detected SNPs relative to the reference genome
    31. `COORD1_SNPs` - Count of filtered SNPs within defined coordinate range
    32. `COORD2_SNPs` - Count of filtered SNPs within defined coordinate range
    33. `corrected_snps` - Count of SNPs corrected when polishing the final assembly
    34. `corrected_indels` - Count of indels corrected when polishing the final assembly
    35. `corrected_Ns` - Count of Ns added to the final assembly due to low coverage
 
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
  - `software_versions.yml`: Captures all software versions used within the workflow, sub-workflows, and modules.
  - `pipeline_report_*`: Reports generated with `--email` / `--email_on_fail` parameters at runtime.
  - `execution_trace_*`: Tracing files that contain information about each process executed in the pipeline.
  - `pipeline_dag_`: Graphical representations of the directed acyclic graph corresponding to the workflow structure defined in the pipeline.
  - `execution_timeline_`: Files that contain information about the execution timeline of tasks completed by the pipeline.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
