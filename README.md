# PolkaPox :microbe: :dna: :accordion:

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)

## Introduction

**PolkaPox** is a workflow for filtering, trimming, QC, reference-based analysis, and de novo assembly of **Illumina short** sequencing reads from orthopoxviruses. 

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

**Read pre-processing**
1. Filter raw reads to capture only Orthopox Virus sequences. ([`Kraken`](https://github.com/DerrickWood/kraken2) + [`SeqTK`](https://github.com/lh3/seqtk))
1. Trim and filter reads for adapter contamination and quality, summarize QC ([`Fastp`](https://github.com/OpenGene/fastp))

**Reference-based analyses**
1. Align reads to reference genome. ([`BWA`](http://bio-bwa.sourceforge.net/))
1. Variant calling and target gene mutation summary. ([`iVar`](https://andersen-lab.github.io/ivar/html/manualpage.html))
1. Consensus sequence generation. ([`iVar`](https://andersen-lab.github.io/ivar/html/manualpage.html))

**Reference-free analyses**
1. De novo assembly. ([`Unicycler`](https://github.com/rrwick/Unicycler))
1. Calculate assembly quality metrics. ([`QUAST`](http://quast.sourceforge.net/))
1. Assembly graph resolution. ([`AssemblyGraph_gfaPy.py`](/bin/AssemblyGraph_gfaPy.py))
1. Align reads to assembled genome. ([`BWA`](http://bio-bwa.sourceforge.net/))
1. Correct assembly errors and ambiguities. ([`iVar`](https://andersen-lab.github.io/ivar/html/manualpage.html))
1. Quantify assembly corrections. ([`MUMmer`](https://github.com/mummer4/mummer#dnadiff))

**Performance summary**
1. Visualize reads metrics and summary of software versions. ([`MultiQC`](http://multiqc.info/))
1. Compile various QC metrics. ([`summarize_qc.py`](/bin/summarize_qc.py))

![workflow diagram](/docs/images/polkapox_workflow.png)

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.
   - *Note*: if running this pipeline on CDC infrascructure (aspen, biolinux), use singularity containers.

3. Clone this repo to your working environment:

   ```console
   git clone git@github.com:CDCgov/polkapox.git
   ```

4. Prepare your samples. You have two options:
> (A) Make a sample sheet which will act as pipeline input. The samplesheet should consist of three columns with sample ID, R1 and R2 specified. If you are using single-end data, only one fastq can be specified and the pipeline will auto detect this. Example:

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
``` 
> (B) Pass a directory of FastQ files as your input and the pipeline will auto-create the samplesheet. See [Usage](../docs/usage.md) for more details.

5. Start running your own analysis!

   ```console
   nextflow run polkapox/main.nf --input {SAMPLESHEET.csv OR input_directory} --outdir {OUTDIR} --fasta {REF.fa} -profile sge,singularity --kraken_db {PATH/TO/DB} --gff {ANNOTATION.gff} --workflow {WORKFLOW} --filter {true/false}
   ```
   
   **note**: If you do not provide `--fasta`, `--gff`, or `--kraken_db`, they will default to a Clade II reference and gff in the `assets` folder of this repo, and a kraken db will be downloaded from `s3://io-pe1-prod-ncezid-oamd-nextstrain/polkapox/orthopox_kdb/`. If you do not specify `--filter` then it will default to `true`. See `nextflow.config` for details.  Add `--file_levels {top (default)/nested}` if passing a directory as input. See [usage](/docs/usage.md) for details.

## Pipeline configuration

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

## Running on CDC cluster

For use on an HPC cluster (rosalind, aspen) the workflow can be run as a job by the following steps:

1) After logging in, activate the nextflow environment:

    ```console
    module load nextflow
    ```

2) Point to the `cdc.config` file, which contains custom profiles for the CDC HPC clusters. Submit individual processes as jobs to the scheduler using a profile defined in `cdc.config`.  For example, to run a job on rosalind:

    ```consol
    nextflow run main.nf --input {SAMPLESHEET.csv} --outdir {OUTDIR} --fasta {REF.fa} -profile rosalind,singularity --kraken_db {PATH/TO/DB} --gff {ANNOTATION.gff} -config /scicomp/reference/nextflow/configs/cdc.config
    ```

## Documentation
### Basic usage
The **PolkaPox** pipeline requires two inputs: [1] a samplesheet (or input directory) `--input` and [2] a `--workflow` definition. Additional details are provided in [usage](/docs/usage.md).  

### Output
Pipeline outputs are organized into sub-directories for each step in the selected workflow. All paths are relative to the top-level results directory `--outdir` and additional details are provided in [output](/docs/output.md).  

```
${outdir}/
  ├── bwa
  ├── bandage
  ├── fastp
  ├── final_assembly
  ├── graph_recon
  ├── graph_recon_mapping
  ├── ivar
  ├── ivar_variants
  ├── kraken2
  ├── multiqc
  ├── mummer
  ├── pipeline_info
  ├── quast
  ├── sample_summary.tsv
  ├── samtools
  ├── seqtk
  ├── unicycler
  └── variant_summaries
```
### De novo assembly
The **PolkaPox** pipeline includes de novo assembly optimized for the linear genome architecture of orthopoxviruses. Additional details are provided in [de novo](/docs/denovo.md).

## Credits

Contributors:\
Kyle O'Connell | Michael Weigand | Jessica Rowell | Shatavia Morrison\
Kristen Knipe | Ethan Hetrick | Crystal Gigante | Lynsey Kovar\
Hunter Seabolt | Dhwani Batra | Daisy McGrath\
Yesh Kulasekarapandian | Jason Caravas

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
---
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/polkapox/blob/master/DISCLAIMER.md) 
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).

## Repository Metadata
```
Organization: NCEZID-OAMD
contact email: ncezid_shareit@cdc.gov
exemption status: NA
exemption justification: NA
description fields: Nextflow workflow for the assembly of Orthopox virus sequences
```
