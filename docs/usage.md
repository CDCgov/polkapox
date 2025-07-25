# PolkaPox: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

## Samplesheet or input directory input

*You will need to either provide a samplesheet, or a path to a directory with your fastq files*

> A) Provide a samplesheet: You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```console
--input '[path to samplesheet file]'
```

> B) Provide a path to FastQ files: You will need to pass the path to a folder containing your samples, as shown below: 

```console
--indir '[path to fastq files]' --file_levels '[top (default)/nested]'
```

> For this option, PolkaPox currently supports two file organization structures. 
> Option 1: Your FastQ files can be organized in the format commonly used by CDC core genomic sequencing lab, in which FastQ files are nested one level under a subfolder.  The subfolder will become the sample name. Here is an example:

```
/[path-to-directory]/NovaSeq/200710_A01000_0210_GF4RBVWJY8_MLA16554/   
│
└───MLA16554A10_2022_660_GA
│   │   MLA16554A10_2022_660_GA_S10_R1_001.fq.gz
│   │   MLA16554A10_2022_660_GA_S10_R2_001.fq.gz
│
└───MLA16554A11_2022_662_GA
│   │   MLA16554A11_2022_662_GA_S10_R1_001.fq.gz
│   │   MLA16554A11_2022_662_GA_S10_R1_001.fq.gz
```

> If your file structure looks as above, you would specify the input directory like this:
```console
--indir /[path-to-directory]/NovaSeq/200710_A01000_0210_GF4RBVWJY8_MLA16554/ --file_levels nested
```

> And the pipeline would automatically create samplesheet like this:
```console
sample,fastq_1,fastq_2
MLA16554A10_2022_660_GA,/[path-to-directory]/NovaSeq/200710_A01000_0210_GF4RBVWJY8_MLA16554/MLA16554A10_2022_660_GA_S10_R1_001.fq.gz,/[path-to-directory]/NovaSeq/200710_A01000_0210_GF4RBVWJY8_MLA16554/MLA16554A10_2022_660_GA_S10_R2_001.fq.gz
MLA16554A11_2022_662_GA,/[path-to-directory]/NovaSeq/200710_A01000_0210_GF4RBVWJY8_MLA16554/MLA16554A11_2022_662_GA_S10_R1_001.fq.gz,/[path-to-directory]/NovaSeq/200710_A01000_0210_GF4RBVWJY8_MLA16554/MLA16554A11_2022_662_GA_S10_R2_001.fq.gz
```  

> Option 2: Your FastQ files can be directly under the input folder you pass.  The sample name will be the filename without the `R1_001` or `R2_001` portion.  For single-end FastQ file, the pipeline expects the `R1_001` suffix. Here is an example:

```
/[path-to-directory]/mpox_samples_2023/
│   MLA16554A10_2022_660_GA_S10_R1_001.fq.gz
│   MLA16554A10_2022_660_GA_S10_R2_001.fq.gz
│   MLA16554A11_2022_662_GA_S10_R1_001.fq.gz
│   MLA16554A11_2022_662_GA_S10_R2_001.fq.gz
```

> If your file structure looks as above, you would specify the input directory like this:
```console
--indir /[path-to-directory]/mpox_samples_2023/
```

> The pipeline would automatically create samplesheet like this:
```console
sample,fastq_1,fastq_2
MLA16554A10_2022_660_GA_S10,/[path-to-directory]/mpox_samples_2023/MLA16554A10_2022_660_GA_S10_R1_001.fq.gz,/[path-to-directory]/mpox_samples_2023/MLA16554A10_2022_660_GA_S10_R2_001.fq.gz
MLA16554A11_2022_662_GA_S10,/[path-to-directory]/mpox_samples_2023/MLA16554A11_2022_662_GA_S10_R1_001.fq.gz,/[path-to-directory]/mpox_samples_2023/MLA16554A11_2022_662_GA_S10_R2_001.fq.gz
```  

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once (e.g. to increase sequencing depth). The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

If you provide a samplesheet, the pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below. If you have single-end data, make sure you still include the first 3 columns, and then only list the forward read followed by a comma. See an example at [samplesheet.test-single.csv](../assets/samplesheet.test-single.csv).

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.test.csv) has been provided with the pipeline. 

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once (e.g. to increase sequencing depth). The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Single-end vs. Paired-end reads

The pipeline assumes paired-end reads by default.  You should specify  `--paired False` for single-end reads when you run the pipeline with a directory path as input (`--indir`). 
If you're providing a samplesheet, the pipeline will automatically detect paired-end vs. single-end reads based on the samplesheet, as described above.

## Entrypoint/Subworkflow options
This workflow contains several subworkflows that allow you as the user to specify which components of the workflow you would like to run. You must select one subworkflow using the `--workflow` parameter. Options include `full`, `filter_reads`, `ref_based`, and `denovo`. Read filtering is always run, but can optionally be run alone. 

## Kraken2 DB options
The default kraken2 database was assembled from 15 Orthopox viruses including Monkeypox virus Clade II (MA001_USA_002), Monkeypox virus Clade I (NC_003310), Variola major (NC_001611), Variola minor (DQ441419), Borealpox (MN240300), Camelpox (NC_003391), Cowpox (NC_003663), Ectromelia (NC_004105), Vaccinia (NC_006998), Taterapox (NC_008291), Racoonpox (NC027213), Volepox (NC_031033), Skunkpox (NC_031038), Akhmeta (NC_055230), Horsepox (NC_066642), and the human t2t reference (GCF_009914755). We also include the path to smaller database with only Monkeypox virus Clade II (MA001_USA_002), and human genome build GRCh38 (GCF_000001405.40) which works quite well on Monkeypox virus Clade II samples. 

By default the pipeline will keep reads pertaining to any of these NCBI taxon IDs which are defined in the file `assets/kraken2_tax_ids.txt`. You can modify this file and point to it for filtering with the `--kraken2_tax_ids` parameter.

You can also recreate this database with the following commands:

Install mamba
```
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh
bash Mambaforge-$(uname)-$(uname -m).sh -b -p $HOME/mambaforge
```
Add mamba to your PATH
```
export PATH=$PATH:/your/mamba/directory/path/bin
```
Install Kraken2
```
mamba install -c bioconda kraken2 -y
```
Copy all your fasta files to a directory, in this case, Fasta_dir
Now use the `add to library`` command to add your fasta files to the kraken database
```
for fasta in Fasta_dir/* ; 
  do kraken2-build --add-to-library $fasta --db orthopox_kdb --threads $THREADS; 
done
```
Download the reference NCBI taxonomy files
```
kraken2-build --download-taxonomy --db orthopox_kdb --threads $THREADS
```
Now build the database
```
kraken2-build --build --db orthopox_kdb --threads $THREADS
```

## Running the pipeline

The typical command for running the pipeline with a sample sheet is as follows:

```console
nextflow run polkapox --input {SAMPLESHEET.csv} --outdir {OUTDIR} --genome {REF.fa} -profile <docker, singularity, test etc> --kraken_db {PATH/TO/DB} --workflow {WORKFLOW} --filter {true/false}
```

To run with an input directory, run as: 

```console
nextflow run polkapox --indir {PATH/TO/DIR} --outdir {OUTDIR} --genome {REF.fa} -profile <docker, singularity, test etc> --kraken_db {PATH/TO/DB} --workflow {WORKFLOW} --filter {true/false}
```

To run with SRA accessions: 

```console
nextflow run polkapox --sra --sra_ids {PATH/TO/SRA/FILE} --outdir {OUTDIR} --genome {REF.fa} -profile <docker, singularity, test etc> --kraken_db {PATH/TO/DB} --workflow {WORKFLOW} --filter {true/false}
```

Note that the pipeline will create the following files in your working directory:

```console
work                # Directory containing the nextflow working files
<OUTIDR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull polkapox
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [polkapox releases page](https://github.com/CDCgov/polkapox/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
  - This profile will run subsampled data which will publish contigs, if you want to test a full assembly run test_full

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/software/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
