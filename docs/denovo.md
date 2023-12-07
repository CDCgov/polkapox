# PolkaPox: De novo Assembly

## Introduction

This document describes the de novo assembly steps of the pipeline. Additional details regarding the output files are provided in [output](/docs/output.md).

## Overview

De novo genome assembly of filtered and trimmed sequencing reads consists of the following steps:  
1. Unicycler draft assembly
1. Assembly path reconstruction
1. Sequence polishing

## Assembly path reconstruction

A unique component of the de novo assembly pipeline within PolkaPox is the reconstruction of a linear orthopox genome from linkages between draft contigs assembled from Unicycler. Briefly, the assembly graph is read with [`gfapy`](https://gfapy.readthedocs.io/en/latest/), a library for parsing and manipulating GFA formatted objects, and then traversed in search of the longest paths starting from the largest contig. Because orthopox genomes are flanked by inverted terminal repeats (ITRs), which become collapsed during de novo assembly, the longest paths in both directions should end at the same contigs. These paths, and the relative coverage depth of contigs containing repeatitive sequences, are used to infer the complete, linear genome sequence. See [`mpxv-AssemblyGraph_gfaPy.py`](/bin/mpxv-AssemblyGraph_gfaPy.py) for more details.

![Assembly graph diagram](/docs/images/bandage_polkapox.png)*Basic assembly graph reconstruction process.*

Once the draft contigs are combined into a linear assembly, the sequence is further 'polished' by mapping filtered reads with BWA and computing the resulting consensus with iVar.

## Limitations

Graph reconstruction makes various assumptions and if Unicycler produces a graph different from expectations, an assembly path will not be produced. In such instances, only the draft multi-contig assembly is reported [`${outdir}/final_assembly/${sample}_draft.fa`].

| Example 1: Broken path | Example 2: Too many contig links|
|-----|-----|
|*WARNING: only 1 path between longest contig and the last contig in the longest path* |*WARNING: Multiple links connect to 2 or more contigs in assembly* |
|![Bad assembly graph 1](/docs/images/bandage_broke.png) |![Bad assembly graph 2](/docs/images/bandage_links.png) |

## Tips for manual resolution

If PolkaPox cannot reconstruct a linear genome sequence, it may still be possible to manually resolve the assembly through additional steps. For example:  
1. Visually inspect and manipulate the graph using a software tool like [Bandage](https://github.com/rrwick/Bandage).  
1. Align the draft contigs to a suitable reference genome and manually scaffold them together.
1. Try another assembly tool, using the filtered and trimmed reads [`${outdir}/fastp/${sample}.fastp.fastq.gz`] as input. 
1. Scaffold your contigs with additional long-read data (e.g. ONT) 

Once you have produced a suitable scaffold, consider mapping the invidual reads back to your assembly to confirm. For example, run PolkaPox with `--workflow refbased` and provide your new assembly FASTA as the reference.