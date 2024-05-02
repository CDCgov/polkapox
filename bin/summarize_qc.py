#!/usr/bin/env python3

import pandas as pd
import os, sys
from glob import glob
import json
import argparse
import subprocess
import logging

logger = logging.getLogger()

def parse_args():

    parser = argparse.ArgumentParser(
        description="Aggregate QC from multiQC directory.",
        epilog="Example: python summarize_qc.py --analysis_dir multiqc_dir/",
    )
    parser.add_argument(
        "--analysis_dir",
        metavar="ANALYSIS_DIR",
        help="Analysis directory")
    parser.add_argument(
        "--samplesheet",
        metavar="SAMPLE_SHEET",
        help="The path to the sample sheet"),
    parser.add_argument(
        "--reference_genome",
        metavar="REFERENCE_GENOME",
        required=False,
        help="The name of the reference genome used for reference-based assembly"),
    parser.add_argument(
        "--kraken_db",
        metavar="KRAKENDB",
        help="The name of the Kraken database"),
    parser.add_argument(
        "--kraken_tax_ids",
        metavar="KRAKENTAXIDS",
        help="The name of the tax IDs used for Seqtk filtering"),
    parser.add_argument(
        "--filter",
        metavar="FILTER",
        help="true or false to filter variants to only params.coord coordinate boundaries"),
    parser.add_argument(
        "--workflow",
        metavar="WORKFLOW",
        help="Nextflow entrypoint, options are filter_reads, ref_based, denovo, or full"),
    parser.add_argument(
        "--coords",
        metavar="COORDS",
        required=False,
        help="Four number list with start,stop for each of the two loci of interest used for counting variants of interest"),
    parser.add_argument(
        "--locus1",
        metavar="LOCUS1",
        required=False,
        help="Feature name of the first coordinates used for filtering"),
    parser.add_argument(
        "--locus2",
        metavar="LOCUS2",
        required=False,
        help="Feature name of the second coordinates used for filtering"),
    
    return parser.parse_args()

def get_raw_filt_counts(search_dir, sample):
    """ Get raw filter counts from fastp output
    :param search_dir: work directory
    :param sample: sample name
    :returns: total raw reads and total filtered reads
    :rtype: tuple
    """
    p = "**/*{}.fastp.json".format(sample)
    fastp_file = glob(os.path.join(search_dir, p), recursive=True)[0]

    f = open(fastp_file)
    try:
        data = json.load(f)
        logger.info(f"Data loaded from {fastp_file}")
    except:
        logger.error(f"Could not load data from {fastp_file}")
        sys.exit(1)

    total_raw_reads = data['summary']['before_filtering']['total_reads']
    total_filtered_reads = data['summary']['after_filtering']['total_reads']

    return total_raw_reads, total_filtered_reads

def get_kraken_stats(search_dir, sample, kraken_db, kraken_tax_ids):
    """ Get stats from kraken output
    :param search_dir: work directory
    :param sample: sample name
    :param kraken_db: path to kraken database
    :kraken_tax_ds: path to kraken tax ids
    :returns: kraken statistics
    :rtype: tuple
    """

    p = "**/*{}.kraken2.classifiedreads.txt".format(sample)
    kraken_reads = glob(os.path.join(search_dir, p), recursive=True)[0]
    try:
        k_data = pd.read_csv(kraken_reads, delim_whitespace=True, usecols=[0,1,2,3,4], header=None)
        logger.info(f"Data load from {kraken_reads}")
    except:
        logger.error(f"Unable to load data from {kraken_reads}")
        sys.exit(1)
    
    total = len(k_data)

    s = "**/{}*.opxreads.txt".format(sample)
    ortho_reads = glob(os.path.join(search_dir, s), recursive=True)[0]

    if os.path.getsize(ortho_reads) > 0:
        s_data = pd.read_csv(ortho_reads, delim_whitespace=True, header=None)
        opx = len(s_data)
    else:
        opx = 0

    human = (k_data[2] == 9606).sum()
    unclass = (k_data[2] == 0).sum()
    
    with open(kraken_tax_ids) as l:
        lines = [line.strip() for line in l.readlines()]
    k_tax_ids = ', '.join(lines)

    if opx == 0:
        opx_perc = 0
    else:
        opx_perc = round(((opx/total) * 100),2)
    human_perc = round(((human/total) * 100),2)
    unclass_perc = round(((unclass/total) * 100),2)

    return total, opx_perc, human_perc, unclass_perc, kraken_db, k_tax_ids

def get_flagstat_denovo(search_dir, sample):
    """ Get samtools flagstat results from denovo mapping
    :param search_dir: work directory
    :param sample: sample name
    :returns: denovo mapping statistics
    :rtype: tuple
    """
    #if path not exist output unknown
    stats=[]
    p = "{}.denovo.flagstat".format(sample)
    if os.path.exists(p):
        logger.info(f"Path {p} found")
        fh=open(p,'r')
        for line in fh.readlines()[0:2]:
            reads=line.split(' ')[0]
            stats.append(reads)
        total_reads_denovo=round(int(stats[0]),2)
        mapped_reads_denovo=round(int(stats[1]),2)
        percent_mapped_denovo=round((mapped_reads_denovo/total_reads_denovo)*100,2)
    else:
        logger.info(f"Path {p} does not exist")
        total_reads_denovo='NaN'
        mapped_reads_denovo='NaN'
        percent_mapped_denovo='NaN'
    
    return total_reads_denovo, mapped_reads_denovo, percent_mapped_denovo

def get_cov_stats(search_dir, sample, reference):
    """ Get coverage stats from samtools output
    :param search_dir: work directory
    :param sample: sample name
    :param reference: path to reference genome
    :returns: coverage statistics
    :rtype: tuple
    """
    p = "**/*{}.depth.tsv".format(sample)
    depth_file = glob(os.path.join(search_dir, p), recursive=True)[0]
    
    if os.path.getsize(depth_file) > 0:
        logger.info(f"Retrieving depth statistics from {depth_file}")
        dp = pd.read_csv(depth_file, sep="\t", header = None)
        avg_dp = int(round(dp[2].mean(),2))
        dp_gt_20 = int(round(len(dp[dp[2] >= 20]),2))
    else:
        logger.info(f"{depth_file} has size zero")
        dp = 0
        avg_dp = 0
        dp_gt_20 = 0
    ref_genome = reference.split('/')[-1]

    return avg_dp, dp_gt_20, ref_genome

def get_gfa_stats(search_dir, sample):
    """ Get stats from Unicycler graph assembly output
    :param search_dir: work directory
    :param sample: sample name
    :returns:gfaResults
    :rtype: list
    """
    p = "**/*{}.assembly.log".format(sample)
    #print(p)
    gfa_log = glob(os.path.join(search_dir, p), recursive=True)
    if gfa_log:
       gfa_log = glob(os.path.join(search_dir, p), recursive=True)[0]
        
       f = open(gfa_log)
       parsed_json = json.load(f)    

       if len(parsed_json.keys()) == 11:
           notes = 'GFA step complete'
           successCode_step9 = list(parsed_json)[-2]
           successCode_step10 = list(parsed_json)[-1]
           final_order_orientation_copy_number = parsed_json[successCode_step9]['output']['final_order_orientation_copy_number']
           final_sequence_length = parsed_json[successCode_step9]['output']['final_sequence_length']
           status = parsed_json[successCode_step10]['status']
           final_itr_length = parsed_json[successCode_step10]['output']['final_itr_length']
           gfaResults = [final_order_orientation_copy_number, float(final_sequence_length), float(final_itr_length), status, notes]
       else:
           failCode_stepN = list(parsed_json)[-1]
           status = parsed_json[failCode_stepN]['status']
           statusReport = 'FAIL'
           gfaResults = ['Unknown','Unknown','Unknown', statusReport, status]
    else:
        status = 'Unicycler-GFA Log No Exist'
        statusReport = 'FAIL'
        gfaResults = ['Unknown','Unknown','Unknown', statusReport, status]
    return gfaResults

def fix_names(df):
    """ Standardize the raw column names
    :param df: input dataframe with columns to standardize
    :returns: a dataframe with standardized columns
    :rtype: Dataframe object
    """
    raw_col_names = ['Sample',
                     'Samtools_mqc-generalstats-samtools-flagstat_total',
                     'Samtools_mqc-generalstats-samtools-mapped_passed',
                     'QUAST_mqc-generalstats-quast-N50',
                     'QUAST_mqc-generalstats-quast-Total_length',
                     'Kraken_mqc-generalstats-kraken-Homo_sapiens',
                     'Kraken_mqc-generalstats-kraken-Monkeypox_virus',
                     'Kraken_mqc-generalstats-kraken-Top_5',
                     'Kraken_mqc-generalstats-kraken-Unclassified',
                     'fastp_mqc-generalstats-fastp-pct_duplication',
                     'fastp_mqc-generalstats-fastp-after_filtering_q30_rate',
                     'fastp_mqc-generalstats-fastp-after_filtering_q30_bases',
                     'fastp_mqc-generalstats-fastp-after_filtering_gc_content',
                     'fastp_mqc-generalstats-fastp-pct_surviving',
                     'fastp_mqc-generalstats-fastp-pct_adapter']
    
    for col_name in raw_col_names[:]:
        # check for human, mpox virus top hit column. add any missing columns from raw_col_names list
        if col_name not in df.columns:
            # if the missing column is either human, mpox virus top hit, remove it from the raw_col_names list
            if col_name == 'Kraken_mqc-generalstats-kraken-Homo_sapiens' or col_name == 'Kraken_mqc-generalstats-kraken-Monkeypox_virus':
                raw_col_names.remove(col_name)
            else:
                # add any other missing columns and give a value of NA
                df[col_name] = "NA"
    
    df = df[raw_col_names]

    if len(raw_col_names) == 14:
        for i in range(0,len(raw_col_names)):
            if not df.columns[i] == raw_col_names[i]:
                logger.error(f"Unexpected column names in multiqc/*general_stats.txt file")
                sys.exit(1)
    else:
        logger.error("Unexpected column names in multiqc/*general_stats.txt file")
        sys.exit(1)
    
    return df

def get_total_snps(search_dir, sample):
    """ Grab the tsv files from ivar output (ivar_variants dir) and get total number of SNPs
    :param search_dir: work directory
    :param sample: sample name
    :returns: total number of snps
    :rtype: int
    """
    p = "**/*{}.ivar.tsv".format(sample)
    snpFile = glob(os.path.join(search_dir, p), recursive=True)[0]
    try:
        df = pd.read_csv(snpFile, sep="\t", header = 0)
        logger.info(f"Successfully read {snpFile}")
    except:
        logger.error(f"Unable to read {snpFile}")
    total_snps = len(df)

    return total_snps

def get_snp_metadata(search_dir, sample, coords):
    """ Get number of SNPs for input coordinates from ivar_summary file (in variant_summaries dir)
    :param search_dir: work directory
    :param sample: sample name
    :param coords: a set of 2 coordinates
    :returns: total number of snps for each of 2 coordinates
    :rtype: tuple
    """    
    C1_count = 0
    C2_count = 0
    coord2_start = coords.split(',')[2]
    p = "**/*{}_ivar_summary.txt".format(sample)
    try:
        snpMeta = glob(os.path.join(search_dir, p), recursive=True)[0]
        logger.info(f"Found {snpMeta}")
    except:
        logger.error(f"Unable to open {snpMeta}")
        sys.exit(1)
    for line in open(snpMeta,'r').readlines():
        if line.split('\t')[1] < coord2_start:
            C1_count += 1
        else:
            C2_count += 1

    return C1_count, C2_count

def get_polish_stats(sample):
    """ Parse the mummer report file and return # corrected SNPs and indels
    :param sample: sample name
    :returns: total number of snps and indels
    :rtype: tuple
    """
    SNPs = None
    Indels = None
    infile = sample + '.report'

    if os.path.exists(infile):
        logger.info(f"{infile} found")
        with open(infile, 'r') as file:
            for line in file:
                if line.startswith("TotalSNPs"):
                    parts = line.split()
                    if len(parts) >= 2:
                        SNPs = round(int(parts[1]),2)
                elif line.startswith("TotalIndels"):
                    parts = line.split()
                    if len(parts) >= 2:
                        Indels = round(int(parts[1]),2)
                    break  # Stop reading after finding TotalIndels
    else:
        logger.info(f"{infile} does not exist")
        SNPs = 'NaN'
        Indels = 'NaN'
    return SNPs, Indels

def count_ns_in_pileup(sample):
    """ Get the consensus SNPs from final mpileup file 
    :param sample: sample name
    :returns: total number of snps 
    :rtype: int
    """
    p = "{}.final.mpileup".format(sample)
    try:
        command = f"cat {p} | awk '($4 < 20){{count++}} END {{print count}}'"
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logger.info(f"Successfully parsed {p}")
    except:
        logger.error(f"Unable to parse {p}")
    if result.returncode == 0:
        return str(result.stdout.strip())
    else:
            return None
    
def main():
    args = parse_args()

    summary_file = glob(os.path.join(args.analysis_dir, '**/*general_stats.txt'), recursive=True)[0]
    sample_file = args.samplesheet
    logger.info(f"Running summary_qc.py for {args.workflow}")
    try:
        summary = pd.read_csv(summary_file, delimiter = "\t")
    except:
        logger.error(f"Unable to open {summary_file}")
        sys.exit(1)
    try:
        samples = pd.read_csv(sample_file, delimiter = ",")
    except:
        logger.error(f"Unable to open {sample_file}")
        sys.exit(1)

    ### Check number of entries per dataframe
    summaryIds = list(summary[summary.columns[0]])
    samplesIds = list(samples[samples.columns[0]])

    sample_summary_notProcessed = set(samplesIds) - set(summaryIds)
    
    if not len(list(sample_summary_notProcessed)) == 0:
         rowMissing = len(summary.columns)-1
         for i in sample_summary_notProcessed:
               naHold = ['NaN'] * rowMissing
               no_sample_data = [i] + naHold
               summary.loc[len(summary)] = no_sample_data
               summary.fillna(value='NaN',inplace=True)

    new_col_names = ['sample', 'reads_total_bwa', 'reads_mapped_bwa', 'n50_unicycler','assembly_length_unicycler','top_taxa_percent_kraken','top5_taxa_percent_kraken','unclassified_percent_kraken','percent_duplication_fastp','q30_rate_postfilter_fastp','q30_bases_postfilter_fastp','gc_content_postfilter_fastp','percent_reads_passed_fastp','percent_adapter_fastp']
    fixed_summary = fix_names(summary)
    fixed_summary.set_axis(new_col_names, axis=1, inplace=True)    

    if args.workflow == 'denovo' or args.workflow == 'full':
        contig_files = glob(os.path.join(args.analysis_dir, '**/*quast_num_contigs_1.txt'), recursive=True)
        if contig_files:
            contig_file = contig_files[0]
            logger.info(f"{contig_file} found")
            contigs = pd.read_csv(contig_file, delimiter = "\t")
            ### Check number of entries per dataframe
            contigsIds = list(contigs[contigs.columns[0]])
            sample_contig_notProcessed = set(samplesIds) - set(contigsIds)
            
            if not len(list(sample_contig_notProcessed)) == 0:
                rowMissing = len(contigs.columns)-1
                for i in sample_contig_notProcessed:
                    naHold = ['NaN'] * rowMissing
                    no_sample_data = [i] + naHold
                    contigs.loc[len(contigs)] = no_sample_data
                    contigs.fillna(value='NaN',inplace=True)

            contigs['n_contigs_unicycler'] = contigs.sum(numeric_only=True,axis = 1)
            fixed_summary = fixed_summary.merge(contigs[['Sample', 'n_contigs_unicycler']],left_on='sample', right_on='Sample').drop('Sample', axis = 1)
        else:
            # If Unicycler run failed, config_files don't get created
            logger.info(f"{contig_file} not found")
            fixed_summary['n_contigs_unicycler'] = 'NaN'    

    if args.workflow == 'ref_based' or args.workflow == 'full':
        fixed_summary['percent_mapped_bwa'] = round(((fixed_summary['reads_mapped_bwa'] / fixed_summary['reads_total_bwa']) * 100),2)
    
    summary_full = fixed_summary.reindex(fixed_summary.columns.tolist() + ['raw_read_count_fastp', 'filtered_read_count_fastp', 'total_raw_reads', 'opx_percent_kraken','human_percent_kraken', 'unclass_percent_kraken', 'average_depth_bwa', 'count_20xdepth_bwa'], axis=1)
    summary_full.to_csv("before_adding_values.csv")

    # get data from other non multiqc input files:
    for sample in summary_full['sample']:
        summary_full.loc[summary_full['sample'] == sample, ['opx_read_count_kraken', 'filtered_read_count_fastp']] = get_raw_filt_counts(args.analysis_dir, sample)
        summary_full.loc[summary_full['sample'] == sample, ['total_raw_reads', 'opx_percent_kraken', 'human_percent_kraken', 'unclass_percent_kraken', 'kraken_db','kraken_tax_ids']] = get_kraken_stats(args.analysis_dir, sample, args.kraken_db, args.kraken_tax_ids)
        if args.workflow == 'ref_based':
            summary_full.loc[summary_full['sample'] == sample, ['average_depth_bwa', 'count_20xdepth_bwa','reference_genome']] = get_cov_stats(args.analysis_dir, sample, args.reference_genome)
            summary_full.loc[summary_full['sample'] == sample, ['total_snps']] = get_total_snps(args.analysis_dir,sample)    
            summary_full.loc[summary_full['sample'] == sample, ['corrected_Ns']] = count_ns_in_pileup(sample)    
            #logic to handle if a user is not interested in variant filtering
            if args.filter == 'true':
                summary_full.loc[summary_full['sample'] == sample, [f'{args.locus1}_SNPs',f'{args.locus2}_SNPs']] = get_snp_metadata(args.analysis_dir, sample, args.coords)
        elif args.workflow == 'denovo':
            summary_full.loc[summary_full['sample'] == sample, ['orientation_copy_number' ,'sequence_length', 'itr_length', 'gfa_status', 'gfa_notes']] = get_gfa_stats(args.analysis_dir, sample)
            summary_full.loc[summary_full['sample'] == sample, ['total_reads_denovo','mapped_reads_denovo','percent_mapped_denovo']] = get_flagstat_denovo(args.analysis_dir, sample)
            summary_full.loc[summary_full['sample'] == sample, ['corrected_snps','corrected_indels']] = get_polish_stats(sample)
            summary_full.loc[summary_full['sample'] == sample, ['corrected_Ns']] = count_ns_in_pileup(sample)    
        elif args.workflow == 'full':
            summary_full.loc[summary_full['sample'] == sample, ['average_depth_bwa', 'count_20xdepth_bwa', 'reference_genome']] = get_cov_stats(args.analysis_dir, sample, args.reference_genome)
            summary_full.loc[summary_full['sample'] == sample, ['total_snps']] = get_total_snps(args.analysis_dir,sample)    
            summary_full.loc[summary_full['sample'] == sample, ['corrected_snps','corrected_indels']] = get_polish_stats(sample)    
            summary_full.loc[summary_full['sample'] == sample, ['corrected_Ns']] = count_ns_in_pileup(sample)    
            summary_full.loc[summary_full['sample'] == sample, ['orientation_copy_number' ,'sequence_length', 'itr_length', 'gfa_status', 'gfa_notes']] = get_gfa_stats(args.analysis_dir, sample)
            summary_full.loc[summary_full['sample'] == sample, ['total_reads_denovo','mapped_reads_denovo','percent_mapped_denovo']] = get_flagstat_denovo(args.analysis_dir, sample)
            if args.filter == 'true':
                summary_full.loc[summary_full['sample'] == sample, [f'{args.locus1}_SNPs',f'{args.locus2}_SNPs']] = get_snp_metadata(args.analysis_dir, sample, args.coords)

    # final column order
    if args.workflow == 'filter_reads':
        summary_full = summary_full[['sample','total_raw_reads','opx_read_count_kraken','opx_percent_kraken','human_percent_kraken','unclass_percent_kraken','kraken_db','kraken_tax_ids','filtered_read_count_fastp','percent_reads_passed_fastp','percent_adapter_fastp','gc_content_postfilter_fastp','q30_rate_postfilter_fastp','percent_duplication_fastp']]
    elif args.workflow == 'ref_based':
        if args.filter == 'true':
            summary_full = summary_full[['sample','reference_genome','total_raw_reads','opx_read_count_kraken','opx_percent_kraken','human_percent_kraken','unclass_percent_kraken','kraken_db','kraken_tax_ids','filtered_read_count_fastp','percent_reads_passed_fastp','percent_adapter_fastp','gc_content_postfilter_fastp','q30_rate_postfilter_fastp','percent_duplication_fastp','reads_mapped_bwa','percent_mapped_bwa','average_depth_bwa','count_20xdepth_bwa','total_snps',f'{args.locus1}_SNPs',f'{args.locus2}_SNPs']]
        if args.filter == 'false':
            summary_full = summary_full[['sample','reference_genome','total_raw_reads','opx_read_count_kraken','opx_percent_kraken','human_percent_kraken','unclass_percent_kraken','kraken_db','kraken_tax_ids','filtered_read_count_fastp','percent_reads_passed_fastp','percent_adapter_fastp','gc_content_postfilter_fastp','q30_rate_postfilter_fastp','percent_duplication_fastp','reads_mapped_bwa','percent_mapped_bwa','average_depth_bwa','count_20xdepth_bwa','total_snps']]
    elif args.workflow == 'denovo':
            summary_full = summary_full[['sample','total_raw_reads','opx_read_count_kraken','opx_percent_kraken','human_percent_kraken','unclass_percent_kraken','kraken_db','kraken_tax_ids','filtered_read_count_fastp','percent_reads_passed_fastp','percent_adapter_fastp','gc_content_postfilter_fastp','q30_rate_postfilter_fastp','percent_duplication_fastp','n_contigs_unicycler','assembly_length_unicycler','n50_unicycler','mapped_reads_denovo','percent_mapped_denovo','orientation_copy_number','sequence_length','itr_length','gfa_status','gfa_notes','corrected_snps','corrected_indels','corrected_Ns']]
    elif args.workflow == 'full':
        if args.filter == 'true':
            summary_full = summary_full[['sample','reference_genome','total_raw_reads','opx_read_count_kraken','opx_percent_kraken','human_percent_kraken','unclass_percent_kraken','kraken_db','kraken_tax_ids','filtered_read_count_fastp','percent_reads_passed_fastp','percent_adapter_fastp','gc_content_postfilter_fastp','q30_rate_postfilter_fastp','percent_duplication_fastp','reads_mapped_bwa','percent_mapped_bwa','average_depth_bwa','count_20xdepth_bwa','n_contigs_unicycler','assembly_length_unicycler','n50_unicycler','mapped_reads_denovo','percent_mapped_denovo','orientation_copy_number','sequence_length','itr_length','gfa_status','gfa_notes','total_snps',f'{args.locus1}_SNPs',f'{args.locus2}_SNPs','corrected_snps','corrected_indels','corrected_Ns']]
        if args.filter == 'false':
            summary_full = summary_full[['sample','reference_genome','total_raw_reads','opx_read_count_kraken','opx_percent_kraken','human_percent_kraken','unclass_percent_kraken','kraken_db','kraken_tax_ids','filtered_read_count_fastp','percent_reads_passed_fastp','percent_adapter_fastp','gc_content_postfilter_fastp','q30_rate_postfilter_fastp','percent_duplication_fastp','reads_mapped_bwa','percent_mapped_bwa','average_depth_bwa','count_20xdepth_bwa','n_contigs_unicycler','assembly_length_unicycler','n50_unicycler','mapped_reads_denovo','percent_mapped_denovo','orientation_copy_number','sequence_length','itr_length','gfa_status','gfa_notes','total_snps','corrected_snps','corrected_indels','corrected_Ns']]

    summary_full.to_csv("sample_summary.tsv", sep = "\t", index = False)
    logger.info(f"Summary results successfully written to sample_summary.tsv")

if __name__ == "__main__":
    main()
