#!/usr/bin/env python3

import os,sys
import argparse
import logging
import re
from glob import glob
from pathlib import Path

# Script to generate sample sheet from a directory of fastq files
# Paired-end files must end in (1 or R1) and (2 or R2)
# Sample sheet format is sample,fastq_1,fastq_2 (paired-end) or sample,fastq (single-end)
logger = logging.getLogger()

def parse_args():
    """ Set up arguments
    :returns: an argparse.Namespace object with arguments 
    :rtype: argparse.Namespace
    """

    parser = argparse.ArgumentParser(
        description="Generate a sample sheet from a directory of fastq files.",
        epilog="Example: python3 create_samplesheet.py -d ./my_samples_folder/",
    )
    parser.add_argument(
        '-d', '--indir',
	default=None,
	required=False,
        metavar="SAMPLES_DIR",
        help="Samples directory")
    parser.add_argument(
        '--single',
        action='store_true',
        required=False,
        help="Flag for single-end reads")
    parser.add_argument(
        '-o', '--outdir',
        default=None,
	    required=False,
        metavar="OUTPUT_DIR",
        help="Path to project dir where samplesheet will be saved")
    parser.add_argument(
        '--project_name',
        default=None,
	    required=False,
        metavar="PROJECT_NAME",
        help="Project name which will be the prefix of samplesheet name")
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args()

def check_path(input_path):
    """ Check if path is valid and if yes, return absolute path
        Add a trailing slash '/' if there isn't one
    :param path: a file path
    :returns: a validated, absolute file path
    :rytpe: str
    """
    if not os.path.isdir(input_path):
        logger.error(f"The given input {input_path} is not a valid directory")
        sys.exit(1)
    else:
        input_path = os.path.realpath(input_path) + '/'
    return input_path

def check_if_R1(file_path):
    """ Return true if file has R1 in the name
    :param path: a file path
    :returns: true if valid filenames, else false
    :rytpe: bool
    """
    return "_R1" in file_path

def find_file_depth(samples_dir):
    """ Return the file depth
    :param samples_dir: a file path
    :returns: dictionary of fastq file levels under samples_dir
    :rtype: dict
    """
    extensions = (".fastq",".fq",".fastq.gz",".fq.gz")
    file_depth = {}
    for dirpath, _, files in os.walk(samples_dir):
        relative_path = os.path.relpath(dirpath, samples_dir)
        depth = 0 if relative_path == "." else relative_path.count(os.path.sep) + 1
        # Ignore files nested >1 level down because these are not handled by rest of code
        if depth > 1:
            continue
        for file in files:
            if file not in file_depth and file.endswith(extensions):
                file_depth[file] = depth   
    return file_depth

def handle_file_depth(depth_dict):
    """ Return whether files are directly under sampls_dir or nested 1 level down
    :param depth_dict: a dictionary with filepaths as keys and file level as value
    :returns: 0 if directly under samples_dir, 1 if nested one level down
    :rtype: boolean
    """
    unique_depths = set(depth_dict.values())
    if unique_depths == {0}:
        nested = 0 # directly under input dir
    elif unique_depths == {1}:
        nested = 1 # nested one level down
    else:
        # Store err val if files are at mixed depths or unexpected val in dict
        nested = 9
    return nested

def remove_id(sample_name):
    """ Removes the sample identifier common in CDC Core Sequencing Lab samples
        Identifier format: a leading underscore, following by 10 numeric digits, at the end of the sample name string
    :param sample_name: the sample name as a string
    :returns: the sample name without the identifier
    :rtype: str
    """
    id_pattern = re.compile(r'_\d{10}$')
    new_sample_name = re.sub(id_pattern, '', sample_name)
    return new_sample_name

def list_samples(samples_dir, single=False):
    """ Return a list of all samples in a directory
    :param samples_dir: a file path
    :returns: a list of paths to all fastq/fq files in samples_dir
    :rtype: list
    """
    samples_dir = os.path.realpath(samples_dir)
    extensions = (".fastq",".fq",".fastq.gz",".fq.gz")
    seqfiles = [] 
    for filename in os.listdir(samples_dir):
        if single and "_R2" in filename:
            log.error("single flag is set to {single} but input directory contains R2 files")
        # Check for fastq files directly under samples_dir
        if filename.endswith(extensions) and "_R1_" in filename:
            seqfiles.append(os.path.join(samples_dir, filename))       
        else:
            # Check for fastq files nested one level down
            subdir = os.path.join(samples_dir, filename)
            if os.path.isdir(subdir):
                for subfilename in os.listdir(subdir):
                    if subfilename.endswith(extensions) and "_R1_" in subfilename:
                        seqfiles.append(os.path.join(subdir, subfilename))
    return seqfiles

def create_samplesheet(samples_list, outdir, outfile, nested, single=False):
    """ Create a samplesheet from files in a directory
    :param samples_dir: a directory of fastq files
    :param samples_list: a list of all samples in samples_dir 
    :returns:
    :rtype:
    """
    #todo: Add a date to the samplesheet name?    
    if single:
        # process single-end files
        with open(f'{outdir}{outfile}', 'w') as samplesheet:
            samplesheet.write('sample,fastq_1,fastq_2\n')
            for sample in sorted(samples_list):
                if nested == 1:
                    s_name = Path(sample).parent.name
                    s_name = remove_id(s_name)
                else:
                    s_name = os.path.splitext(Path(sample).stem)[0].replace('_R1_001','')
                samplesheet.write(f'{s_name},{sample},\n')
    else:
        # process as paired-end files
        with open(f'{outdir}{outfile}', 'w') as samplesheet:
            samplesheet.write('sample,fastq_1,fastq_2\n')
            for sample in sorted(samples_list):
                if nested == 1:
                    s_name = s_name = Path(sample).parent.name
                    s_name = remove_id(s_name)
                else:
                    s_name = os.path.splitext(Path(sample).stem)[0].replace('_R1_001','')
                r2 = sample.replace('R1','R2')
                samplesheet.write(f'{s_name},{sample},{r2}\n')

def main():
    args = parse_args()
    input_dir = check_path(args.indir)
    output_dir = check_path(args.outdir)
    if args.project_name:
        outfile = f"{args.project_name}_samplesheet.csv"
    else:
        outfile = f"samplesheet.csv"
    if args.single:
        samples_list = list_samples(input_dir, single=True)
    else:
        samples_list = list_samples(input_dir)

    # Find whether fastq files are nested 1 level down or directly under input dir
    depths_dict = find_file_depth(input_dir)
    nested = handle_file_depth(depths_dict)
    # Error if files are nested in a way the script doesn't handle
    if nested == 9:
        logger.error(f"fastq files in {input_dir} must be either directly under input folder or nested one level down")
        sys.exit(1)
    if args.single:
        create_samplesheet(samples_list,output_dir, outfile, nested, single=True)
    else:
        create_samplesheet(samples_list, output_dir, outfile, nested)
    logger.info(f"{outfile} created. Please review and be sure it is correct.")

if __name__ == '__main__':
    main()
