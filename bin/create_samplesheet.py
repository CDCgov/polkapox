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
logger = logging.getLogger(__name__)

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
        '--file_levels',
        choices=['nested', 'all', 'top'],
        default='nested',
        required=False,
        metavar="FILE_LEVELS",
        help="Option for creating a sample sheet: 'nested' (default) for only nested files,\n"
             "'all' for all files in the directory, 'top' for only top-level files"
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

def are_files_nested(samples_dir):
    """ Return true if there are nested fastq files, else return false
    :param samples_dir: a file path
    :returns: true if there are fastq files nested one level under samples_dir
    :rtype: boolean
    """
    extensions = (".fastq",".fq",".fastq.gz",".fq.gz")
    for dirpath, _, files in os.walk(samples_dir):
        if dirpath == samples_dir:
            continue
        for file in files:
            if file.endswith(extensions):
                return True
    return False

def are_files_top_level(samples_dir):
    """ Return true if there are fastq files directly under samples_dir, else return false.
    :param samples_dir: a directory path
    :returns: true if there are files with specified extensions directly under samples_dir
    :rtype: boolean
    """
    extensions = (".fastq",".fq",".fastq.gz",".fq.gz")
    for file in os.listdir(samples_dir):
        if any(file.endswith(ext) for ext in extensions):
            return True
    return False

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

def list_samples(samples_dir, file_levels, single=False):
    """ Return a list of all samples in a directory
    :param samples_dir: a file path
    :returns: a list of paths to all fastq/fq files in samples_dir
    :rtype: list
    """
    samples_dir = os.path.realpath(samples_dir)
    extensions = (".fastq",".fq",".fastq.gz",".fq.gz")
    seqfiles = {} 
    for filename in os.listdir(samples_dir):
        if single and "_R2" in filename:
            logger.error("single flag is set to {single} but input directory contains R2 files")
        if file_levels == 'top' or file_levels == 'all':
            # Check for fastq files directly under samples_dir
            if filename.endswith(extensions) and ("_R1" in filename or "_1." in filename):
                #seqfiles.append(os.path.join(samples_dir, filename))
                sample_path = os.path.join(samples_dir, filename)
                s_name = os.path.splitext(Path(sample_path).stem)[0]
                s_name = re.sub(r'_R1_001$', '', s_name) # remove R1_001 at the end of filename, if it exists
                s_name = re.sub(r'_R1$', '', s_name) # remove R1 at the end of filename, if it exists
                s_name = re.sub(r'_1$', '', s_name) # remove a _1 only if it occurs at end of filename
                s_name = remove_id(s_name)
                seqfiles[s_name] = sample_path
        if file_levels == 'nested' or file_levels == 'all':
            # Check for fastq files nested one level down
            subdir = os.path.join(samples_dir, filename)
            if os.path.isdir(subdir):
                for subfilename in os.listdir(subdir):
                    if subfilename.endswith(extensions) and "_R1" in subfilename:
                        #seqfiles.append(os.path.join(subdir, subfilename))
                        sample_path = os.path.join(subdir, subfilename)
                        s_name = Path(sample_path).parent.name # use the subdirectory name as the sample name
                        s_name = remove_id(s_name)
                        seqfiles[s_name] = sample_path
    return seqfiles


def create_samplesheet(samples_dict, outdir, outfile, single=False):
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
            samplesheet.write('sample,fastq_1\n')
            for sample_name, sample_path in sorted(samples_dict.items()):
                samplesheet.write(f'{sample_name},{sample_path},\n')
    else:
        # process as paired-end files
        with open(f'{outdir}{outfile}', 'w') as samplesheet:
            samplesheet.write('sample,fastq_1,fastq_2\n')
            for sample_name, sample_path in sorted(samples_dict.items()):
                r2 = sample_path.replace('R1','R2')
                samplesheet.write(f'{sample_name},{sample_path},{r2}\n')

def main():
    args = parse_args()
    logfile = f"{os.path.basename(__file__).strip('.py')}.log"
    logging.basicConfig(filename=logfile, 
                        level=logging.DEBUG)
    input_dir = check_path(args.indir)
    output_dir = check_path(args.outdir)
    if args.project_name:
        outfile = f"{args.project_name}_samplesheet.csv"
    else:
        outfile = f"samplesheet.csv"
    
    if args.file_levels != 'all':
        if args.file_levels == 'nested':
            # check if there are nested files
            if are_files_nested(input_dir):
                file_levels = 'nested'
            # if no nested files, check if top level files
            elif are_files_top_level(input_dir):
                logger.info(f"--file_levels is set to {args.file_levels} but {input_dir} does not contain nested files. Running on top level files.")
                file_levels = 'top'
            else:
                # if not top level or nested files, return error
                logger.error(f"{input_dir} doesn't contain files with fastq, fq, fastq.gz, or fq.gz extensions")
                sys.exit(1)
        elif args.file_levels == 'top':
            if are_files_top_level(input_dir):
                file_levels = 'top'
            else:
                logger.error(f"--file_levels is set to {args.file_levels} but {input_dir} doesn't contain top-level files with fastq, fq, fastq.gz, or fq.gz extensions.")
                sys.exit(1)   
    else:
        file_levels = 'all'


    # Get the list of samples
    if args.single:
        samples_dict = list_samples(input_dir, file_levels = file_levels, single=True)
    else:
        samples_dict = list_samples(input_dir, file_levels = file_levels)

    # Error if files are nested in a way the script doesn't handle
    if args.single:
        create_samplesheet(samples_dict, output_dir, outfile, single=True)
    else:
        create_samplesheet(samples_dict, output_dir, outfile)
    logger.info(f"{outfile} created. Please review and be sure it is correct.")

if __name__ == '__main__':
    main()
