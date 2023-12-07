#!/usr/bin/env python

import argparse
import os
import subprocess as sp

#define the user inputs 
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: path to directory with sequence data")
    return parser.parse_args()

def generate_samples():
	#write the outfile for all samples
	outfile='all_samples.vcf.summary.txt'
	if os.path.exists(outfile):
		os.remove(outfile)
	fh_out=open(outfile,'w')
	fh_out.write(
        'run' + '\t' + \
        'sample'+ '\t' + \
        'coord' + '\t' + \
        'refBase' + '\t' + \
        'altBase' + '\t' + \
        'totalDP' + '\t' + \
        'altDP' + '\t' + \
        'altQual' + '\t' + \
        'altFreq' + '\t' + \
        'altPval' + '\t' + \
        'GFFfeature' + '\t' + \
    	'RefCodon' + '\t' + \
        'RefAA' + '\t' + \
        'AltCodon' + '\t' + \
        'AltAA' + '\n')
	for f in os.listdir('.'):
		if f.endswith('_ivar_summary.txt'):
			fh=open(f,'r')
			for line in fh.readlines()[1:]:
				line=line.strip()
				fh_out.write(line+'\n')

def main():
	#define the arguments
	args = get_args()
	combine_summaries()

if __name__ == '__main__':
    main()

