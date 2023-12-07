#!/usr/bin/env python

import argparse
import os
import subprocess as sp

#define the user inputs 
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tsv", required=True, help="REQUIRED: The path to the ivar tsv output.")
    parser.add_argument("-s", "--sample", required=True, help="REQUIRED: The sample name.")
    parser.add_argument("-af", "--af_cutoff", required=True, help="REQUIRED: Lower bounds of allele frequency cutoff for variant filtering")
    parser.add_argument("-DP", "--DP_cutoff", required=True, help="REQUIRED: Lower bounds of total read depth at site for variant filtering")
    parser.add_argument("-AD", "--altDP_cutoff", required=True, help="REQUIRED: Lower bounds of alternative allele read depth at site for variant filtering")
    parser.add_argument("-aq", "--altQ_cutoff", required=True, help="REQUIRED: Lower bounds of quality score for variant site for variant filtering")
    parser.add_argument("-r", "--run", required=True, help="REQUIRED: The name of the sequencing run, likely the same as the outdir name")
    parser.add_argument("-fs", "--F13L_start", required=True, help="REQUIRED: The starting coordinates for gene F13L")
    parser.add_argument("-fe", "--F13L_end", required=True, help="REQUIRED: The ending coordinates for gene F13L")
    parser.add_argument("-es", "--E9L_start", required=True, help="REQUIRED: The starting coordinates for gene E9L")
    parser.add_argument("-ee", "--E9L_end", required=True, help="REQUIRED: The ending coordinates for gene E9L")

    return parser.parse_args()

def summarize_tsv(tsv,run,sample,F13L_start,F13L_end,E9L_start,E9L_end,af_cutoff,DP_cutoff,altDP_cutoff,altQ_cutoff):
    
    # define outfile and write the header
    file_out=sample + '_ivar_summary.txt'
    fh_out=open(file_out,'w')
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
    #now open the ivar tsv file
    fh=open(tsv,'r')
    sample=tsv.split('.')[0]
    #iterate through and define the params
    for line in fh.readlines()[1:]:
        line=line.strip()
        i=line.split('\t')
        coords=i[1]
        refBase=i[2]
        altBase=i[3]
        totalDP=i[11]
        altDP=i[7]
        altQual=i[9]
        altFreq=i[10]
        altPval=i[12]
        GFFfeature=i[14]
        RefCodon=i[15]
        RefAA=i[16]
        AltCodon=i[17]
        AltAA=i[18]
        '''
        apply position filtering for F13L
        default values for filtering are: 
            if > or = to 
            5x total depth
            3x alt allele depth
            1% allele freq
            20 quality score
        '''
        #if GFFfeature== 'cds-QNP13638.1' \
        if int(coords) >= int(F13L_start) and int(coords) <= int(F13L_end) \
            and int(altDP) >= int(altDP_cutoff) \
            and int(totalDP) >= int(DP_cutoff) \
            and float(altFreq) >= float(af_cutoff) \
            and int(altQual) >= int(altQ_cutoff):
            fh_out.write(
                run + '\t' + \
                sample + '\t' + \
                coords + '\t' + \
                refBase + '\t' + \
                altBase + '\t' + \
                totalDP + '\t' + \
                altDP + '\t' + \
                altQual + '\t' + \
                altFreq + '\t' + \
                altPval + '\t' + \
                GFFfeature + '\t' + \
                RefCodon + '\t' + \
                RefAA + '\t' + \
                AltCodon + '\t' + \
                AltAA + '\n')
        # apply position filtering for E9L
        # same filtering regime as above
        #elif GFFfeature == 'cds-QNP13652.1' \
        elif int(coords) >= int(E9L_start) and int(coords) <= int(E9L_end) \
            and int(altDP) >= int(altDP_cutoff) \
            and int(totalDP) >= int(DP_cutoff) \
            and float(altFreq) >= float(af_cutoff) \
            and int(altQual) >= int(altQ_cutoff):
            fh_out.write(
                run + '\t' + \
                sample + '\t' + \
                coords + '\t' + \
                refBase + '\t' + \
                altBase + '\t' + \
                totalDP + '\t' + \
                altDP + '\t' + \
                altQual + '\t' + \
                altFreq + '\t' + \
                altPval + '\t' + \
                GFFfeature + '\t' + \
                RefCodon + '\t' + \
                RefAA + '\t' + \
                AltCodon + '\t' + \
                AltAA + '\n')

def main():
    args = get_args()
    summarize_tsv(args.tsv,args.run,args.sample,\
                  args.F13L_start,args.F13L_end,\
                  args.E9L_start,args.E9L_end,\
                  args.af_cutoff,args.DP_cutoff,\
                  args.altDP_cutoff,args.altQ_cutoff)

if __name__ == '__main__':
    main()
