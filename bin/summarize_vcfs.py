#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import subprocess as sp

#define the user inputs 
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v1", "--vcf1", required=True, help="REQUIRED: The path to the ivar vcf.")
    parser.add_argument("-v2", "--vcf2", required=True, help="REQUIRED: The full path to the lofreq vcf.")
    parser.add_argument("-s", "--sample", required=True, help="REQUIRED: The sample name.")
    parser.add_argument("-exc", "--exclude_coords", required=True, help="REQUIRED: The full path to list of coordinates to exclude, which is a simple text file with coords like 145689")
    parser.add_argument("-af", "--af_cutoff", required=True, help="REQUIRED: Lower bounds of allele frequency cutoff for variant filtering")
    parser.add_argument("-min", "--min_reads", required=True, help="REQUIRED: Minimum reads required to retain variant when filtering")
    parser.add_argument("-r", "--run", required=True, help="REQUIRED: The name of the sequencing run, likely the same as the outdir name")
    return parser.parse_args()

def summarize_vcfs(vcf1,vcf2,sample,exclude_coords,af_cutoff,min_reads,run):
    
    #set outfiles
    summary_tmp="vcf.summary.tmp"
    fh_tmp=open(summary_tmp,'w')

    #create lists for coords to exclude
    fh_coords=open(exclude_coords,'r')
    exclude=[]
    for line in fh_coords:
        line=line.strip()
        exclude.append(line)

    #create lists for unique sites for each caller
    vcfs = [vcf1,vcf2]
    iv_uniq=[]
    lf_uniq=[]

    #now let's parse the vcfs and write to a tmp file
    for v in vcfs:
        vcf=open(v,'r')
        caller=v.split('.')[1]
        for line in vcf:
            if not line.startswith('#'):
                line=line.split('\t')
                coord=line[1]
                if coord not in exclude:
                    if caller == "ivar":
                        ref_al=line[3]
                        alt_al=line[4]
                        depth=int(line[9].split(":")[4])
                        freq=line[9].split(":")[7]
                        uniq=sample+':'+coord
                        iv_uniq.append(uniq)
                        if float(freq) >= float(af_cutoff) and depth > int(min_reads):
                            fh_tmp.write(
                            caller +'\t'+ \
                            run +'\t'+ \
                            sample+'\t'+ \
                            str(coord)+'\t'+\
                            ref_al+'\t'+\
                            alt_al+'\t'+\
                            str(depth)+'\t'+\
                            str(freq))
                        elif caller == "lofreq":
                            coord=line[1]
                            ref_al=line[3]
                            alt_al=line[4]
                            info=line[7]
                            '''
                            In lofreq the depth for alt allele is given in the INFO col as DP4
                            where DP4=forward ref allele, reverse ref allele, forward alt allele, reverse alt allele
                            so the count we want is the sum of the third and fourth value
                            '''
                            depth=str(int(info.split(';')[3].split('=')[1].split(',')[2])+\
                                int(info.split(';')[3].split('=')[1].split(',')[3]))
                            freq=info.split(';')[1].split('=')[1]
                            uniq=sample+':'+coord
                            lf_uniq.append(uniq)
                            #filter out everything below min af cutoff
                            if float(freq) >= float(af_cutoff) and depth > int(min_reads):
                                fh_tmp.write(
                                caller +'\t'+ \
                                run +'\t'+ \
                                sample+'\t'+ \
                                coord+'\t'+\
                                ref_al+'\t'+\
                                alt_al+'\t'+\
                                depth+'\t'+\
                                freq+'\t'+'\n'
                                )

    #open the tmp table for reading
    fh_tmp.close()
    fh_tmp=open(summary_tmp,'r')

    #open the final output files
    summary=sample+".summary_per_caller.txt"
    final=sample+'.vcf.summary.txt'

    fh_summary=open(summary,'w')
    fh_out=open(final,'w')
    fh_out.write('Caller'+'\t'+\
        "Run"+'\t'+\
        "sample"+'\t'+\
        "Coordinate"+'\t'+\
        "Ref_Allele"+'\t'+\
        "Alt_Allele"+'\t'+\
        "Alt_Read_Depth"+'\t'+\
        "Alt_Freq"+"\n")

    #print to summary out file
    total_vars=len(set(iv_uniq+lf_uniq))
    unique_ivar = ((set(iv_uniq)-set(lf_uniq)))
    unique_lofreq = ((set(lf_uniq)-set(iv_uniq)))

    fh_summary.write('# potentially drug resistant variants'+'\t'+\
                str(total_vars) + '\n' + \
                'variants unique to iVar'+'\t'+ \
                str(len(unique_ivar))+'\n'\
                "variants unique to LoFreq"+'\t'+ \
                str(len(unique_lofreq))+'\n'
                )
                
    #count vars per locus
    F13L_end=40204
    E9L_start=51427

    F13L=[]
    E9L=[]
    indV=[]
    indV_ivar=[]
    indV_lofreq=[]

    #now write the variant table out
    for line in fh_tmp:
        line=line.strip()
        caller=line.split('\t')[0]
        run=line.split('\t')[1]
        sample=line.split('\t')[2]
        coord=line.split('\t')[3]
        if caller == "ivar":
            fh_out.write(line+'\n')
            indV.append(sample)
            indV_ivar.append(sample)
            if int(coord) <= int(F13L_end):
                F13L.append(coord)
            elif int(coord) >= int(E9L_start):
                E9L.append(coord)
        elif caller == "lofreq":
            if sample+':'+coord in unique_lofreq:
                fh_out.write(line+'\n')
                indV.append(sample)
                indV_lofreq.append(sample)
    fh_summary.write('number of variant sites in F13L'+'\t'+\
                str(len(set(F13L)))+'\n'+\
                'number of variant sites in E9L'+'\t'+ \
                str(len(set(E9L))))

    #clean up
    fh_out.close()
    fh_summary.close()

def main():
	#define the arguments
	args = get_args()
	summarize_vcfs(args.vcf1,args.vcf2,args.sample,args.exclude_coords,args.af_cutoff,args.min_reads,args.run)

if __name__ == '__main__':
    main()
