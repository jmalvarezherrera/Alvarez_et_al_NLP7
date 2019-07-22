#!/usr/bin/python

''' 
This script sorts all the bam files and runs GenomicFeaures script to create a read count file 'counts_combined.csv'. Also creates another output file
'filtered_counts_combined.csv' that has low expressed genes removed. A gene will be discarded, if none of the
treatment groups has median reads greater than 10.

Description: 
Creates and submit scratch files for all the .bam files in current directory

Copy do_readcounts.py, discard_lowexpressed.py, genomicfeautures.R and designfile.txt to the current folder with the .bam files. 
A sample_designfile.txt is provided. 

USAGE: python do_readcounts.py -d designfile.txt -p /path/rnaseqscripts -c 10

rnaseqscripts folder contains the Araport11 gff file. 

1) A log file is included
2) Run genomicfeatures.R to counts reads per gene
3) Included a python script low-expressed genes to filter genes based on treatment groups 
4) User can provide cutoff for discarding low-expressed genes
5) Can now work on full path for rnaseqcripts folder
6) Uses python3
'''

import sys,os,shutil
import glob
from subprocess import call
import argparse as agp

def do_readcounts(designfile,scriptpath,cutoff):

        logfile= open('log_readcounts.txt','wb')

	if scriptpath[-1]=='/':
                scriptpath= scriptpath[:-1]

	if not cutoff:
                cutoff=10

	datadir= os.getcwd()
        print 'Current working directory= ',datadir
	sbatch= open(datadir+'/sbatch_readcounts.sh','w')
	sbatch.write('#!/bin/sh\n')
	sbatch.write('#\n')
        sbatch.write('#SBATCH --verbose\n')
        sbatch.write('#SBATCH --job-name=run_readcounts\n')
        sbatch.write('#SBATCH --output=readcount_%j.out\n')
        sbatch.write('#SBATCH --error=readcount_%j.err\n')
        sbatch.write('#SBATCH --time=02:00:00\n')
        sbatch.write('#SBATCH --nodes=1\n')
        sbatch.write('#SBATCH --cpus-per-task=1\n')
        sbatch.write('#SBATCH --mem=20GB\n')
        sbatch.write('\n')
        sbatch.write('cd '+datadir+'\n')
        sbatch.write('module load r/intel/3.3.2\n')
	sbatch.write('module load samtools/intel/1.3.1\n')
	sbatch.write('module load python3/intel/3.5.3\n')
	sbatch.write('mkdir allbamfiles\n')
	with open(designfile) as design:
		next(design)
		for val_d in design:
			sbatch.write('cp '+(val_d.split('\t')[1]+'/accepted_hits.bam ')+('allbamfiles/'+val_d.split('\t')[1]+'_accepted_hits.bam')+'\n')
			sbatch.write('samtools sort '+('allbamfiles/'+val_d.split('\t')[1]+'_accepted_hits.bam')+' -o '+('allbamfiles/'+val_d.split('\t')[1]+'_accepted_hits.sorted.bam')+'\n')
	sbatch.write('rm allbamfiles/*_accepted_hits.bam\n')
	sbatch.write('cp '+scriptpath+'/genomicfeatures.R . \n')
	sbatch.write('Rscript genomicfeatures.R '+scriptpath+' \n')
	sbatch.write('cp '+scriptpath+'/discard_lowexpressed.py . \n')
	sbatch.write(''.join(['python3 discard_lowexpressed.py -c counts_combined.csv -d ',designfile,' -p ',str(cutoff),'\n']))
	sbatch.write('\n')
	sbatch.close()
	call(["sbatch", "sbatch_readcounts.sh"])
	shutil.copyfileobj(open("sbatch_readcounts.sh",'rb'), logfile)


if __name__=='__main__':
        parser= agp.ArgumentParser()
        parser.add_argument('-d',help='Provide design file for your data',required=True)
	parser.add_argument('-p',help='Provide path for scripts',required=True)
        parser.add_argument('-c',help='Provide Median cutoff for discarding low-expressed genes')
	args = parser.parse_args()
        do_readcounts(args.d,args.p,args.c)
