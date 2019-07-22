#!/usr/bin/python

''' 
Script to create scratch files and submit jobs for tophat mapping of the RNA-seq data.
Description- Creates and submit pbs files for all the .fq files in current directory

Copy do_tophatmapping.py to the current folder.

USAGE: python do_tophatmapping.py -p /path/rnaseqscripts
rnaseqscripts folder that contains the indexes and Araport11 gff file

1) Sample name is tophat output directory
2) Can also handle / after rnascripts path
3) Included log file 
4) Time and memory options
5) Log file keeps record of sample to job id mapping
6) Accepts user provided gff file
7) If the tophat output dir already exists for a sample and has align_summary file, it does not rerun the job for this sample
'''

import sys,os
import shutil
import glob
import subprocess as sub
import argparse as agp

def domapping(scriptpath,gfffile,time,memory):

	logfile= open('log_tophat.txt','wb')

	if scriptpath[-1]=='/':
		scriptpath= scriptpath[:-1]

	datadir= os.getcwd()
	print 'Current working directory= ',datadir
	
	if not gfffile:
		gfffile=scriptpath+'/Araport11_GFF3_genes_transposons.201606.gff'
	if not time:
		time='03:00:00'
	if not memory:
		memory='20'

	for myfile in glob.glob('*.fq'):
		print 'Sample= ',myfile
		sbatch= open(datadir+'/sbatch_tophatmapping.sh','w')
		sbatch.write('#!/bin/sh\n')
		sbatch.write('#\n')
		sbatch.write('#SBATCH --verbose\n')
		sbatch.write('#SBATCH --job-name=run_Tophat\n')
		sbatch.write('#SBATCH --output=tophat_%j.out\n')
		sbatch.write('#SBATCH --error=tophat_%j.err\n')
		sbatch.write('#SBATCH --time='+time+'\n')
		sbatch.write('#SBATCH --nodes=1\n')
		sbatch.write('#SBATCH --cpus-per-task=4\n')
		sbatch.write('#SBATCH --mem='+str(memory)+'GB\n')
		sbatch.write('\n')
		sbatch.write('cd '+datadir+'\n')
		sbatch.write('module load bowtie2/intel/2.2.9\n')
		sbatch.write('module load tophat/intel/2.1.1\n')
		resdir= myfile.split('.')[0]
		#print 'resdir= ',resdir
		flag=0
		# checks if resdir already exists and has align summary file (which indicates this sample already has output). Don't rerun
		if resdir:
			filesindir= glob.glob(resdir+"/*.txt")
			for val_fl in filesindir:
				if 'align_summary.txt' in val_fl:
					flag=1
		# if the align_summary.txt exists, don't run tophat command for this sample. Continue with the next sample.
		if flag==1:
			print('\n Warning: Tophat not reunning for sample ',myfile)
			print('output directory already exists for this sample\nIf you would like to rerun tophat for this sample, delete the directory\n')
			logfile.write('\n Warning: Tophat not reunning for sample '+myfile+' \n')
			logfile.write('Output directory already exists for this sample. If you would like to rerun tophat for this sample, delete the directory\n')
			continue	
		sbatch.write('mkdir '+resdir+'\n')
		mycmd= ''.join(['tophat -o ',resdir,' -p 4 -G ',gfffile,' ',scriptpath,'/bowtie_index_Ath/TAIR10_Ath_Allchr ',myfile.strip()])
		sbatch.write(mycmd)
		sbatch.write('\n')
		sbatch.close()
		shutil.copyfileobj(open("sbatch_tophatmapping.sh",'rb'), logfile)
		p= sub.Popen(["sbatch", "sbatch_tophatmapping.sh"], stdout=sub.PIPE, stderr=sub.PIPE)
		output, errors = p.communicate()
		logfile.write('Sample '+myfile.strip()+', '+output.strip())
		logfile.write('\n\n')
		
if __name__=='__main__':
	parser= agp.ArgumentParser()
	parser.add_argument('-p', help='Provide path for scripts', required=True)
	parser.add_argument('-g', help='Provide your own  gff file: default-Araport11_GFF3_genes_transposons.201606.gff')
	parser.add_argument('-t', help='Time Limit: default 3hrs- specify in the format -t 03:00:00')
	parser.add_argument('-m', help='Memory Limit: default 20GB- specify in the format -m 20')
	args = parser.parse_args()
	domapping(args.p,args.g,args.t,args.m)

