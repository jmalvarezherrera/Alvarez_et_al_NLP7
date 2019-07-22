#!/usr/bin/python
''' 
This script is to rename files, to create sbatch files and submit job for quality filtering of data.
Description- Creates and submit pbs files for all the .fq files in current directory.
Copy qualityfiltering.py and designfile.txt to the current folder with the fast.gz files. A sample_designfile.txt is provided. 
USAGE: python qualityfiltering.py -d designfile.txt -f fastq.gz

1) Can work with fastq.gz or fastq files
2) Changed the design format
3) Submits one job per sample
4) Checks if all the files (with user specified format) exists in the design file. If not the program exits without submitting the job
5) Time and memory options
6) Log file keeps record of sample to job id mapping
7) If qualityfiltered folder already exists, the program exits
'''

import sys,os
import shutil
import glob
import subprocess as sub
import argparse as agp


def qualityfilteringPBS(designfile, fileformat, time, memory):

	renamedict= dict()
	with open(designfile) as design:
    		next(design)
		for val_d in design:
			renamedict[val_d.split('\t')[0].strip()]= val_d.split('\t')[1].strip()
	
	if not time:
		time='03:00:00'
	if not memory:
		memory=10
	
	logfile= open('log_qualfilter.txt','wb')	

        datadir= os.getcwd()
        print('Current working directory= ',datadir)
	
	# Verify if design file is correct-
        ## check if all sample files with specified format are included in the design file.
	## if not then display error not all the files exist and exit the program.
	allfiles= glob.glob('*.'+fileformat)
	if set(allfiles)==set(list(renamedict.keys())):
		if not os.path.exists('qualityfiltered'):
                	os.makedirs('qualityfiltered')
        	else:
                	sys.exit('Error- qualityfiltered folder already exists!\n Program Exit')
		for myfile in glob.glob('*.'+fileformat):
                	print('myfile= ',myfile)
			sbatch= open(datadir+'/sbatch_qualitytrim.sh','w')
        		sbatch.write('#!/bin/sh\n')
        		sbatch.write('#\n')
        		sbatch.write('#SBATCH --verbose\n')
        		sbatch.write('#SBATCH --job-name=qualityfilter\n')
        		sbatch.write('#SBATCH --output=qualfilter_'+renamedict[myfile].strip()+'_%j.out\n')
        		sbatch.write('#SBATCH --error=qualfilter_'+renamedict[myfile].strip()+'_%j.err\n')
        		sbatch.write('#SBATCH --time='+time+'\n')
        		sbatch.write('#SBATCH --nodes=1\n')
        		sbatch.write('#SBATCH --cpus-per-task=1\n')
        		sbatch.write('#SBATCH --mem='+str(memory)+'GB\n')
       			sbatch.write('\n')
        		sbatch.write('module load fastx_toolkit/intel/0.0.14\n')
        		sbatch.write('cd '+datadir+'\n')
			cmd_rename= ''.join(['mv ',myfile.strip(),' ',renamedict[myfile].strip()+'.'+fileformat])
			sbatch.write(cmd_rename+'\n')
			if '.gz' in fileformat: # check if file is zipped, unzip if zipped
				sbatch.write('gunzip '+(renamedict[myfile].strip()+'.'+fileformat)+'\n')
                		myfile_unzip= (renamedict[myfile]+'.'+fileformat)[:-3]
			else:
				myfile_unzip= renamedict[myfile]+'.'+fileformat
			trimclipout= 'qualityfiltered/'+renamedict[myfile].strip().split('.')[0]+'.clip.fq'
                	mycmd1= ''.join(['fastx_clipper -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -l 20 -v -i ',myfile_unzip.strip(),' -o ',trimclipout,'\n'])
                	sbatch.write(mycmd1)
			sbatch.write('\n')
			sbatch.close()
			shutil.copyfileobj(open("sbatch_qualitytrim.sh",'rb'), logfile)
			p= sub.Popen(["sbatch", "sbatch_qualitytrim.sh"], stdout=sub.PIPE, stderr=sub.PIPE)
			output, errors = p.communicate()
			logfile.write('Sample '+renamedict[myfile].strip()+', '+output.strip())
			logfile.write('\n\n')
	else:
		print('\nError: Can not submit jobs for this step')
		print('File(s) missing from the design file= ',list(set(allfiles)-set(list(renamedict.keys()))))		
		print('Include these file(s) in the design file or remove from the current folder\n')
if __name__=='__main__':
        parser= agp.ArgumentParser()
	parser.add_argument('-d', help='Provide design file for your data', required=True)
	parser.add_argument('-f', help='Provide your file format', required=True)
	parser.add_argument('-t', help='Time Limit (default 3hrs, specify in the format 03:00:00)')
	parser.add_argument('-m', help='Memory Limit (default 10GB, specify in the format 10)')
	args = parser.parse_args()
	qualityfilteringPBS(args.d,args.f,args.t,args.m)

