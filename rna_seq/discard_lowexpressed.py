#!/usr/bin/python3

'''
This script discards the low expressed genes. It is executed as part of the do_readcounts.py script 

'''

import sys,os
import glob
from subprocess import call
import argparse as agp
import pandas as pd
from collections import defaultdict

def discard_lowexp(countcombined, designfile, cutoff):

	treat_dict= dict()
	group_dict= defaultdict(list)

	with open(designfile,'r') as design:
		next(design)
		for val_d in design:
			treat_dict[val_d.split('\t')[1]]= val_d.split('\t')[2]
			group_dict[val_d.split('\t')[1]]= val_d.split('\t')[2]+'_'+val_d.split('\t')[2]
	
	#print('group_dict= ',group_dict)

	pdframe= pd.read_table(countcombined, sep=',', index_col=0, header=0)
	#print('pdframe= ',pdframe.columns)
	
	pdframe.rename(columns=lambda x: x.split('_accepted')[0], inplace=True)
	#print('after= ',pdframe.columns)

	maskdf = pdframe.groupby(group_dict, axis=1).median().gt(int(cutoff)).any(1)
	
	outputfile= 'filtered_counts_combined.csv'
	pdframe[maskdf].to_csv(outputfile,sep=',')
	#outputfile.close()
		
if __name__=='__main__':
	parser= agp.ArgumentParser()
	parser.add_argument('-c',help='counts file',required=True)
	parser.add_argument('-d',help='design file',required=True)
	parser.add_argument('-p',help='Median cutoff for discarding low expressed genes')
	args = parser.parse_args()
	discard_lowexp(args.c, args.d, args.p)

