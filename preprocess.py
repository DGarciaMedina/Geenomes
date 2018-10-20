'''
given file path, get list of names and list of np arrays containing data for each cell line
'''
import numpy as np
import pandas as pd

def process_data(file):
	df=pd.read_csv(file,sep=',')
	gene_names = df.values[:,0]
	data = df.values[:,1:]
	all_data = []
	for i in range(data.shape[1]/3):
		all_data.append(data[:,3*i:3*i+3])
	
	return gene_names, all_data
