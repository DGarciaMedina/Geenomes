import numpy as np
from numpy import genfromtxt
import csv
import pandas as pd

in_file = "RPKMs.csv"

def main():
	genes, data = process_data(in_file)
	print genes

def process_data(file):
	df=pd.read_csv(file,sep=',')
	gene_names = df.values[:,0]
	data = df.values[:,1:]
	all_data = []
	for i in range(data.shape[1]/3):
		all_data.append(data[:,3*i:3*i+2])
	
	return gene_names, all_data
		
if __name__ == '__main__':
	main()
