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
    for i in range(int(data.shape[1]/3)):
        all_data.append(data[:,3*i:3*i+3])
    
    return gene_names, all_data

def pick_data(file, cell_lines):
    df=pd.read_csv(file,sep=',')
    gene_names = df.values[:,0]
    line_names = list(df.columns.values)
    #print(line_names)
    data = df.values[:,1:]
    all_data = []
    for i in range(int(data.shape[1]/3)):
        #print(line_names[3*i+1][:4])
        if (line_names[3*i+1][:4]) in cell_lines:
            all_data.append(data[:,3*i:3*i+3])
    
    return gene_names, all_data
    

if __name__ == "__main__":
    in_file = "RPKMs.csv"
    genes, data = pick_data(in_file, ["M263"])