import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from bioservices.kegg import KEGG

np.set_printoptions(threshold=np.nan)

data = []

df = pd.read_csv('RPKMs.csv',delimiter=",") 

k = KEGG()

#for i in range(100):
#    print(i,"****")
#    print("//\n",k.get_pathway_by_gene(str(df["symbol"][i]), "hsa"))
    
    
def search_pathways_4_list(list_of_genes):
    
    matrix = [[0 for j in range(len(list_of_genes))] for i in range(0)]
    list_of_pathways = []
    dict_of_genes = {}
    
    for i,gene in enumerate(list_of_genes):
        try:
            pathways = k.get_pathway_by_gene(gene, "hsa")
            
            if pathways != None:
                pathways = pathways.values()
                
                dict_of_genes[gene] = list(pathways)
            
                for pathway in pathways:
                    
                    if pathway not in list_of_pathways:
                        list_of_pathways.append(pathway)
                        
                        matrix.append([0 for j in range(len(list_of_genes))])
                        
                        
                    matrix[list_of_pathways.index(pathway)][i] = 1
            
                
        except Exception as e:
            print(e)
            
    return matrix,list_of_pathways,dict_of_genes


list_of_genes = df["symbol"][:10]

matrix, list_of_pathways, dict_of_genes = search_pathways_4_list(list_of_genes)

print(np.array(matrix))
print(dict_of_genes)


 