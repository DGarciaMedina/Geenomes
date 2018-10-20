import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from bioservices.kegg import KEGG

data = []

df = pd.read_csv('RPKMs.csv',delimiter=";") 

k = KEGG()

#for i in range(100):
#    print(i,"****")
#    print("//\n",k.get_pathway_by_gene(str(df["symbol"][i]), "hsa"))
    
    
def search_pathways_4_list(list_of_genes):
    for i,gene in enumerate(list_of_genes):
        try:
            print(i,"****")
            print("//\n",k.get_pathway_by_gene(gene, "hsa"))
        except Exception as e:
            print(e)


search_pathways_4_list(df["symbol"][:100])



 