import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from bioservices.kegg import KEGG

data = []

df = pd.read_csv('RPKMs.csv',delimiter=";") 
    
#print(df["M381.Control"][0:3])

k = KEGG()

#k.list("organism")
#
#print(k.organismIds())

for i in range(100):
    print(i,"****")
#    print(k.find("hsa", str(df["symbol"][i])))
    print("//\n",k.get_pathway_by_gene(str(df["symbol"][i]), "hsa"))

#print(k.find("hsa", "a1bg"))




 