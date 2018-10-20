from preprocess import process_data
import numpy as np
import pickle

in_file = "RPKMs.csv"
thresh = 0.2

def test():

	genes, data = process_data(in_file)
	print genes
	print data

def main():
	
	all_res_genes = []
	
	genes, data = process_data(in_file)
	for i in range(len(data)):
		res_idx = get_res_genes(data[i])
		res_genes = []
		for j in range(len(res_idx)):
			res_genes.append(genes[j])
		print "resistant genes found: {}".format(len(res_genes))
		all_res_genes.append(res_genes)
	
	res_gene_list = all_res_genes[0]
	for i in range(len(all_res_genes)):
		res_gene_list = [x for x in res_gene_list if x in all_res_genes[i]]
		
	pickle.dump(res_gene_list, "res_gene_list.txt")    
	print "Pickled as res_gene_list.txt"
        
	print res_gene_list
	print len(res_gene_list)
	
def get_res_genes(cell_line):
	grad1 = cell_line[:,1] - cell_line[:,0]
	grad2 = cell_line[:,2] - cell_line[:,1]
	grad1_idx = [idx for idx in range(len(grad1)) if grad1[idx]/cell_line[idx,1]>thresh]
	grad2_idx = [idx for idx in range(len(grad2)) if grad2[idx]/cell_line[idx,2]<-thresh]
	res_idx = [x for x in grad1_idx if x in grad2_idx]
	
	return res_idx

if __name__ == "__main__":
	main()
