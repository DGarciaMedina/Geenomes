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

list_of_genes = ['A1BG', 'A1CF', 'A2M', 'A2ML1', 'A2MP1', 'A3GALT2', 'A4GALT', 'A4GNT', 'AA06', 'AAAS', 'AACS', 'AACSP1', 'AADAC', 'AADACL2', 'AADACL3', 'AADACL4', 'AADACP1', 'AADAT', 'AAED1', 'AAGAB', 'AAK1', 'AAMDC', 'AAMP', 'AANAT', 'AAR2', 'AARD', 'AARS', 'AARS2', 'AARSD1', 'AASDH', 'AASDHPPT', 'AASS', 'AATBC', 'AATF', 'AATK', 'ABAT', 'ABCA1', 'ABCA10', 'ABCA11P', 'ABCA12', 'ABCA13', 'ABCA17P', 'ABCA2', 'ABCA3', 'ABCA4', 'ABCA5', 'ABCA6', 'ABCA7', 'ABCA8', 'ABCA9', 'ABCB1', 'ABCB10', 'ABCB11', 'ABCB4', 'ABCB5', 'ABCB6', 'ABCB7', 'ABCB8', 'ABCB9', 'ABCC1', 'ABCC10', 'ABCC11', 'ABCC12', 'ABCC13', 'ABCC2', 'ABCC3', 'ABCC4', 'ABCC5', 'ABCC6', 'ABCC6P1', 'ABCC6P2', 'ABCC8', 'ABCC9', 'ABCD1', 'ABCD2', 'ABCD3', 'ABCD4', 'ABCE1', 'ABCF1', 'ABCF2', 'ABCF3', 'ABCG1', 'ABCG2', 'ABCG4', 'ABCG5', 'ABCG8', 'ABHD1', 'ABHD10', 'ABHD11', 'ABHD12', 'ABHD12B', 'ABHD13', 'ABHD14A', 'ABHD14A-ACY1', 'ABHD14B', 'ABHD15', 'ABHD16A', 'ABHD16B', 'ABHD17A', 'ABHD17B', 'ABHD17C', 'ABHD18', 'ABHD2', 'ABHD3', 'ABHD4', 'ABHD5', 'ABHD6', 'ABHD8', 'ABI1', 'ABI2', 'ABI3', 'ABI3BP', 'ABL1', 'ABL2', 'ABLIM1', 'ABLIM2', 'ABLIM3', 'ABO', 'ABR', 'ABRA', 'ABRACL', 'ABT1', 'ABTB1', 'ABTB2', 'ACAA1', 'ACAA2', 'ACACA', 'ACACB', 'ACAD10', 'ACAD11', 'ACAD8', 'ACAD9', 'ACADL', 'ACADM', 'ACADS', 'ACADSB', 'ACADVL', 'ACAN', 'ACAP1', 'ACAP2', 'ACAP3', 'ACAT1', 'ACAT2', 'ACBD3', 'ACBD4', 'ACBD5', 'ACBD6', 'ACBD7', 'ACCS', 'ACCSL', 'ACD', 'ACE', 'ACE2', 'ACER1', 'ACER2', 'ACER3', 'ACHE', 'ACIN1', 'ACKR1', 'ACKR2', 'ACKR3', 'ACKR4', 'ACLY', 'ACMSD', 'ACO1', 'ACO2', 'ACOD1', 'ACOT1', 'ACOT11', 'ACOT12', 'ACOT13', 'ACOT2', 'ACOT4', 'ACOT6', 'ACOT7', 'ACOT8', 'ACOT9', 'ACOX1', 'ACOX2', 'ACOX3', 'ACOXL', 'ACP1', 'ACP2', 'ACP5', 'ACP6', 'ACP7', 'ACPP', 'ACPT', 'ACR', 'ACRBP', 'ACRC', 'ACRV1', 'ACSBG1', 'ACSBG2', 'ACSF2', 'ACSF3', 'ACSL1', 'ACSL3', 'ACSL4', 'ACSL5', 'ACSL6', 'ACSM1', 'ACSM2A', 'ACSM2B', 'ACSM3', 'ACSM4', 'ACSM5', 'ACSM6', 'ACSS1', 'ACSS2', 'ACSS3', 'ACTA1', 'ACTA2', 'ACTB', 'ACTBL2', 'ACTC1', 'ACTG1', 'ACTG1P17', 'ACTG1P20', 'ACTG1P4', 'ACTG2', 'ACTL10', 'ACTL6A', 'ACTL6B', 'ACTL7A', 'ACTL7B', 'ACTL8', 'ACTL9', 'ACTN1', 'ACTN2', 'ACTN3', 'ACTN4', 'ACTR10', 'ACTR1A', 'ACTR1B', 'ACTR2', 'ACTR3', 'ACTR3B', 'ACTR3BP2', 'ACTR3BP5', 'ACTR3C', 'ACTR5', 'ACTR6', 'ACTR8', 'ACTRT1', 'ACTRT2', 'ACTRT3', 'ACVR1', 'ACVR1B', 'ACVR1C', 'ACVR2A', 'ACVR2B', 'ACVRL1', 'ACY1', 'ACY3', 'ACYP1', 'ACYP2', 'ADA', 'ADAD1', 'ADAD2', 'ADAL', 'ADAM10', 'ADAM11', 'ADAM12', 'ADAM15', 'ADAM17', 'ADAM18', 'ADAM19', 'ADAM1A', 'ADAM2', 'ADAM20', 'ADAM20P1', 'ADAM21', 'ADAM21P1', 'ADAM22', 'ADAM23', 'ADAM28', 'ADAM29', 'ADAM30', 'ADAM32', 'ADAM33', 'ADAM3A', 'ADAM5', 'ADAM6', 'ADAM7', 'ADAM8', 'ADAM9', 'ADAMDEC1', 'ADAMTS1', 'ADAMTS10', 'ADAMTS12', 'ADAMTS13', 'ADAMTS14', 'ADAMTS15', 'ADAMTS16', 'ADAMTS17', 'ADAMTS18', 'ADAMTS19', 'ADAMTS2', 'ADAMTS20', 'ADAMTS3', 'ADAMTS4', 'ADAMTS5', 'ADAMTS6', 'ADAMTS7', 'ADAMTS7P1', 'ADAMTS8', 'ADAMTS9', 'ADAMTS9-AS2', 'ADAMTSL1', 'ADAMTSL2', 'ADAMTSL3', 'ADAMTSL4', 'ADAMTSL5', 'ADAP1', 'ADAP2', 'ADAR', 'ADARB1', 'ADARB2', 'ADAT1', 'ADAT2', 'ADAT3', 'ADCK1', 'ADCK2', 'ADCK3', 'ADCK4', 'ADCK5', 'ADCY1', 'ADCY10', 'ADCY10P1', 'ADCY2', 'ADCY3', 'ADCY4', 'ADCY5', 'ADCY6', 'ADCY7', 'ADCY8', 'ADCY9', 'ADCYAP1', 'ADCYAP1R1', 'ADD1', 'ADD2', 'ADD3', 'ADGB', 'ADGRA1', 'ADGRA2', 'ADGRA3', 'ADGRB1', 'ADGRB2', 'ADGRB3', 'ADGRD1', 'ADGRD2', 'ADGRE1', 'ADGRE2', 'ADGRE3', 'ADGRE4P', 'ADGRE5', 'ADGRF1', 'ADGRF2', 'ADGRF3', 'ADGRF4', 'ADGRF5', 'ADGRG1', 'ADGRG2', 'ADGRG3', 'ADGRG4', 'ADGRG5', 'ADGRG6', 'ADGRG7', 'ADGRL1', 'ADGRL2', 'ADGRL3', 'ADGRL4', 'ADGRV1', 'ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7', 'ADHFE1', 'ADI1', 'ADIG', 'ADIPOQ', 'ADIPOR1', 'ADIPOR2', 'ADIRF', 'ADK', 'ADM', 'ADM2', 'ADM5', 'ADNP', 'ADNP2', 'ADO', 'ADORA1', 'ADORA2A', 'ADORA2B', 'ADORA3', 'ADPGK', 'ADPRH', 'ADPRHL1', 'ADPRHL2', 'ADPRM', 'ADRA1A', 'ADRA1B', 'ADRA1D', 'ADRA2A', 'ADRA2B', 'ADRA2C', 'ADRB1', 'ADRB2', 'ADRB3', 'ADRBK1', 'ADRBK2', 'ADRM1', 'ADSL', 'ADSS', 'ADSSL1', 'ADTRP', 'AEBP1', 'AEBP2', 'AEN', 'AES', 'AFAP1', 'AFAP1L1', 'AFAP1L2', 'AFF1', 'AFF2', 'AFF3', 'AFF4', 'AFG3L1P', 'AFG3L2', 'AFM', 'AFMID', 'AFP', 'AFTPH', 'AGA', 'AGAP1', 'AGAP11', 'AGAP2', 'AGAP3', 'AGAP4', 'AGAP5', 'AGAP6', 'AGAP7P', 'AGAP9', 'AGBL1', 'AGBL2', 'AGBL3', 'AGBL4', 'AGBL5', 'AGER', 'AGFG1', 'AGFG2', 'AGGF1', 'AGK', 'AGL', 'AGMAT', 'AGMO', 'AGO1', 'AGO2', 'AGO3', 'AGO4', 'AGPAT1', 'AGPAT2', 'AGPAT3', 'AGPAT4', 'AGPAT4-IT1', 'AGPAT5', 'AGPS', 'AGR2', 'AGR3', 'AGRN', 'AGRP', 'AGT', 'AGTPBP1', 'AGTR1', 'AGTR2', 'AGTRAP', 'AGXT', 'AGXT2', 'AHCTF1', 'AHCTF1P1', 'AHCY', 'AHCYL1', 'AHCYL2', 'AHDC1', 'AHI1', 'AHNAK', 'AHNAK2', 'AHR', 'AHRR', 'AHSA1', 'AHSA2', 'AHSG', 'AHSP', 'AICDA', 'AIDA', 'AIF1', 'AIF1L', 'AIFM1', 'AIFM2', 'AIFM3', 'AIG1', 'AIM1', 'AIM1L', 'AIM2', 'AIMP1', 'AIMP2', 'AIP', 'AIPL1', 'AIRE', 'AIRN', 'AJAP1', 'AJUBA', 'AK1', 'AK2', 'AK3', 'AK4', 'AK5', 'AK7', 'AK8', 'AK9', 'AKAP1', 'AKAP10', 'AKAP11', 'AKAP12', 'AKAP13', 'AKAP14', 'AKAP17A', 'AKAP2', 'AKAP3', 'AKAP4', 'AKAP5', 'AKAP6', 'AKAP7', 'AKAP8', 'AKAP8L', 'AKAP9', 'AKIP1', 'AKIRIN1', 'AKIRIN2', 'AKNA', 'AKNAD1', 'AKR1A1', 'AKR1B1', 'AKR1B10', 'AKR1B15', 'AKR1C1', 'AKR1C2', 'AKR1C3', 'AKR1C4', 'AKR1C6P', 'AKR1C8P', 'AKR1D1', 'AKR1E2', 'AKR7A2', 'AKR7A2P1', 'AKR7A3', 'AKR7L', 'AKT1', 'AKT1S1', 'AKT2', 'AKT3', 'AKTIP', 'ALAD', 'ALAS1', 'ALAS2', 'ALB', 'ALCAM', 'ALDH16A1', 'ALDH18A1', 'ALDH1A1', 'ALDH1A2', 'ALDH1A3', 'ALDH1B1', 'ALDH1L1', 'ALDH1L1-AS2', 'ALDH1L2', 'ALDH2', 'ALDH3A1', 'ALDH3A2', 'ALDH3B1', 'ALDH3B2', 'ALDH4A1', 'ALDH5A1', 'ALDH6A1', 'ALDH7A1', 'ALDH8A1', 'ALDH9A1', 'ALDOA', 'ALDOB', 'ALDOC', 'ALG1', 'ALG10', 'ALG10B', 'ALG11', 'ALG12', 'ALG13', 'ALG14', 'ALG1L', 'ALG1L2', 'ALG1L9P', 'ALG2', 'ALG3', 'ALG5', 'ALG6', 'ALG8', 'ALG9', 'ALK', 'ALKBH1', 'ALKBH2', 'ALKBH3', 'ALKBH4', 'ALKBH5', 'ALKBH6', 'ALKBH7', 'ALKBH8', 'ALLC', 'ALMS1', 'ALMS1P1', 'ALOX12', 'ALOX12B', 'ALOX12P2', 'ALOX15', 'ALOX15B', 'ALOX15P1', 'ALOX5', 'ALOX5AP', 'ALOXE3', 'ALPI', 'ALPK1', 'ALPK2', 'ALPK3', 'ALPL', 'ALPP', 'ALPPL2', 'ALS2', 'ALS2CL', 'ALS2CR11', 'ALS2CR12', 'ALX1', 'ALX3', 'ALX4', 'ALYREF', 'AMACR', 'AMBN', 'AMBP', 'AMBRA1', 'AMD1', 'AMDHD1', 'AMDHD2', 'AMELX', 'AMELY', 'AMER1', 'AMER2', 'AMER3', 'AMFR', 'AMH', 'AMHR2', 'AMIGO1', 'AMIGO2', 'AMIGO3', 'AMMECR1', 'AMMECR1L', 'AMN', 'AMN1', 'AMOT', 'AMOTL1', 'AMOTL2', 'AMPD1', 'AMPD2', 'AMPD3', 'AMPH', 'AMT', 'AMTN', 'AMY1A', 'AMY2A', 'AMY2B', 'AMZ1', 'AMZ2', 'AMZ2P1', 'ANAPC1', 'ANAPC10', 'ANAPC11', 'ANAPC13', 'ANAPC15', 'ANAPC16', 'ANAPC1P1', 'ANAPC2', 'ANAPC4', 'ANAPC5', 'ANAPC7', 'ANG', 'ANGEL1', 'ANGEL2', 'ANGPT1', 'ANGPT2', 'ANGPT4', 'ANGPTL1', 'ANGPTL2', 'ANGPTL3', 'ANGPTL4', 'ANGPTL5', 'ANGPTL6', 'ANGPTL7', 'ANGPTL8', 'ANHX', 'ANK1', 'ANK2', 'ANK3', 'ANKAR', 'ANKDD1A', 'ANKDD1B', 'ANKEF1', 'ANKFN1', 'ANKFY1', 'ANKH', 'ANKHD1', 'ANKHD1-EIF4EBP3', 'ANKIB1', 'ANKK1', 'ANKLE1', 'ANKLE2', 'ANKMY1', 'ANKMY2', 'ANKRA2', 'ANKRD1', 'ANKRD10', 'ANKRD11', 'ANKRD12', 'ANKRD13A', 'ANKRD13B', 'ANKRD13C', 'ANKRD13D', 'ANKRD16', 'ANKRD17', 'ANKRD18A', 'ANKRD18B', 'ANKRD18DP', 'ANKRD19P', 'ANKRD2', 'ANKRD20A11P', 'ANKRD20A12P', 'ANKRD20A19P', 'ANKRD20A3', 'ANKRD20A4', 'ANKRD20A5P', 'ANKRD20A8P', 'ANKRD20A9P', 'ANKRD22', 'ANKRD23', 'ANKRD24', 'ANKRD26', 'ANKRD26P1', 'ANKRD26P3', 'ANKRD27', 'ANKRD28', 'ANKRD29', 'ANKRD30A', 'ANKRD30B', 'ANKRD30BL', 'ANKRD30BP2', 'ANKRD30BP3', 'ANKRD31', 'ANKRD33', 'ANKRD33B', 'ANKRD34A', 'ANKRD34B', 'ANKRD34C', 'ANKRD35', 'ANKRD36', 'ANKRD36B', 'ANKRD36BP1']

matrix, list_of_pathways, dict_of_genes = search_pathways_4_list(list_of_genes)

def plot_bar_chart_for_pathways(matrix,list_of_pathways):
    
    y = [0]*len(matrix)
    
    for i,pathway in enumerate(matrix):
        y[i] = sum(matrix[i])
        
    plt.bar(np.arange(len(y)),y)
    plt.xticks(np.arange(len(y)),list_of_pathways)
    
plot_bar_chart_for_pathways(matrix,list_of_pathways)

print(np.array(matrix))
print(dict_of_genes)


 