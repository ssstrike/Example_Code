##Shawn Striker::TGen::6/29/2021##
#Context Notes
#After performing superenhancer identification and gene annotation we have two output files per cell type (BIN67,COV434,SCCOHT1), one for each method of gene annotation we used (GeneMapper,ChipPeakanno). We need to pull in the annotated genes we want from each and filter using COSMIC T1 oncogenes. Finally perform some simple intersections between each of the three cell types.

import pandas as pd
import os

#Run with or without the oncogene filter
oncogeneFilter=True

#change to directory which contains subdirectories holding GeneMapper_output and ChipPeakanno_output. Also contains oncogene list.
os.chdir("C:/Users/sstriker/Documents/SCCOHT/Superenhancers")
#os.chdir("C:/Users/sstriker/Documents/SCCOHT/R/csawV2/significantSEs")

#get oncogene filter list, obtained from COSMIC tier 1 oncogene list.
oncogeneList = []
df2 = pd.read_csv('Census_all_1_5_2021.csv', sep = ',')
for line in df2.index:
    if 'oncogene' in df2.loc[line]['Role in Cancer']: #filter out possible TSGs
        oncogeneList.append(df2.loc[line][df2.columns[0]])

#list of subsirectories which contain mappper files.
mapper = ["GeneMapper_output", "ChipPeakanno_output"]

#list of cell types.
cell = ["BIN67", "COV434","SCCOHT1"]
#cell = ["59M", "CAOV3","COV362","ES2", "OAW28","OVCAR3","OVCAR4"]

#pick correct column to use from geneMapper output.
mapperHeader=['PROXIMAL_GENES']
#mapperHeader=['OVERLAP_GENES', 'PROXIMAL_GENES', 'CLOSEST_GENE','NUM_LOCI']

os.chdir("C:/Users/sstriker/Documents/SCCOHT/Superenhancers/ROSECREAM")
#dictonary holds all the genecalls between each caller.
genecalls={}
for mapperType in mapper:
    for cellType in cell:
        fileListtmp = os.listdir(mapperType)
        fileList = [x for x in fileListtmp if cellType in x]#filter for only one celltype
        for file in fileList:
            prevGeneList = [] #hold already called genes, removes duplicates
            comboName = mapperType[0]+"_"+file.split("_")[1]
            if comboName not in genecalls.keys(): #start building of dict keys
                genecalls[comboName] = []
            if mapperType == 'GeneMapper_output':
                df = pd.read_csv(mapperType+"/"+file, sep = '\t')
                for line in df.index:
                    for headerColumn in mapperHeader:
                        if df.loc[line][headerColumn] == df.loc[line][headerColumn]: #ignore nan
                            for genex in df.loc[line][headerColumn].split(','):
                                if genex in prevGeneList:
                                    continue
                                if genex != genex: #ignore nan
                                    continue
                                if oncogeneFilter:
                                    if genex in oncogeneList:
                                        if genex not in genecalls[comboName]:
                                            genecalls[comboName].append(genex)
                                        prevGeneList.append(genex)
                                else:
                                    if genex not in genecalls[comboName]:
                                        genecalls[comboName].append(genex)
                                    prevGeneList.append(genex)
            else:
                df = pd.read_csv(mapperType+"/"+file, sep = ',')
                for line in df.index:
                    gene = df.loc[line]["symbol"]
                    if gene in prevGeneList:
                        continue
                    if gene != gene:
                        continue
                    if oncogeneFilter:
                        if gene in oncogeneList:
                            if gene not in genecalls[comboName]:
                                genecalls[comboName].append(gene)
                            prevGeneList.append(gene)
                    else:
                        if gene not in genecalls[comboName]:
                            genecalls[comboName].append(gene)
                        prevGeneList.append(gene)

#pull out data into cell type dict.
#----Cell line----#
Cellgenes = {}
for j in cell:
    tmpList = []
    for i in genecalls:
        if j in i:
            for x in genecalls[i]:
                if x not in tmpList:
                    tmpList.append(x)
    Cellgenes[j] = tmpList


#output intersections
d = {'three1' : [Cellgenes['BIN67'],Cellgenes['COV434'],Cellgenes['SCCOHT1']], 'two1' : [Cellgenes['BIN67'],Cellgenes['COV434']], 'two2' : [Cellgenes['BIN67'],Cellgenes['SCCOHT1']], 'two3' : [Cellgenes['COV434'],Cellgenes['SCCOHT1']]}
int2 = []
for d1 in d.keys():
    result = set(d[d1][0]).intersection(*d[d1][1:])
    for y in result:    
        int2.append(y)
int2=list(set(int2))
with open('SE_GeneList_CellLines_OGfilter_Gthan2.txt', 'w') as outfile:
    outfile.write('\n'.join(int2))

resultAll = set(Cellgenes['BIN67']).intersection(Cellgenes['COV434'],Cellgenes['SCCOHT1'])
with open('SE_GeneList_CellLines_OGfilter_intAll3.txt', 'w') as outfile:
    outfile.write('\n'.join(list(resultAll)))
    
result1 = set(Cellgenes['BIN67']+Cellgenes['COV434']+Cellgenes['SCCOHT1'])
with open('SE_GeneList_CellLines_OGfilter_Gthan1.txt', 'w') as outfile:
    outfile.write('\n'.join(list(result1)))

#Final Context Notes
#There have been muliple iterations of this script being used with different inputs and generated outputs. I didn't want to touch it too much so you could get the idea of my scripting potential but I did remove some bits that would have made no sense unless you knew the history.

