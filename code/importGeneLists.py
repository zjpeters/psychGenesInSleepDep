#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 10:28:06 2026

@author: zjpeters
"""
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

rawdata = os.path.join('/','home','zjpeters','Documents','psychGenesInSleepDep','rawdata')
geneListLocation = os.path.join(rawdata,'geneLists')
#%% load hippocampus data
"""
upregulated genes stored in Table S2, p < 0.1
downregulated genes stored in Table S5, p < 0.1
dysregulated genes stored in Table S8, p < 0.05
"""

hippDEGsUp = pd.read_excel(os.path.join(rawdata, 'mouse_hippocampus_acute_sd_degs.xlsx'), sheet_name='Table S2', skiprows=1)

hippDEGsDown = pd.read_excel(os.path.join(rawdata, 'mouse_hippocampus_acute_sd_degs.xlsx'), sheet_name='Table S5', skiprows=1)

hippDEGsDys = pd.read_excel(os.path.join(rawdata, 'mouse_hippocampus_acute_sd_degs.xlsx'), sheet_name='Table S8', skiprows=1)
# combine up and downregulated genes into one dataframe
hippDEGs = pd.concat([hippDEGsUp, hippDEGsDown])

#%% create function for searching list within other list

def compareGeneLists(referenceList, listToCheck):
    referenceList = referenceList.dropna()
    listToCheck = listToCheck.dropna()
    casefoldList = []
    for gene in referenceList:
        casefoldList.append(gene.casefold())

    # look for which DEGs are potentially present in the reference list
    foundGenes = []
    for gene in np.array(listToCheck):
        if gene.casefold() in np.array(casefoldList):
            foundGenes.append(gene)
    return foundGenes

#%% load frontal cortex data

sheetNames = pd.ExcelFile(os.path.join(rawdata, 'mouse_frontal_cortex_acute_sd_degs.xlsx'))
degDict = {}
for cellType in sheetNames.sheet_names:
    degDict[cellType] = pd.read_excel(os.path.join(rawdata, 'mouse_frontal_cortex_acute_sd_degs.xlsx'), sheet_name=cellType)
# remove unlabeled cells from data

#%% load ensembl genes and gene counts from original hippocampus data
bulkGeneList = pd.read_csv(os.path.join(rawdata, 'GSE166831_gene_names.csv'), index_col=0)
#%% create dictionary containing only brain diseases of interest
"""
identify genes present in xenium data that are also present in:
    "A comparison of anatomic and cellular transcriptome structures across 40 human brain diseases"
    https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002058#pbio.3002058.s002
"""

zeighamiGeneList = pd.read_excel(os.path.join(geneListLocation,"journal.pbio.3002058.s002.xlsx"), sheet_name='Associated genes')

disOfInterest = ['Autistic Disorder', 'Bipolar Disorder', 'Depressive Disorder', 'Schizophrenia', 'Sleep disorders']
zeighamiGenes = {}
for disease in disOfInterest:
    zeighamiGenes[disease] = zeighamiGeneList[disease].dropna()
    
zeighamiGenesCasefold = {}
for disease in zeighamiGenes.keys():
    tempGeneList = []
    for gene in np.array(zeighamiGenes[disease]):
        tempGeneList.append(gene.casefold())
    zeighamiGenesCasefold[disease] = tempGeneList

# look for which DEGs are potentially present in disease list
zeighamiDegDict = {}
for disease in zeighamiGenesCasefold.keys():
    zeighamiDegDict[disease] = []
    for deg in np.array(hippDEGs['Gene Name']):
        if deg.casefold() in np.array(zeighamiGenesCasefold[disease]):
            zeighamiDegDict[disease].append(deg)
            
#%% import and check SFARI gene list
sfariGeneList = pd.read_csv(os.path.join(geneListLocation, 'SFARI-Gene_genes_05-01-2026release_06-13-2026export.csv'))

sfariGenes = sfariGeneList['gene-symbol']
    
sfariGenesCasefold = []
for gene in sfariGenes:
    sfariGenesCasefold.append(gene.casefold())

# look for which DEGs are potentially present in disease list
sfariDegList = []

for deg in hippDEGs['Gene Name']:
    if deg.casefold() in np.array(sfariGenesCasefold):
        sfariDegList.append(deg)

#%% import and check schizophrenia gene list
schizGeneList = pd.read_csv(os.path.join(geneListLocation, 'INT-17_SCZ_High_Confidence_Gene_List.csv'))

schizGenes = schizGeneList['sczgenenames'].dropna()
    
schizGenesCasefold = []
for gene in schizGenes:
    schizGenesCasefold.append(gene.casefold())

# look for which DEGs are potentially present in disease list
schizDegList = []

for deg in hippDEGs['Gene Name']:
    if deg.casefold() in np.array(schizGenesCasefold):
        schizDegList.append(deg)
#%% check for overlap with overall list
schizGenesInList = []
for gene in np.array(list(bulkGeneList['Gene_Name'])):
    print(gene.casefold())
    if gene.casefold() in np.array(schizGenesCasefold):
        schizGenesInList.append(gene)
#%% import and check bipolar gene list
bpGeneList = pd.read_excel(os.path.join(geneListLocation, 'NIHMS1687813-supplement-Supplementary_Tables.xlsx'), sheet_name='Table S4', skiprows=1)

bpGenes = bpGeneList['Gene '].dropna()
    
bpCasefold = []
for gene in bpGenes:
    bpCasefold.append(gene.casefold())

# look for which DEGs are potentially present in disease list
bpDegList = []

for deg in allDegs:
    if deg.casefold() in np.array(bpCasefold):
        bpDegList.append(deg)

#%% import and check MDD gene list
mddGeneList = pd.read_excel(os.path.join(geneListLocation, '41588_2026_2638_MOESM4_ESM.xlsx'), sheet_name='Supplementary Table 1')

mddGenes = mddGeneList['Risk gene'].dropna()
    
mddCasefold = []
for gene in mddGenes:
    mddCasefold.append(gene.casefold())

# look for which DEGs are potentially present in disease list
mddDegList = []

for deg in allDegs:
    if deg.casefold() in np.array(mddCasefold):
        mddDegList.append(deg)

#%% write results to excel file

writer = pd.ExcelWriter(os.path.join(derivatives, 'Xenium_SD_DEGs_in_psych_risk_gene_lists.xlsx'))

degListDF = pd.DataFrame(allDegs, columns=['gene-name'])
degListDF.to_excel(writer, sheet_name='Xenium SD unique DEG list', index=False)
for disease in zeighamiDegDict.keys():
    if len(zeighamiDegDict[disease]) > 0:
        degListDF = pd.DataFrame(zeighamiDegDict[disease], columns=['gene-name'])
        degListDF.to_excel(writer, sheet_name=f'Zeighami {disease}', index=False)

degListDF = pd.DataFrame(sfariDegList, columns=['gene-name'])
degListDF.to_excel(writer, sheet_name='SFARI', index=False)

degListDF = pd.DataFrame(schizDegList, columns=['gene-name'])
degListDF.to_excel(writer, sheet_name='Schizophrenia', index=False)

degListDF = pd.DataFrame(bpDegList, columns=['gene-name'])
degListDF.to_excel(writer, sheet_name='Bipolar disorder', index=False)

degListDF = pd.DataFrame(mddDegList, columns=['gene-name'])
degListDF.to_excel(writer, sheet_name='Major depressive disorder', index=False)
writer.close()