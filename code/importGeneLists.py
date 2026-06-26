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
derivatives = os.path.join('/','home','zjpeters','Documents','psychGenesInSleepDep', 'derivatives')
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

#%% load frontal cortex data

sheetNames = pd.ExcelFile(os.path.join(rawdata, 'mouse_frontal_cortex_acute_sd_degs.xlsx'))
degDict = {}
for cellType in sheetNames.sheet_names:
    degDict[cellType] = pd.read_excel(os.path.join(rawdata, 'mouse_frontal_cortex_acute_sd_degs.xlsx'), sheet_name=cellType)


#%% load ensembl genes and gene counts from original hippocampus data
bulkGeneList = pd.read_csv(os.path.join(rawdata, 'GSE166831_gene_names.csv'), index_col=0)
#%% create function for searching list within other list

def compareGeneLists(referenceList, listToCheck):
    """
    Takes to Pandas formatted lists of genes and looks for genes from the 
    listToCheck inside of the referenceList. The two lists must have the same
    sort of gene labeling, i.e. both Ensembl IDs, or both gene names, otherwise
    the search will not work properly. Can use code in 'convertEnsemblToGeneName.py'
    to convert from Ensembl to gene name.    

    Parameters
    ----------
    referenceList : pandas series
        The list that will be searched for matching the listToCheck.
    listToCheck : pandas series
        Will look for genes from the listToCheck within the referenceList.

    Returns
    -------
    list
        List of genes from listToCheck that are present in referenceList.
    """
    # remove unlabeled cells from data
    referenceList = referenceList.dropna()
    listToCheck = listToCheck.dropna()
    # convert gene list to lowercase for search
    casefoldList = []
    for gene in referenceList:
        casefoldList.append(gene.casefold())
    # look for which DEGs are potentially present in the reference list
    foundGenes = []
    for gene in np.array(listToCheck):
        if gene.casefold() in np.array(casefoldList):
            foundGenes.append(gene)
    return foundGenes

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
    
zeighamiInHippDegs = dict.fromkeys(disOfInterest)
for disease in disOfInterest:
    foundGenes = compareGeneLists(hippDEGs['Gene Name'], zeighamiGenes[disease])
    zeighamiInHippDegs[disease] = foundGenes

writer = pd.ExcelWriter(os.path.join(derivatives, 'zeighami_DEG_lists.xlsx'))

for disease in disOfInterest:
    zInHipDegDiseaseList = pd.DataFrame(zeighamiInHippDegs[disease])
    zInHipDegDiseaseList.to_excel(writer, sheet_name=disease, index=False)
writer.close()
#%% import and check SFARI gene list
sfariGeneList = pd.read_csv(os.path.join(geneListLocation, 'SFARI-Gene_genes_05-01-2026release_06-13-2026export.csv'))

sfariInHippDegs = compareGeneLists(hippDEGs['Gene Name'], sfariGeneList['gene-symbol'])
sfariDF = pd.DataFrame(sfariInHippDegs)
sfariDF.to_excel(os.path.join(derivatives, 'sfari_hipp_DEG_lists.xlsx'))

sfariInCortDegs = dict.fromkeys(degDict.keys())
for cortLayer in degDict.keys():
    sfariInCortDegs[cortLayer] = compareGeneLists(degDict[cortLayer]['Gene_Name'], sfariGeneList['gene-symbol'])
    
writer = pd.ExcelWriter(os.path.join(derivatives, 'sfari_cort_DEG_lists.xlsx'))

for cortLayer in degDict.keys():
    sfariInCortDegList = pd.DataFrame(sfariInCortDegs[cortLayer])
    sfariInCortDegList.to_excel(writer, sheet_name=cortLayer, index=False)
writer.close()
#%% import and check schizophrenia gene list
schizGeneList = pd.read_csv(os.path.join(geneListLocation, 'INT-17_SCZ_High_Confidence_Gene_List.csv'))

schizGenesInHippDegs = compareGeneLists(hippDEGs['Gene Name'], schizGeneList['sczgenenames'])
schizDF = pd.DataFrame(schizGenesInHippDegs)
schizDF.to_excel(os.path.join(derivatives, 'schiz_hipp_DEG_lists.xlsx'))

schizInCortDegs = dict.fromkeys(degDict.keys())
for cortLayer in degDict.keys():
    schizInCortDegs[cortLayer] = compareGeneLists(degDict[cortLayer]['Gene_Name'], schizGeneList['sczgenenames'])
    
writer = pd.ExcelWriter(os.path.join(derivatives, 'schiz_cort_DEG_lists.xlsx'))

for cortLayer in degDict.keys():
    schizInCortDegList = pd.DataFrame(schizInCortDegs[cortLayer])
    schizInCortDegList.to_excel(writer, sheet_name=cortLayer, index=False)
writer.close()
#%% import and check bipolar gene list
bpGeneList = pd.read_excel(os.path.join(geneListLocation, 'NIHMS1687813-supplement-Supplementary_Tables.xlsx'), sheet_name='Table S4', skiprows=1)

bpGenesInHippDegs = compareGeneLists(hippDEGs['Gene Name'], bpGeneList['Gene '])
bpDF = pd.DataFrame(bpGenesInHippDegs)
bpDF.to_excel(os.path.join(derivatives, 'bp_hipp_DEG_lists.xlsx'))

bpInCortDegs = dict.fromkeys(degDict.keys())
for cortLayer in degDict.keys():
    bpInCortDegs[cortLayer] = compareGeneLists(degDict[cortLayer]['Gene_Name'], bpGeneList['Gene '])
    
writer = pd.ExcelWriter(os.path.join(derivatives, 'bp_cort_DEG_lists.xlsx'))

for cortLayer in degDict.keys():
    bpInCortDegList = pd.DataFrame(bpInCortDegs[cortLayer])
    bpInCortDegList.to_excel(writer, sheet_name=cortLayer, index=False)
writer.close()
#%% import and check MDD gene list
mddGeneList = pd.read_excel(os.path.join(geneListLocation, '41588_2026_2638_MOESM4_ESM.xlsx'), sheet_name='Supplementary Table 1')

mddGenesInHippDegs = compareGeneLists(hippDEGs['Gene Name'], mddGeneList['Risk gene'])
mddDF = pd.DataFrame(mddGenesInHippDegs)
mddDF.to_excel(os.path.join(derivatives, 'mdd_hipp_DEG_lists.xlsx'))

mddInCortDegs = dict.fromkeys(degDict.keys())
for cortLayer in degDict.keys():
    mddInCortDegs[cortLayer] = compareGeneLists(degDict[cortLayer]['Gene_Name'], mddGeneList['Risk gene'])
    
writer = pd.ExcelWriter(os.path.join(derivatives, 'mdd_cort_DEG_lists.xlsx'))

for cortLayer in degDict.keys():
    mddInCortDegList = pd.DataFrame(mddInCortDegs[cortLayer])
    mddInCortDegList.to_excel(writer, sheet_name=cortLayer, index=False)
writer.close()
