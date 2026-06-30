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
import scipy.stats as sp_stats

rawdata = os.path.join('/','home','zjpeters','Documents','psychGenesInSleepDep','rawdata')
derivatives = os.path.join('/','home','zjpeters','Documents','psychGenesInSleepDep', 'derivatives')
geneListLocation = os.path.join(rawdata,'geneLists')

#%% load hippocampus data
"""
upregulated genes stored in Table S2, p < 0.1
downregulated genes stored in Table S5, p < 0.1
dysregulated genes stored in Table S8, p < 0.05
"""

hippDEGsUp = pd.read_excel(os.path.join(rawdata, 'mouse_hippocampus_acute_sd_bulkRNASeq_degs.xlsx'), sheet_name='Table S2', skiprows=1)

hippDEGsDown = pd.read_excel(os.path.join(rawdata, 'mouse_hippocampus_acute_sd_bulkRNASeq_degs.xlsx'), sheet_name='Table S5', skiprows=1)

hippDEGsDys = pd.read_excel(os.path.join(rawdata, 'mouse_hippocampus_acute_sd_bulkRNASeq_degs.xlsx'), sheet_name='Table S8', skiprows=1)
# combine up and downregulated genes into one dataframe
hippDEGs = pd.concat([hippDEGsUp, hippDEGsDown])

# load ensembl gene names from original hippocampus data
bulkGeneList = pd.read_csv(os.path.join(rawdata, 'GSE166831_gene_names.csv'), index_col=0)
#%% load frontal cortex data

sheetNames = pd.ExcelFile(os.path.join(rawdata, 'mouse_frontal_cortex_acute_sd_snRNASeq_degs.xlsx'))
degDictFCsnRNASeq = {}
for cellType in sheetNames.sheet_names:
    degDictFCsnRNASeq[cellType] = pd.read_excel(os.path.join(rawdata, 'mouse_frontal_cortex_acute_sd_snRNASeq_degs.xlsx'), sheet_name=cellType)

cortDEGs = pd.read_excel(os.path.join(rawdata, 'mouse_frontal_cortex_acute_sd_bulkRNASeq_degs.xlsx'), sheet_name='DGE')

#%% create function for searching list within other list

def compareGeneLists(referenceList, listToCheck):
    """
    Takes to Pandas formatted lists of genes and looks for genes from the 
    listToCheck inside of the referenceList. The two lists must have the same
    sort of gene labeling, i.e. both Ensembl IDs, or both gene names, otherwise
    the search will not work properly. Can use code in 'convertEnsemblToGeneName.py'
    to convert from Ensembl to gene name.    
    
    Calculates the survival function of the hypergeometric distribution assuming
    population size (N) is number of genes in bulk, total number of successes 
    (K) is the number of genes from listToCheck in bulk list, quantity drawn
    per trial (n) is the number of genes from referenceList, and the number of
    observed successes (k) is the overlap between listToCheck and referenceList

    Parameters
    ----------
    referenceList : pandas series
        The list that will be searched for matching the listToCheck.
    listToCheck : pandas series
        Will look for genes from the listToCheck within the referenceList.

    Returns
    -------
    foundGenes : list
        List of genes from listToCheck that are present in referenceList.
    p_val : float
        Result of hypergeometric distribution survival function
    """
    # remove unlabeled cells from data
    referenceList = referenceList.dropna()
    listToCheck = np.array(listToCheck.dropna())
    # convert gene list to lowercase for search
    casefoldList = []
    for gene in referenceList:
        casefoldList.append(gene.casefold())
    casefoldList = np.array(casefoldList)
    # look for which DEGs are potentially present in the reference list
    foundGenes = []
    for gene in listToCheck:
        if gene.casefold() in casefoldList:
            foundGenes.append(gene)
    # calculate p-value for overlap 
    # import info from bulk list
    bulkList = bulkGeneList['Gene_Name'].dropna()
    bulkCasefold = []
    for gene in bulkList:
        bulkCasefold.append(gene.casefold())
    bulkCasefold = np.array(bulkCasefold)
    genesInBulkList = []
    for gene in listToCheck:
        if gene.casefold() in bulkCasefold:
            genesInBulkList.append(gene)
    # 'population size', i.e. total number of genes in bulk RNASeq
    N = len(bulkList)
    # sample size being 'drawn', i.e. number of DEGs 
    n = len(referenceList)
    # number of total 'successes' in list, i.e. genes associated with the given disorder
    K = len(genesInBulkList)
    # number of 'successes' drawn, i.e. the number of overlapping genes
    k = len(foundGenes)
    # print(N, n, K, k)
    # scipy uses different variable naming conventions
    p_val = sp_stats.hypergeom(M=N, n=K, N=n).sf(k-1)
    print(p_val)
    return foundGenes, p_val

#%% create dictionary of all disease gene lists

sfariGeneList = pd.read_csv(os.path.join(geneListLocation, 'SFARI-Gene_genes_05-01-2026release_06-13-2026export.csv'))
schizGeneList = pd.read_csv(os.path.join(geneListLocation, 'INT-17_SCZ_High_Confidence_Gene_List.csv'))
bpGeneList = pd.read_excel(os.path.join(geneListLocation, 'NIHMS1687813-supplement-Supplementary_Tables.xlsx'), sheet_name='Table S4', skiprows=1)
mddGeneList = pd.read_excel(os.path.join(geneListLocation, '41588_2026_2638_MOESM4_ESM.xlsx'), sheet_name='Supplementary Table 1')
adGeneList = pd.read_csv(os.path.join(geneListLocation, 'agora_ad_gene-list.csv'))
bp2GeneList = pd.read_excel(os.path.join(geneListLocation, '41586_2024_8468_MOESM4_ESM.xlsx'), sheet_name='S31', skiprows=16)

# use general names, append PMID if list comes from specific paper
diseaseGeneLists = dict.fromkeys(list(['SFARI', 'Schiz_psychencode', 'Bipolar_34002096', 'Bipolar_39843750', 'MDD_42271086', 'Alzheimers_agora']))
# create dictionary with key value for various things needed for analysis and saving
for i in diseaseGeneLists.keys():
    diseaseGeneLists[i] = dict.fromkeys(['Gene_name', 'Shared_genes_cort', 'Shared_genes_hipp', 'p-value_cort', 'p-value_hipp'])
diseaseGeneLists['SFARI']['Gene_name'] = sfariGeneList['gene-symbol']
diseaseGeneLists['Schiz_psychencode']['Gene_name'] = schizGeneList['sczgenenames']
diseaseGeneLists['Bipolar_34002096']['Gene_name'] = bpGeneList['Gene ']
diseaseGeneLists['Bipolar_39843750']['Gene_name'] = bp2GeneList['GENE']
diseaseGeneLists['MDD_42271086']['Gene_name'] = mddGeneList['Risk gene']
diseaseGeneLists['Alzheimers_agora']['Gene_name'] = adGeneList['Gene Symbol']
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
    foundGenes, p_val = compareGeneLists(hippDEGs['Gene Name'], zeighamiGenes[disease])
    zeighamiInHippDegs[disease] = foundGenes

writer = pd.ExcelWriter(os.path.join(derivatives, 'hippocampus_bulk', 'zeighami_hipp_bulkRNASeq_DEG_lists.xlsx'))

for disease in disOfInterest:
    zInHipDegDiseaseList = pd.DataFrame(zeighamiInHippDegs[disease], columns=['Gene Name'])
    zInHipDegDiseaseList.to_excel(writer, sheet_name=disease, index=False)
writer.close()

zeighamiInCortBulkDegs = dict.fromkeys(disOfInterest)
for disease in disOfInterest:
    foundGenes, p_val = compareGeneLists(cortDEGs['Gene_Name'], zeighamiGenes[disease])
    zeighamiInCortBulkDegs[disease] = foundGenes

writer = pd.ExcelWriter(os.path.join(derivatives, 'frontal_cortical_bulk', 'zeighami_cort_bulkRNASeq_DEG_lists.xlsx'))

for disease in disOfInterest:
    zInCortDegDiseaseList = pd.DataFrame(zeighamiInCortBulkDegs[disease], columns=['Gene Name'])
    zInCortDegDiseaseList.to_excel(writer, sheet_name=disease, index=False)
writer.close()

#%% loop over disease gene lists and perform gene list comparisons

for disease in diseaseGeneLists.keys():
    diseaseGeneLists[disease]['Shared_genes_hipp'], diseaseGeneLists[disease]['p-value_hipp'] = compareGeneLists(hippDEGs['Gene Name'], diseaseGeneLists[disease]['Gene_name'])
    diseaseGeneLists[disease]['Shared_genes_cort'], diseaseGeneLists[disease]['p-value_cort'] = compareGeneLists(cortDEGs['Gene_Name'], diseaseGeneLists[disease]['Gene_name'])
