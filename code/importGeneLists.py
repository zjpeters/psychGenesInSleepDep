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

#%% import and check SFARI gene list
sfariGeneList = pd.read_csv(os.path.join(geneListLocation, 'SFARI-Gene_genes_05-01-2026release_06-13-2026export.csv'))

sfariInHippDegs, pValSfariHipp = compareGeneLists(hippDEGs['Gene Name'], sfariGeneList['gene-symbol'])
sfariHippDF = pd.DataFrame(sfariInHippDegs, columns=['Gene Name'])
sfariHippDF.to_excel(os.path.join(derivatives, 'hippocampus_bulk', 'sfari_hipp_bulkRNASeq_DEG_lists.xlsx'), index=False)

sfariInCortDegs, pValSfariCort = compareGeneLists(cortDEGs['Gene_Name'], sfariGeneList['gene-symbol'])
sfariCortDF = pd.DataFrame(sfariInCortDegs, columns=['Gene Name'])
sfariCortDF.to_excel(os.path.join(derivatives, 'frontal_cortical_bulk', 'sfari_cort_bulkRNASeq_DEG_lists.xlsx'), index=False)

# looks at snRNASeq
sfariInCortsnRNASeqDegs = dict.fromkeys(degDictFCsnRNASeq.keys())
for cortLayer in degDictFCsnRNASeq.keys():
    sfariInCortsnRNASeqDegs[cortLayer], p_val = compareGeneLists(degDictFCsnRNASeq[cortLayer]['Gene_Name'], sfariGeneList['gene-symbol'])
    
writer = pd.ExcelWriter(os.path.join(derivatives, 'frontal_cortical_snRNASeq', 'sfari_cort_snRNASeq_DEG_lists.xlsx'))

for cortLayer in degDictFCsnRNASeq.keys():
    sfariInCortDegList = pd.DataFrame(sfariInCortsnRNASeqDegs[cortLayer], columns=['Gene Name'])
    sfariInCortDegList.to_excel(writer, sheet_name=cortLayer, index=False)
writer.close()
#%% import and check schizophrenia gene list
schizGeneList = pd.read_csv(os.path.join(geneListLocation, 'INT-17_SCZ_High_Confidence_Gene_List.csv'))

schizGenesInHippDegs, pValSchizHipp = compareGeneLists(hippDEGs['Gene Name'], schizGeneList['sczgenenames'])
schizHippDF = pd.DataFrame(schizGenesInHippDegs, columns=['Gene Name'])
schizHippDF.to_excel(os.path.join(derivatives, 'hippocampus_bulk', 'schiz_hipp_bulkRNASeq_DEG_lists.xlsx'), index=False)

schizGenesInCortDegs, pValSchizCort = compareGeneLists(cortDEGs['Gene_Name'], schizGeneList['sczgenenames'])
schizCortDF = pd.DataFrame(schizGenesInCortDegs, columns=['Gene Name'])
schizCortDF.to_excel(os.path.join(derivatives, 'frontal_cortical_bulk', 'schiz_cort_bulkRNASeq_DEG_lists.xlsx'), index=False)

schizInCortsnRNASeqDegs = dict.fromkeys(degDictFCsnRNASeq.keys())
for cortLayer in degDictFCsnRNASeq.keys():
    schizInCortsnRNASeqDegs[cortLayer], p_val = compareGeneLists(degDictFCsnRNASeq[cortLayer]['Gene_Name'], schizGeneList['sczgenenames'])
    
writer = pd.ExcelWriter(os.path.join(derivatives, 'frontal_cortical_snRNASeq', 'schiz_cort_snRNASeq_DEG_lists.xlsx'))

for cortLayer in degDictFCsnRNASeq.keys():
    schizInCortDegList = pd.DataFrame(schizInCortsnRNASeqDegs[cortLayer], columns=['Gene Name'])
    schizInCortDegList.to_excel(writer, sheet_name=cortLayer, index=False)
writer.close()
#%% import and check bipolar gene list
bpGeneList = pd.read_excel(os.path.join(geneListLocation, 'NIHMS1687813-supplement-Supplementary_Tables.xlsx'), sheet_name='Table S4', skiprows=1)

bpGenesInHippDegs, pValBPHipp = compareGeneLists(hippDEGs['Gene Name'], bpGeneList['Gene '])
bpHippDF = pd.DataFrame(bpGenesInHippDegs, columns=['Gene Name'])
bpHippDF.to_excel(os.path.join(derivatives, 'hippocampus_bulk', 'bp_hipp_bulkRNASeq_DEG_lists.xlsx'), index=False)

bpGenesInCortDegs, pValBPCort = compareGeneLists(cortDEGs['Gene_Name'], bpGeneList['Gene '])
bpCortDF = pd.DataFrame(bpGenesInCortDegs, columns=['Gene Name'])
bpCortDF.to_excel(os.path.join(derivatives, 'frontal_cortical_bulk', 'bp_cort_bulkRNASeq_DEG_lists.xlsx'), index=False)

bpInCortsnRNASeqDegs = dict.fromkeys(degDictFCsnRNASeq.keys())
for cortLayer in degDictFCsnRNASeq.keys():
    bpInCortsnRNASeqDegs[cortLayer], p_val = compareGeneLists(degDictFCsnRNASeq[cortLayer]['Gene_Name'], bpGeneList['Gene '])
    
writer = pd.ExcelWriter(os.path.join(derivatives, 'frontal_cortical_snRNASeq', 'bp_cort_snRNASeq_DEG_lists.xlsx'))

for cortLayer in degDictFCsnRNASeq.keys():
    bpInCortDegList = pd.DataFrame(bpInCortsnRNASeqDegs[cortLayer], columns=['Gene Name'])
    bpInCortDegList.to_excel(writer, sheet_name=cortLayer, index=False)
writer.close()
#%% import and check MDD gene list
mddGeneList = pd.read_excel(os.path.join(geneListLocation, '41588_2026_2638_MOESM4_ESM.xlsx'), sheet_name='Supplementary Table 1')

mddGenesInHippDegs, pValMddHipp = compareGeneLists(hippDEGs['Gene Name'], mddGeneList['Risk gene'])
mddHippDF = pd.DataFrame(mddGenesInHippDegs, columns=['Gene Name'])
mddHippDF.to_excel(os.path.join(derivatives, 'hippocampus_bulk', 'mdd_hipp_bulkRNASeq_DEG_lists.xlsx'), index=False)

mddGenesInHippDegs, pValMddCort = compareGeneLists(cortDEGs['Gene_Name'], mddGeneList['Risk gene'])
mddCortDF = pd.DataFrame(mddGenesInHippDegs, columns=['Gene Name'])
mddCortDF.to_excel(os.path.join(derivatives, 'frontal_cortical_bulk', 'mdd_cort_bulkRNASeq_DEG_lists.xlsx'), index=False)

mddInCortsnRNASeqDegs = dict.fromkeys(degDictFCsnRNASeq.keys())
for cortLayer in degDictFCsnRNASeq.keys():
    mddInCortsnRNASeqDegs[cortLayer], p_val = compareGeneLists(degDictFCsnRNASeq[cortLayer]['Gene_Name'], mddGeneList['Risk gene'])
    
writer = pd.ExcelWriter(os.path.join(derivatives, 'frontal_cortical_snRNASeq', 'mdd_cort_snRNASeq_DEG_lists.xlsx'))

for cortLayer in degDictFCsnRNASeq.keys():
    mddInCortDegList = pd.DataFrame(mddInCortsnRNASeqDegs[cortLayer], columns=['Gene Name'])
    mddInCortDegList.to_excel(writer, sheet_name=cortLayer, index=False)
writer.close()

#%% import and check AD gene list
adGeneList = pd.read_csv(os.path.join(geneListLocation, 'agora_ad_gene-list.csv'))

adGenesInHippDegs, pValAdHipp = compareGeneLists(hippDEGs['Gene Name'], adGeneList['Gene Symbol'])
adHippDF = pd.DataFrame(adGenesInHippDegs, columns=['Gene Name'])
adHippDF.to_excel(os.path.join(derivatives, 'hippocampus_bulk', 'ad_hipp_bulkRNASeq_DEG_lists.xlsx'), index=False)

adGenesInHippDegs, pValAdCort = compareGeneLists(cortDEGs['Gene_Name'], adGeneList['Gene Symbol'])
adCortDF = pd.DataFrame(adGenesInHippDegs, columns=['Gene Name'])
adCortDF.to_excel(os.path.join(derivatives, 'frontal_cortical_bulk', 'ad_cort_bulkRNASeq_DEG_lists.xlsx'), index=False)

adInCortsnRNASeqDegs = dict.fromkeys(degDictFCsnRNASeq.keys())
for cortLayer in degDictFCsnRNASeq.keys():
    adInCortsnRNASeqDegs[cortLayer], p_val = compareGeneLists(degDictFCsnRNASeq[cortLayer]['Gene_Name'], adGeneList['Gene Symbol'])
    
writer = pd.ExcelWriter(os.path.join(derivatives, 'frontal_cortical_snRNASeq', 'ad_cort_snRNASeq_DEG_lists.xlsx'))

for cortLayer in degDictFCsnRNASeq.keys():
    adInCortDegList = pd.DataFrame(adInCortsnRNASeqDegs[cortLayer], columns=['Gene Name'])
    adInCortDegList.to_excel(writer, sheet_name=cortLayer, index=False)
writer.close()

#%% import and check second BD list gene list
bp2GeneList = pd.read_excel(os.path.join(geneListLocation, '41586_2024_8468_MOESM4_ESM.xlsx'), sheet_name='S31', skiprows=16)

bp2GenesInHippDegs, pValBPHipp = compareGeneLists(hippDEGs['Gene Name'], bp2GeneList['GENE'])
bp2HippDF = pd.DataFrame(bp2GenesInHippDegs, columns=['Gene Name'])
bp2HippDF.to_excel(os.path.join(derivatives, 'hippocampus_bulk', 'bp2_hipp_bulkRNASeq_DEG_lists.xlsx'), index=False)

bp2GenesInCortDegs, pValBPCort = compareGeneLists(cortDEGs['Gene_Name'], bp2GeneList['GENE'])
bp2CortDF = pd.DataFrame(bp2GenesInCortDegs, columns=['Gene Name'])
bp2CortDF.to_excel(os.path.join(derivatives, 'frontal_cortical_bulk', 'bp2_cort_bulkRNASeq_DEG_lists.xlsx'), index=False)

bp2InCortsnRNASeqDegs = dict.fromkeys(degDictFCsnRNASeq.keys())
for cortLayer in degDictFCsnRNASeq.keys():
    bp2InCortsnRNASeqDegs[cortLayer], p_val = compareGeneLists(degDictFCsnRNASeq[cortLayer]['Gene_Name'], bp2GeneList['GENE'])
    
writer = pd.ExcelWriter(os.path.join(derivatives, 'frontal_cortical_snRNASeq', 'bp2_cort_snRNASeq_DEG_lists.xlsx'))

for cortLayer in degDictFCsnRNASeq.keys():
    bp2InCortDegList = pd.DataFrame(bp2InCortsnRNASeqDegs[cortLayer], columns=['Gene Name'])
    bp2InCortDegList.to_excel(writer, sheet_name=cortLayer, index=False)
writer.close()  
#%% calculate hypergeometric distribution of acute SD genes in disease list
# sfariGenesInBulk = compareGeneLists(bulkGeneList['Gene_Name'], sfariGeneList['gene-symbol'])

# # 'population size', i.e. total number of genes in bulk RNASeq
# N = len(bulkGeneList['Gene_Name'].dropna())
# # sample size being 'drawn', i.e. number of DEGs 
# n = len(hippDEGs)
# # number of total 'successes' in list, i.e. genes associated with the given disorder
# K = len(sfariGenesInBulk)
# # number of 'successes' drawn, i.e. the number of overlapping genes
# k = len(sfariInHippDegs)

# # probability that the number of overlapping DEGs with the gene list is by chance
# p_sfari_hipp = sp_stats.hypergeom(M=51718, n=198, N=1146).sf(14-1)
# print(p_sfari_hipp)
# print(N, n, K, k)