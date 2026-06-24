#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 05:21:25 2026

@author: zjpeters
"""
import os
import pandas as pd
import numpy as np
import requests
import json
rawdata = os.path.join('/','home','zjpeters','Documents','psychGenesInSleepDep','rawdata')

hippBulkData = pd.read_csv(os.path.join(rawdata, 'GSE166831_hipp_counts_all_samples_NSD_and_SD.csv'))
#%% write function to convert Ensembl ID to gene ID
geneList = ["ENSMUSG00000082108","ENSMUSG00000022485","ENSMUSG00000050621"]
jsonGeneList = json.dumps({'api': 1, 'ids': geneList})
def ensemblIDToGeneName(ensemblID):
    urlBase = 'http://biotools.fr/mouse/ensembl_symbol_converter/?api=1&id='
    # ensemblID = 'ENSMUSG00000082108'
    getRequest = requests.get(f"{urlBase}{ensemblID}")
    jsonRequest = getRequest.json()
    return jsonRequest[ensemblID]

ensemblIDToGeneName('ENSMUSG00000082108')

#%% remove decimals from ensembl ids

convertedGeneList = []
for gene in hippBulkData['EnsemblID']:
    convertedGeneList.append(gene.split('.')[0])
#%% Use longer list to generate
"""
website can only handle ~290 genes at a time, so breaking it into smaller 
pieces is required to make it work
"""
urlBase = 'http://biotools.fr/mouse/ensembl_symbol_converter/'
nGenes = len(convertedGeneList)
nRepeats = int(np.ceil(nGenes / 275))
renamedGeneList = {}
listIdx = np.linspace(0, nGenes, nRepeats, dtype='int32')
startingIdx = 0
for i in listIdx:
    if i == 0:
        continue
    currGeneList = convertedGeneList[startingIdx:i]
    jsonData = {'api': 1, 'ids': json.dumps(currGeneList)}
    postRequest = requests.get(urlBase, jsonData)
    currDict = postRequest.json()
    renamedGeneList = renamedGeneList | currDict
    # print(startingIdx, i)
    startingIdx = i
    
renamedDF = pd.DataFrame.from_dict(renamedGeneList, orient="index", columns=['Gene_Name'])

renamedDF.to_csv(os.path.join(rawdata, 'GSE166831_gene_names.csv'))
