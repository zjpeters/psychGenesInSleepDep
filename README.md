# Overview

This repository uses two sources of differentially expressed genes (DEGs) in the mouse hippocampus and frontal cortex to look for overlap with lists of genes known to be related to the following psychiatric disorders: autism, schizophrenia, bipolar disorder, and major depressive disorder. 

# Analysis performed

- A list of genes associated one of the disorders (e.g. bipolar disorder) is compared to a list of DEGs (e.g. acute SD hippocampal bulk RNASeq DEGs) to look for overlap between the two lists and a spreadsheet of the overlapping genes is generated and output into the `derivatives` folder

- The hypergeometric distribution is calculated using the following assumptions (variable letters used are from wikipedia for hypergeometric mean, not from scipy): 

    - $N$, population size, i.e. total number of genes in bulk analysis

    - $K$, number of total successes in dataset, i.e. the number of genes in disease list

    - $n$, quantity drawn with each trial, i.e. the number of genes in the DEG list

    - $k$, number of observed successes, i.e. overlap between disease and DEG lists

# DEG data sources

## ["A global transcriptional atlas of the effect of acute sleep deprivation in the mouse frontal cortex" ](https://pmc.ncbi.nlm.nih.gov/articles/PMC11402219/)

- DEG list for bulk RNA-seq: `sourcedata/ScienceDirect_files_29Jun2026_07-17-43.237/1-s2.0-S2589004224019771-mmc8.xlsx`

    - Headers: [Gene_Stable_ID	Gene_Stable_ID_Version	Gene_Name	Gene_Description	log10mean	log2FC	pvalue	qvalue]

    - Copied to `rawdata/mouse_frontal_cortex_acute_sd_bulkRNASeq_degs.xlsx`

- DEG list for snRNA-seq: `sourcedata/ScienceDirect_files_29Jun2026_07-17-43.237/1-s2.0-S2589004224019771-mmc2.xlsx`
    
    - Headers: [Gene_Stable_ID	Gene_Stable_ID_Version	Gene_Name	Gene_Description	logFC	logCPM	LR	PValue	FDR]

    - Copied to `rawdata/mouse_frontal_cortex_acute_sd_snRNASeq_degs.xlsx

- Simultaneous bulk and single nuclear RNA-seq

- Did not show major differences in gene expression in non-neuronal cell types

- Most cell types showed downregulation rather than upregulation

- Majority of DEGs are present in glutamatergic neurons, while many of those in GABAergic neurons are shared with glutamatergic

## ["Altered hippocampal transcriptome dynamics following sleep deprivation"](https://link.springer.com/article/10.1186/s13041-021-00835-1#Sec21)

- Genes significantly upregulated by sleep deprivation list in spreadsheet `sourcedata/13041_2021_835_MOESM2_ESM.xlsx`, sheet 'Table S2, with headers: [Ensembl ID	Gene Name	Gene Type	LogFC	LogCPM	F	P-value	FDR	Effect Size]

- Genes significantly downregulated by sleep deprivation list in spreadsheet `sourcedata/13041_2021_835_MOESM2_ESM.xlsx`, sheet 'Table S5, with headers: [Ensembl ID	Gene Name	Gene Type	LogFC	LogCPM	F	P-value	FDR	Effect Size]

# Psychiatric disorder gene list sources

## "A comparison of anatomic and cellular transcriptome structures across 40 human brain diseases"

[A comparison of anatomic and cellular transcriptome structures across 40 human brain diseases](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002058#pbio.3002058.s002)

- Filename: journal.pbio.3002058.s002.xlsx

- Gene list curated from Zeighami, et al 2023

- Includes 40 different brain diseases, gene lists are on "Associated genes" sheet

- 206 genes for autism, 384 genes for bipolar, 251 genes for depressive disorder, 733 genes for schizophrenia, 7 genes for sleep disorders

## SFARI autism gene list

- Filename: SFARI-Gene_genes_05-01-2026release_06-13-2026export.csv

https://gene.sfari.org/database/gene-scoring/

- Columns include: status, gene-symbol, gene-name, ensembl-id chromosome genetic-category

- 1277 genes

## Psychencode schizophrenia high confidence gene list

- Filename: INT-17_SCZ_High_Confidence_Gene_List.csv

http://resource.psychencode.org/

- Columns include: scgenenames, ensembl_names

- 321 genes

## Bipolar_34002096

- [Genome-wide association study of over 40,000 bipolar disorder cases provides new insights into the underlying biology](https://pmc.ncbi.nlm.nih.gov/articles/PMC8192451/)

Bipolar disorder list

- Filename: NIHMS1687813-supplement-Supplementary_Tables.xlsx

- Gene list curated from Mullins, et al 2021:

- Looks at GWAS to find genes associated with bipolar disorder

- Sheet named "Table S4" contains list of 161 significant genes

## Bipolar_39843750

- [Genomics yields biological and phenotypic insights into bipolar disorder](https://www.nature.com/articles/s41586-024-08468-9)

- Looked into gene list from supplementary table `rawdata/geneLists/41586_2024_8468_MOESM4_ESM.xlsx` sheet S31

- https://pgc.unc.edu/for-researchers/download-results/

## MDD_42271086, Major depressive disorder list

- [Linking GWAS risk genes to transcriptional features of major depressive disorder via in vivo Perturb-seq](https://www.nature.com/articles/s41588-026-02638-3)

- Filename: 41588_2026_2638_MOESM4_ESM.xlsx

- Gene list curated from Zhang, et al 2026:

- Uses a list of MDD risk genes to look at in vivo functions in mice

- Sheet named "Supplementary Table 1" contains list of 112 risk genes

## [AGORA Alzheimer's Gene List](https://agora.adknowledgeportal.org/genes/nominated-targets)

- Spreadsheet stored in `rawdata/geneLists/agora_ad_gene-list.csv`

- column headers: [Gene Symbol	Nominations	Year First Nominated	Nominating Teams	Cohort Study	Program	Input Data	Pharos Class]

- From agora website:

    "More than 900 targets have been contributed to Agora. These nominations have been provided by researchers from the National Institute on Aging's Accelerating Medicines Partnership in Alzheimer's Disease (AMP-AD) consortium and Target Enablement to Accelerate Therapy Development for Alzheimer's Disease (TREAT-AD) centers, as well as other research teams. The focus of Agora is to host targets that have been identified using computational analyses of high-dimensional genomic, proteomic and/or metabolomic data derived from human samples. These nominations have been released early in the research discovery cycle in order to facilitate evaluation by the broader research community."

# Future directions

- MDDOmics for MDD data

    - [csuligroup](https://www.csuligroup.com/CellCom/Home)

[The National Institute on Aging Genetics of Alzheimer's Disease Data Storage Site](https://www.niagads.org/)

- [Alzheimer's disease gene list](https://agora.adknowledgeportal.org/)

    
    