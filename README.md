# Overview

This repository uses two sources of differentially expressed genes in the mouse hippocampus and frontal cortex to look for overlap with lists of genes known to be related to the following psychiatric disorders: autism, schizophrenia, bipolar disorder, and major depressive disorder. 

# DEG data sources

["A global transcriptional atlas of the effect of acute sleep deprivation in the mouse frontal cortex" ](https://pmc.ncbi.nlm.nih.gov/articles/PMC11402219/)

- DEG list in spreadsheet `sourcedata/mmc2.xlsx` with the headers: [Gene_Stable_ID	Gene_Stable_ID_Version	Gene_Name	Gene_Description	logFC	logCPM	LR	PValue	FDR]

- Copied to `rawdata/mouse_frontal_cortex_acute_sd_degs.xlsx

["Altered hippocampal transcriptome dynamics following sleep deprivation"](https://link.springer.com/article/10.1186/s13041-021-00835-1#Sec21)

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

## Bipolar disorder list

- Filename: NIHMS1687813-supplement-Supplementary_Tables.xlsx

- Gene list curated from Mullins, et al 2021:

https://pmc.ncbi.nlm.nih.gov/articles/PMC8192451/

- Looks at GWAS to find genes associated with bipolar disorder

- Sheet named "Table S4" contains list of 161 significant genes

## Major depressive disorder list

- Filename: 41588_2026_2638_MOESM4_ESM.xlsx

- Gene list curated from Zhang, et al 2026:

https://www.nature.com/articles/s41588-026-02638-3

- Uses a list of MDD risk genes to look at in vivo functions in mice

- Sheet named "Supplementary Table 1" contains list of 112 risk genes

## Future directions

- MDDOmics for MDD data

    - [csuligroup](https://www.csuligroup.com/CellCom/Home)

[The National Institute on Aging Genetics of Alzheimer's Disease Data Storage Site](https://www.niagads.org/)

