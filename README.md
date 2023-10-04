# snRNA-sequencing-data-analysis-pipeline
### This repository contains the standard workflow steps that can analyse single cell RNA sequencing data by clustering the genes, estimating the differentially expressed genes and performing functional enrichment analysis of each cluster.  


## Project: Single-cell profiling of breast cancer T cells
### Data: trm_counts.tsv.gz
### Data description: Cells from two samples were aggregated using UMI downsampling as implemented in Cellranger. The genes detected in less than three cells were filtered out. Multiplet captures were removed by excluding cells with total UMI count > 7500 or cells with detected number of genes > 2500. Mitochondrial filtering was also carried out in the dataset, despite code being present to do it as well. 
### Data source: https://singlecell.broadinstitute.org/single_cell/study/SCP2331/single-cell-profiling-of-breast-cancer-t-cells-reveals-a-tissue-resident-memory-subset-associated-with-improved-prognosis#study-download

## Workflow steps to analyze the single cell RNA sequencing data:
1. Open the downloaded single cell sequencing data using the read_tsv() or read_csv() function.  


