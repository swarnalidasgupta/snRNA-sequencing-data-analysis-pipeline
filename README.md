# snRNA-sequencing-data-analysis-pipeline
### This repository contains the standard workflow steps that can analyse single cell RNA sequencing data by clustering the genes, estimating the differentially expressed genes and performing functional enrichment analysis of each cluster.  


## Project: Single-cell profiling of breast cancer T cells
### Data: trm_counts.tsv.gz
### Data description: Cells from two samples were aggregated using UMI downsampling as implemented in Cellranger. The genes detected in less than three cells were filtered out. Multiplet captures were removed by excluding cells with total UMI count > 7500 or cells with detected number of genes > 2500. Mitochondrial filtering was also carried out in the dataset, despite code being present to do it as well. 
### Data source: https://singlecell.broadinstitute.org/single_cell/study/SCP2331/single-cell-profiling-of-breast-cancer-t-cells-reveals-a-tissue-resident-memory-subset-associated-with-improved-prognosis#study-download

## Workflow steps to analyze the single cell RNA sequencing data:
1. Open the downloaded single cell sequencing data using the read_tsv() or read_csv() function.  Study the structure of the dataset using the str() function.
2. Create a Seurat object using the CreateSeuratObject() function.
3. Now, we need to do quality control and filtering. We need to filter out cells which have higher mitochondrial gene percentage, because this signifies higher metabolic activity, so higher cell death. We need to look at the number of features or number of genes in the cell. If there is a very high or low number of genes in the cell, these cells need to be filtered out. View the metadata in the Seurat file. Create a new column which calculates the mitochondrial percentage. The function PercentageFeatureSet() can be used to match the pattern of (“^MT-”) to this new column.
4. View the features/genes and RNA count in a scatterplot using the function FeatureScatter().
5. Now, we must filter the data depending on the nFeatures_RNA and the nCount_RNA. Filter the cells using subset() to contain cells with gene features greater than 200 but less than 2500.
6. Now we must normalize the data. For that, we use NormalizeData() function of Seurat. We need to normalize the data so that we get the observations in relative measures so that we can compare across different cells.
7. The next step is to identify the variable features, i.e. finding out the differentially expressed genes. We can use the FindVariableFeatures() function for this. Then we can view the top 10 differentially expressed genes using the VariableFeatures() function, which returns a list.
8. Then there is scaling. Scaling is important to remove noise due to variations. The noise needs to be removed so that clustering does not occur because of these variations. The function ScaleData() can be used for this.
9. Linear dimensionality reduction is the next step. Dimensionality refers to the features that describe each data point in a dataset. It is important to do principal component analysis to reduce the dimensionality of the data so that it can be plotted in a 2D or 3D scale. This needs to be done while preserving the maximum of the variation data. Reducing the dimensionality makes it easier to cluster the cells into groups. So, first calculate PCA by using the function RunPCA. The PCA needs to be run on the variable features. After that, we need to print out the first 10 PCA values  as well as the top 10 gene features that contribute to each principal component. Then, we visualize the dimensionality of the data using an elbow plot, mainly to decide how many dimensions we want to cluster in.  
10. Clustering: This is one of the main steps. The main aim is to cluster the cells into groups of similar gene expression profile. For this, we can use functions like FindNeighbours() and FindClusters().  The FindNeighbors function in Seurat finds nearest neighbors for each cell based on the first 15 principal components (dims = 1:15). It helps identify cells that are close to each other in the reduced-dimensional space. The FindClusters() function in Seurat performs cell clustering on the Seurat object by clustering groups together cells that have similar gene expression profiles, helping to identify distinct cell populations or cell types within the data. The resolution parameter specifies the granularity of clustering. Higher values lead to more fine-grained clusters, while lower values lead to larger clusters. The result of this operation is stored in the Seurat object, with each cell assigned to a specific cluster. The DimPlot() function can be used to visualize the clusters.
11. Non linear dimensionality reduction: Linear dimensionality reduction techniques like Principal Component Analysis (PCA) and linear t-SNE are powerful for capturing linear relationships in data. However, real-world data often contains nonlinear relationships and structures that LDR methods might not fully capture. NLDR methods UMAP (Uniform Manifold Approximation and Projection) are designed to uncover nonlinear structures in the data.
12. Finding the differentially expressed genes within the clusters: We can use the function FindAllMarkers(). These marker genes are genes that exhibit significant differences in expression between clusters and can provide insights into the characteristics and identities of the cell populations within your dataset.
13. The next step is functional enrichment analysis of each cluster. Functional enrichment analysis of each cluster based on differentially expressed genes (DEGs) is a crucial step in understanding the biological significance of the gene expression changes within each cluster. We can perform this analysis using tools and packages like GO (Gene Ontology) analysis, KEGG (Kyoto Encyclopedia of Genes and Genomes) pathway analysis, or other annotation databases.


