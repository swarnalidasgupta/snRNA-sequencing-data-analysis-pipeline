# Assignment for cluster analysis and determination of differentially expressed genes
# Workflow for snRNA sequencing data analysis using Seurat

# Load all required libraries
library(Seurat)
library(hdf5r)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)

# Load the downloaded  file
datafile <- read_tsv("/Users/Swarnali_1/Downloads/trm_counts.tsv.gz")
str(datafile)

# Observe the structure of the dataset
str(datafile)

# Create a Seurat file
trm_counts.seurat.object <- CreateSeuratObject(counts = datafile, project = "trm")

#1. Quality control 
View(trm_counts.seurat.object@meta.data)
trm_counts.seurat.object[["percent_mt"]] <- PercentageFeatureSet(trm_counts.seurat.object, pattern = "^MT-")
View(trm_counts.seurat.object@meta.data)
# View the data in a plot to observe the distribution of the RNA transcripts and the genes
FeatureScatter(trm_counts.seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")

#2. Filtering the data
trm_counts.seurat.object <- subset(trm_counts.seurat.object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

#3. Normalising the data
trm_counts.seurat.object <- NormalizeData(trm_counts.seurat.object)
str(trm_counts.seurat.object)

# 4. Identifying the variable features
trm_counts.seurat.object <- FindVariableFeatures(trm_counts.seurat.object, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(trm_counts.seurat.object), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(trm_counts.seurat.object)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 5. Scaling
all.genes <- rownames(trm_counts.seurat.object)
trm_counts.seurat.object <- ScaleData(trm_counts.seurat.object, features = all.genes)

# 6. Linear dimensionality reduction (PCA)
trm_counts.seurat.object <- RunPCA(trm_counts.seurat.object, features = VariableFeatures(object = trm_counts.seurat.object))
# visualize PCA results
print(trm_counts.seurat.object[["pca"]], dims = 1:10, nfeatures = 10)
# DimHeatmap(gse_counts.seurat.object, dims = 1, cells = 500, balanced = TRUE)
# determine dimensionality of the data
ElbowPlot(trm_counts.seurat.object)

# 7. Clustering
trm_counts.seurat.object <- FindNeighbors(trm_counts.seurat.object, dims = 1:15)
trm_counts.seurat.object <- FindClusters(trm_counts.seurat.object, resolution = 0.5)
View(trm_counts.seurat.object@meta.data)

DimPlot(trm_counts.seurat.object, group.by = "RNA_snn_res.1", label = TRUE)

# 8. Non-linear dimensionality reduction 
trm_counts.seurat.object <- RunUMAP(trm_counts.seurat.object, dims = 1:10)
DimPlot(trm_counts.seurat.object, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
View(trm_counts.seurat.object@meta.data)

# 9. Finding the differentially expressed genes in the clusters
markers_data <- FindAllMarkers(trm_counts.seurat.object, logfc.threshold = 0.25, min.pct = 0.25)
View(markers_data)
#write.csv(markers_data, file = "markers_data.csv", row.names = FALSE)

# 10. Functional enrichment analysis of each cluster
# Extract the genes in terms of the clusters
unique_cluster_names <- unique(markers_data$cluster)
# Create a list to contain the enrichment results
cluster_enrichment_results <- list()
# Loop through the unique cluster names to extract out all the DEGs into a vector
for (i in unique_cluster_names){
  cluster_specific_genes <- markers_data$gene[markers_data$cluster == i] 
# Perform functional enrichment analysis for the clusters
    results <- enrichGO(
    gene = cluster_specific_genes,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05
  )
# Store the results in a list
  cluster_enrichment_results[[i]] <- results
}

combined_results <- rbind(cluster_enrichment_results)
output_file <- "cluster_enrichment_results.csv"

# Export the combined results to a CSV file
write.csv(combined_results, file = output_file, row.names = FALSE)
# Print a message indicating that the file has been saved
cat("Cluster enrichment results have been saved to", output_file, "\n")

