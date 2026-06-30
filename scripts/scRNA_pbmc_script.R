#installed and running the necessary libraries
suppressPackageStartupMessages({
  library(dplyr)    
  library(spatstat.core)
  library(Seurat)
  library(patchwork)
  library(DoubletFinder)
  library(SingleR)
  library(enrichR)
  library(CellChat)
  library(SingleCellExperiment)
  library(SeuratWrappers)
  library(tidyverse)
  library(monocle3)
  library(celldex)
})

set.seed(42)

directory <- "~/Desktop/SC/scbi_ds1" #path to the file, change accordingly

#Information of the data sample
Data_Sample <- data.frame(file = c('GSM4138872_scRNA_BMMC_D1T1.rds', 'GSM4138873_scRNA_BMMC_D1T2.rds', 'GSM4138874_scRNA_CD34_D2T1.rds', 'GSM4138875_scRNA_CD34_D3T1.rds' ), 
                          Names = c('BMMC_D1T1', 'BMMC_D1T2', 'CD34_D2T1', 'CD34_D3T1'), 
                          Donor = c('D1', 'D1', 'D2', 'D3'), 
                          Replicate = c('T1', 'T2', 'T1', 'T1'), 
                          Sex = c('F', 'F', 'M', 'M'), 
                          Group = c('BMMC', 'BMMC', 'CD34', 'CD34'))


loadDataSet <- function(directory, Data_Sample, i){        #Loading the object 
  filename <-paste(directory, Data_Sample$file[i], sep="/")
  raw_counts <- readRDS(file = filename)
  Data_Seurat_Sample <- CreateSeuratObject(counts=raw_counts, project=Data_Sample$Names[i], assay = "RNA")   #then Converting the data sample into Seurat Object
  return(Data_Seurat_Sample)
}

list_of_samples <- list()  #Loading the object
for (i in 1: length(Data_Sample$file)){
  Data_Seurat_Sample <- loadDataSet(directory, Data_Sample, i);
  list_of_samples <- c(list_of_samples, Data_Seurat_Sample);
}


#AddMetaData, Name,Donor, Replicate etc.
addMetaData <- function(Data_Seurat_Sample, Data_Sample,i){
  Data_Seurat_Sample$orig.ident <- Data_Sample$Names[i]
  Data_Seurat_Sample$Donor <- Data_Sample$Donor[i]
  Data_Seurat_Sample$Replicate <- Data_Sample$Replicate[i]
  Data_Seurat_Sample$Sex <- Data_Sample$Sex[i]
  Data_Seurat_Sample$Group <- Data_Sample$Group[i]
  return(Data_Seurat_Sample)
}

for (i in 1:length(list_of_samples)) {       #applying metadata to all the samples
  list_of_samples[[i]] <- addMetaData(list_of_samples[[i]], Data_Sample, i)
  sample <-  list_of_samples[[i]]
  print(paste("Sample", i))
  print(paste("Number of cells:", ncol(sample)))     #checking number of cells
  print(paste("Number of genes:", nrow(sample)))     #checking number of genes
  print("View Metadata")
  print(head(sample@meta.data))
  View(sample@meta.data)
  readline(prompt = "To view next line press [ENTER]") #the gene and cell counts can be seen while this code runs
}

#######################################WEEK1##########################################################################


##Mitochondrial Counts
for (i in 1:length(list_of_samples)) {
  sample_name <- Data_Sample$Names[i]
  assign(sample_name, list_of_samples[[i]])
}
#adding mitochondrial percentage to each sample
BMMC_D1T1[['percent.mt']] <- PercentageFeatureSet(list_of_samples[[1]], pattern = "^MT-")
BMMC_D1T2[['percent.mt']] <- PercentageFeatureSet(list_of_samples[[2]], pattern = "^MT-")
CD34_D2T1[['percent.mt']] <- PercentageFeatureSet(list_of_samples[[3]], pattern = "^MT-")
CD34_D3T1[['percent.mt']] <- PercentageFeatureSet(list_of_samples[[4]], pattern = "^MT-")


##Violin Plot to visualise features, counts and mitochondrial counts
VP_BMMC_D1T1 <- VlnPlot(BMMC_D1T1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
print(VP_BMMC_D1T1)
ggsave("~/Desktop/SC/Figures_Assignment_1/VP_BMMC_D1T1.png", plot = VP_BMMC_D1T1, width = 10, height = 10, dpi = 300) #path change accordingly## all paths

VP_BMMC_D1T2 <- VlnPlot(BMMC_D1T2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
print(VP_BMMC_D1T2)
ggsave("~/Desktop/SC/Figures_Assignment_1/VP_BMMC_D1T2.png", plot = VP_BMMC_D1T2, width = 10, height = 10, dpi = 300)

VP_CD34_D2T1 <- VlnPlot(CD34_D2T1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
print(VP_CD34_D2T1)
ggsave("~/Desktop/SC/Figures_Assignment_1/VP_CD34_D2T1.png", plot = VP_CD34_D2T1, width = 10, height = 10, dpi = 300)

VP_CD34_D3T1 <- VlnPlot(CD34_D3T1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
print(VP_CD34_D3T1)
ggsave("~/Desktop/SC/Figures_Assignment_1/VP_CD34_D3T1.png", plot = VP_CD34_D3T1, width = 10, height = 10, dpi = 300)

#FeatureScatter
FS_BMMC_D1T1 <- FeatureScatter(BMMC_D1T1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(FS_BMMC_D1T1)
ggsave("~/Desktop/SC/Figures_Assignment_1/FS_BMMC_D1T1.png", plot = FS_BMMC_D1T1, width = 10, height = 10, dpi = 300)

FS_BMMC_D1T2 <- FeatureScatter(BMMC_D1T2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(FS_BMMC_D1T2)
ggsave("~/Desktop/SC/Figures_Assignment_1/FS_BMMC_D1T2.png", plot = FS_BMMC_D1T2, width = 10, height = 10, dpi = 300)

FS_CD34_D2T1 <- FeatureScatter(CD34_D2T1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(FS_CD34_D2T1)
ggsave("~/Desktop/SC/Figures_Assignment_1/FS_CD34_D2T1.png", plot = FS_CD34_D2T1, width = 10, height = 10, dpi = 300)

FS_CD34_D3T1 <- FeatureScatter(CD34_D3T1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(FS_CD34_D3T1)
ggsave("~/Desktop/SC/Figures_Assignment_1/FS_CD34_D3T1.png", plot = FS_CD34_D3T1, width = 10, height = 10, dpi = 300)

#########################################WEEK2##############################################################################
## preprocessing steps 
preprocess_sample <- function(sample) {
  sample %>%
    subset(subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA < 7500) %>% #subseting the data to filter cells based on the number of detected genes and total counts, to remove too few or too many genes
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%  #normalising gene expression levels to account for differences in sequencing
    FindVariableFeatures() %>% #identifying genes with high variability across cells
    ScaleData(features = VariableFeatures(.)) %>% #standarditing expression of the variable genes
    RunPCA() %>% #to reduce data dimension's
    RunUMAP(dims = 1:20) %>% #dimension reduction and visualisaton
    FindNeighbors(dims = 1:20) %>% #linking cells based on similarity, important to cluster them
    FindClusters(resolution = 0.5) #define cell clusters at a specific resolution
}

set.seed(42)  #setting this seed is enough but wanted to avoid risk of receiving any randome number generation
BMMC_D1T1_preprocessed <- preprocess_sample(BMMC_D1T1) # preprocessing the samples
set.seed(42)
BMMC_D1T2_preprocessed <- preprocess_sample(BMMC_D1T2)
set.seed(42)
CD34_D2T1_preprocessed <- preprocess_sample(CD34_D2T1)
set.seed(42)
CD34_D3T1_preprocessed <- preprocess_sample(CD34_D3T1)

##ElbowPlot to retain principal componenets
EP_BMMC_D1T1 <- ElbowPlot(BMMC_D1T1_preprocessed)
ggsave("~/Desktop/SC/Figures_Assignment_1/EP_BMMC_D1T1.png", plot = EP_BMMC_D1T1, width = 10, height = 10, dpi = 300)
EP_BMMC_D1T2 <- ElbowPlot(BMMC_D1T2_preprocessed)
ggsave("~/Desktop/SC/Figures_Assignment_1/EP_BMMC_D1T2.png", plot = EP_BMMC_D1T2, width = 10, height = 10, dpi = 300)
EP_CD34_D2T1 <- ElbowPlot(CD34_D2T1_preprocessed)
ggsave("~/Desktop/SC/Figures_Assignment_1/EP_CD34_D2T1.png", plot = EP_CD34_D2T1, width = 10, height = 10, dpi = 300)
EP_CD34_D3T1 <- ElbowPlot(CD34_D3T1_preprocessed)
ggsave("~/Desktop/SC/Figures_Assignment_1/EP_CD34_D3T1.png", plot = EP_CD34_D3T1, width = 10, height = 10, dpi = 300)

##Doublet finder steps
#BMMC1_D1T1
sweep.res.list_BMMC_D1T1 <- paramSweep(BMMC_D1T1_preprocessed, PCs = 1:20, sct = FALSE) #possible doublet classification threshold, idntifying optimal parameters for doublet detection
sweep.stats_BMMC_D1T1 <- summarizeSweep(sweep.res.list_BMMC_D1T1, GT = FALSE) #summarising the result from paramSweep
bcmvn_BMMC_D1T1 <- find.pK(sweep.stats_BMMC_D1T1) #to find the best pK value for doublet scoring
ggplot(bcmvn_BMMC_D1T1, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- 0.01 #fixing the threshold because as the optimal pk changes everytime the code runs, and so does the column name

annotations <- BMMC_D1T1_preprocessed@meta.data$seurat_clusters #extracting cluster annotations from Seurat Object's metadata
homotypic.prop <- modelHomotypic(annotations) #calculating the proportion of homotypic doublet within the extracted clusters
nExp_poi <- round(0.07*nrow(BMMC_D1T1_preprocessed@meta.data)) #calculating expectde number of doublets, here 7% of the total cells are the doublet rate, we should usually consider between 5-10% depending on our dataaset.
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) #adjusting the expected number of doublets, by subtracting the estimate of the homotypic doublets.

BMMC_D1T1_DF <- doubletFinder(BMMC_D1T1_preprocessed, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE) #finding doublets
DP_BMMC_D1T1_DF <- DimPlot(BMMC_D1T1_DF, reduction = 'umap', group.by = "DF.classifications_0.25_0.01_364" )

table(BMMC_D1T1_DF@meta.data$DF.classifications_0.25_0.01_364)
BMMC_D1T1_DF_singlet <- subset(BMMC_D1T1_DF, subset = DF.classifications_0.25_0.01_364 == "Singlet") #Capturing the singlets

# Plotting UMAP for singlets only
DP_BMMC_D1T1_DF_singlet <- DimPlot(BMMC_D1T1_DF_singlet, reduction = 'umap')
print(DP_BMMC_D1T1_DF_singlet)
ggsave("~/Desktop/SC/Figures_Assignment_1/DP_BMMC_D1T1_DF_singlet.png", plot = DP_BMMC_D1T1_DF_singlet, width = 10, height = 10, dpi = 300)


#BMMC_D1T2
sweep.res.list_BMMC_D1T2 <- paramSweep(BMMC_D1T2_preprocessed, PCs = 1:20, sct = FALSE)
sweep.stats_BMMC_D1T2 <- summarizeSweep(sweep.res.list_BMMC_D1T2, GT = FALSE)
bcmvn_BMMC_D1T2 <- find.pK(sweep.stats_BMMC_D1T2)
ggplot(bcmvn_BMMC_D1T2, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- 0.01
annotations <- BMMC_D1T2_preprocessed@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.07*nrow(BMMC_D1T2_preprocessed@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

BMMC_D1T2_DF <- doubletFinder(BMMC_D1T2_preprocessed, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
DimPlot(BMMC_D1T2_DF, reduction = 'umap', group.by = "DF.classifications_0.25_0.01_376" )

table(BMMC_D1T2_DF@meta.data$DF.classifications_0.25_0.01_376)
BMMC_D1T2_DF_singlet <- subset(BMMC_D1T2_DF, subset = DF.classifications_0.25_0.01_376 == "Singlet")

DP_BMMC_D1T2_DF_singlet <- DimPlot(BMMC_D1T2_DF_singlet, reduction = 'umap')
print(DP_BMMC_D1T2_DF_singlet)
ggsave("~/Desktop/SC/Figures_Assignment_1/DP_BMMC_D1T2_DF_singlet.png", plot = DP_BMMC_D1T2_DF_singlet, width = 10, height = 10, dpi = 300)


#CD34_D2T1
sweep.res.list_CD34_D2T1 <- paramSweep(CD34_D2T1_preprocessed, PCs = 1:20, sct = FALSE)
sweep.stats_CD34_D2T1 <- summarizeSweep(sweep.res.list_CD34_D2T1, GT = FALSE)
bcmvn_CD34_D2T1 <- find.pK(sweep.stats_CD34_D2T1)
ggplot(bcmvn_CD34_D2T1, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- 0.01
annotations <- CD34_D2T1_preprocessed@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.07*nrow(CD34_D2T1_preprocessed@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

CD34_D2T1_DF <- doubletFinder(CD34_D2T1_preprocessed, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
DimPlot(CD34_D2T1_DF, reduction = 'umap', group.by = "DF.classifications_0.25_0.01_123")

table(CD34_D2T1_DF@meta.data$DF.classifications_0.25_0.01_123)
CD34_D2T1_DF_singlet <- subset(CD34_D2T1_DF, subset = DF.classifications_0.25_0.01_123 == "Singlet")

DP_CD34_D2T1_DF_singlet <- DimPlot(CD34_D2T1_DF_singlet, reduction = 'umap')
print(DP_CD34_D2T1_DF_singlet)
ggsave("~/Desktop/SC/Figures_Assignment_1/DP_CD34_D2T1_DF_singlet.png", plot = DP_CD34_D2T1_DF_singlet, width = 10, height = 10, dpi = 300)


#CD34_D3T1
sweep.res.list_CD34_D3T1 <- paramSweep(CD34_D3T1_preprocessed, PCs = 1:20, sct = FALSE)
sweep.stats_CD34_D3T1 <- summarizeSweep(sweep.res.list_CD34_D3T1, GT = FALSE)
bcmvn_CD34_D3T1 <- find.pK(sweep.stats_CD34_D3T1)
ggplot(bcmvn_CD34_D3T1, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()


pK <- 0.01
annotations <- CD34_D3T1_preprocessed@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.07*nrow(CD34_D3T1_preprocessed@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

CD34_D3T1_DF <- doubletFinder(CD34_D3T1_preprocessed, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
DimPlot(CD34_D3T1_DF, reduction = 'umap', group.by =  "DF.classifications_0.25_0.01_349")

table(CD34_D3T1_DF@meta.data$DF.classifications_0.25_0.01_349)
CD34_D3T1_DF_singlet <- subset(CD34_D3T1_DF, subset =  DF.classifications_0.25_0.01_349 == "Singlet")

DP_CD34_D3T1_DF_singlet <- DimPlot(CD34_D3T1_DF_singlet, reduction = 'umap')
print(DP_CD34_D3T1_DF_singlet)
ggsave("~/Desktop/SC/Figures_Assignment_1/DP_CD34_D3T1_DF_singlet.png", plot = DP_CD34_D3T1_DF_singlet, width = 10, height = 10, dpi = 300)



##Merging without batch correction
merge_without_batch_correction <- function(samples_list) {
  merged <- Reduce(function(x, y) merge(x, y), samples_list) %>% #repeating all the pre-processing steps before merging the samples, one for without batch correction and the other for the batch correction
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:20) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.5)
  
  return(merged)
}


set.seed(42) 
merged_without_batch <- merge_without_batch_correction(
  list(BMMC_D1T1_DF_singlet, BMMC_D1T2_DF_singlet, CD34_D2T1_DF_singlet, CD34_D3T1_DF_singlet)
)

# Visualisation 
DP_merge_without_batch <- DimPlot(merged_without_batch, reduction = 'umap', group.by = "orig.ident") +
  ggtitle("Merged Without Batch Correction")
print(DP_merge_without_batch)
ggsave("~/Desktop/SC/Figures_Assignment_1/DP_merge_without_batch.png", plot = DP_merge_without_batch, width = 10, height = 10, dpi = 300)

##Merging with batch correction
merge_with_batch_correction <- function(samples_list) {
  anchors <- FindIntegrationAnchors(object.list = samples_list, dims = 1:20) #identifying integration anchors between differentsamples to align with the cell population across batches
  merged <- IntegrateData(anchorset = anchors, dims = 1:20) %>% #integrating data using identified anchors and then reducing dimesnions for a single analysis (together)
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:20) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.5)
  
  DefaultAssay(merged) <- "integrated"
  
  return(merged)
}

set.seed(42)
merged_with_batch <- merge_with_batch_correction(
  list(BMMC_D1T1_DF_singlet, BMMC_D1T2_DF_singlet, CD34_D2T1_DF_singlet, CD34_D3T1_DF_singlet)
)

# Visualisation
#Samples
DP_merged_with_batch <- DimPlot(merged_with_batch, reduction = 'umap', group.by = "orig.ident") +
  ggtitle("Merged With Batch Correction")
print(DP_merged_with_batch)
ggsave("~/Desktop/SC/Figures_Assignment_1/DP_merged_with_batch.png", plot = DP_merged_with_batch, width = 10, height = 10, dpi = 300)

#Clusters
DP_merged_with_batch_clusters <- DimPlot(merged_with_batch, reduction = 'umap', group.by = "seurat_clusters") + 
  ggtitle("UMAP with clusters")
print(DP_merged_with_batch_clusters)
ggsave("~/Desktop/SC/Figures_Assignment_1/DP_merged_with_batch_clusters.png", plot = DP_merged_with_batch_clusters, width = 10, height = 10, dpi = 300)




##############################################WEEK3#####################################################################################################################

DefaultAssay(merged_with_batch) <- "RNA" #changing Assay
merged_with_batch[["RNA"]] <- JoinLayers(merged_with_batch[["RNA"]]) #to combine multiple data,counts etc to a single data, counts, etc.
RNA_data <- merged_with_batch[["RNA"]]@layers[["data"]]
RNA_matrix <- as.matrix(RNA_data)
rownames(RNA_matrix) <- rownames(merged_with_batch[["RNA"]])

dim(RNA_matrix)
hpca_ref <- HumanPrimaryCellAtlasData()  #reference for cell annotation
singleR_results <- SingleR(test = RNA_matrix, ref = hpca_ref, labels = hpca_ref$label.main) #cell type identification using tool SingleR by comparing the gene expression matrix 'RNA_matrix'  
merged_with_batch <- AddMetaData(merged_with_batch, metadata = singleR_results$labels, col.name = "SingleR_Labels")

#BROAD- general overview of the cell types in the dataset#main labels
singleR_broad <- SingleR( test = RNA_matrix, ref = hpca_ref, labels = hpca_ref$label.main)
merged_with_batch <- AddMetaData(merged_with_batch, metadata = singleR_broad$labels, col.name = "Broad_Labels")
Heatmap_Broad <- plotScoreHeatmap(singleR_broad, main = "Broad Cell Type Heatmap")
print(Heatmap_Broad)
ggsave("~/Desktop/SC/Figures_Assignment_1/Heatmap_Broad.png", plot = Heatmap_Broad, width = 20, height = 10, dpi = 300)

#FINE- detailed and specific cell type classification
singleR_fine <- SingleR( test = RNA_matrix, ref = hpca_ref, labels = hpca_ref$label.fine)
merged_with_batch <- AddMetaData(merged_with_batch, metadata = singleR_fine$labels, col.name = "Broad_Fine")
Heatmap_Fine <- plotScoreHeatmap(singleR_fine, main = "Fine Cell Type Heatmap")
print(Heatmap_Fine)
ggsave("~/Desktop/SC/Figures_Assignment_1/Heatmap_Fine.png", plot = Heatmap_Fine, width = 20, height = 10, dpi = 300)

#COMBINED- broad and fine, enriching cell type metdata
combined_labels <- paste(singleR_broad$labels, singleR_fine$labels)
merged_with_batch <- AddMetaData(merged_with_batch, metadata = combined_labels, col.name = "Combined_Labels")
singleR_broad_scores <- singleR_broad$scores
singleR_fine_scores <- singleR_fine$scores #heatmap, categorising cells, visualise the confidence and clarity
combined_scores <- cbind(singleR_broad_scores, singleR_fine_scores)
combined_singleR <- singleR_broad
combined_singleR$scores <- combined_scores
Heatmap_combined <- plotScoreHeatmap(combined_singleR, main = "Combined Cell Type Heatmap")
print(Heatmap_combined)
ggsave("~/Desktop/SC/Figures_Assignment_1/Heatmap_combined.png", plot = Heatmap_combined, width = 20, height = 10, dpi = 300)

#Cell_types
DP_cell_types <- DimPlot(merged_with_batch, reduction = "umap", group.by = "SingleR_Labels", label = TRUE) +
  ggtitle("UMAP with Cluster Cell Type Annotations")
print(DP_cell_types)
ggsave("~/Desktop/SC/Figures_Assignment_1/DP_cell_types.png", plot = DP_cell_types, width = 20, height = 10, dpi = 300)

DefaultAssay(merged_with_batch) <- "integrated" #changing to integrated for clustering
merged_with_batch <- FindClusters(merged_with_batch, resolution = 0.4)
DP_Cluster <- DimPlot(merged_with_batch, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Clusters plot")
print(DP_Cluster)
ggsave("~/Desktop/SC/Figures_Assignment_1/DP_Cluster.png", plot = DP_Cluster, width = 20, height = 10, dpi = 300)

DefaultAssay(merged_with_batch) <- "RNA" #changing it back to integrated in order to perform FindAllMarkers
DE_Analysis_Markers <- FindAllMarkers(object = merged_with_batch, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#markers from the table
marker_genes <- c("CD34", "CD38", "Sca1", "Kit") 
marker_genes_1 <- c("CD52", "CSF3R", "ca1", "Kit") 
marker_genes_2 <- c("Flk2", "IL7R", "ELANE", "IL3")
marker_genes_3 <- c("GM-CSF", "M-CSF", "CD19", "CD20") 
marker_genes_4 <- c("SDC1", "IGHA1", "IGLC1", "MZB1") 
marker_genes_5 <- c("JCHAIN", "CD3D", "CD3E", "CD8A", "CD8B") 
marker_genes_6 <- c("CD4", "FCGR3A", "NCAM1", "NKG7", "KLRB1") 
marker_genes_7 <- c("GATA1", "HBB", "HBA1", "HBA2", "IRF8", "IRF4", "IRF7") 
marker_genes_8 <- c("CD1C", "CD207", "ITGAM", "NOTCH2", "SIRPA") 
marker_genes_9 <- c("CD14", "CCL3", "CCL4", "IL1B", "FCGR3A") 
marker_genes_10 <- c("CD68", "S100A12", "GATA2")
#violin plot
VlnPlot(object = merged_with_batch, features = marker_genes, group.by = "seurat_clusters")
VlnPlot(object = merged_with_batch, features = marker_genes_1, group.by = "seurat_clusters")
VlnPlot(object = merged_with_batch, features = marker_genes_2, group.by = "seurat_clusters")
VlnPlot(object = merged_with_batch, features = marker_genes_3, group.by = "seurat_clusters")
VlnPlot(object = merged_with_batch, features = marker_genes_4, group.by = "seurat_clusters")
VlnPlot(object = merged_with_batch, features = marker_genes_5, group.by = "seurat_clusters")
VlnPlot(object = merged_with_batch, features = marker_genes_6, group.by = "seurat_clusters")
VlnPlot(object = merged_with_batch, features = marker_genes_7, group.by = "seurat_clusters")
VlnPlot(object = merged_with_batch, features = marker_genes_8, group.by = "seurat_clusters")
VlnPlot(object = merged_with_batch, features = marker_genes_9, group.by = "seurat_clusters")
VlnPlot(object = merged_with_batch, features = marker_genes_10, group.by = "seurat_clusters")
#Identify Cell types and Assign Cluster Names.
manual_annotations <- c(           #the question was about using abbreviations in brackets, so i did this. not sure if this is what it meant.
  "0" = "HSC(HSC)",
  "1" = "LMPP(LMPP)",
  "2" = "CLP(CLP)",
  "3" = "GMP(GMP)",
  "4" = "CMP(CMP)",
  "5" = "B Cells(B)",
  "6" = "Pre B-cell Progenitors(PreB)",
  "7" = "Plasma Cells(Plasma)",
  "8" = "CD8+ T Cells(CD8)",
  "9" = "CD4+ T Cells(CD4)",
  "10" = "NK Cells(NK)",
  "11" = "Erythrocytes(Ery)",
  "12" = "cDC(cDC)",
  "13" = "pDC(PDC)",
  "14" = "CD14+ Monocytes(Mono CD14)",
  "15" = "CD16+ Monocytes (Mono CD16)",
  "16" = "Basophils(Bas)"
)

#UMAP
names(manual_annotations) <- levels(merged_with_batch)
merged_with_batch <- RenameIdents(merged_with_batch, manual_annotations) #Renamed
DP_manual <- DimPlot(merged_with_batch, reduction = "umap", label = TRUE, pt.size = 0.5)
print(DP_manual)
ggsave("~/Desktop/SC/Figures_Assignment_1/DP_manual.png", plot = DP_manual, width = 20, height = 10, dpi = 300)

#Randomly selecting genes from Table 2 for visualisation
Violin_Plot_3 <- VlnPlot(merged_with_batch, features = c("CD3D", "GATA2", "CD14"), group.by = "SingleR_Labels") +
  theme(
    plot.title = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  )
print(Violin_Plot_3)
ggsave("~/Desktop/SC/Figures_Assignment_1/Violin_Plot_3.png", plot = Violin_Plot_3, width = 20, height = 10, dpi = 300)

#Visualisation using Feature plot
FeaturePlot_3 <- FeaturePlot(merged_with_batch, features = c("CD3D", "GATA2", "CD14"))
print(FeaturePlot_3)
ggsave("~/Desktop/SC/Figures_Assignment_1/FeaturePlot_3.png", plot = FeaturePlot_3, width = 10, height = 10, dpi = 300)

#Barplot
cell_type_counts <- data.frame(Project_data = merged_with_batch$orig.ident, cell_type = Idents(merged_with_batch)) #data frame containing original identifiers and cell types from merged datset

cell_type_proportions <- cell_type_counts %>% #calculating counts and proportions of each celltype conditon
  group_by(Project_data, cell_type) %>%
  summarise(count = n()) %>% #summarise to count the number of cells in each sample
  mutate(proportion = count / sum(count)) %>% #proportion of each cell in each condition
  ungroup()

Bar_Plot <- ggplot(cell_type_proportions, aes(x = Project_data, y = proportion, fill = cell_type)) +
         geom_bar(stat = "identity", position = "fill") +
         labs(
           title = "Cell-type Proportions per Sample",
           x = "Project_data",
           y = "Proportion"
         ) + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
print(Bar_Plot)
ggsave("~/Desktop/SC/Figures_Assignment_1/Bar_Plot.png", plot = Bar_Plot, width = 10, height = 10, dpi = 300)

#Differentially expresssed genes
##DEG_B_T_Cells
b_cells <- subset(merged_with_batch, idents = "B Cells(B)")
cd4_t_cells <- subset(merged_with_batch, idents = "CD4+ T Cells(CD4)") 
cd8_t_cells <- subset(merged_with_batch, idents = "CD8+ T Cells(CD8)")
t_cells <- subset(merged_with_batch, idents = c("CD4+ T Cells(CD4)", "CD8+ T Cells(CD8)")) #merging the two subsets of t cells into 1

DEG_B_T_Cells <- FindMarkers(merged_with_batch, ident.1 = "B Cells(B)", ident.2 = c("CD4+ T Cells(CD4)", "CD8+ T Cells(CD8)"), min.pct = 0.25, logfc.threshold = 0.25) #performing differential gene analysis between two cell types
DEG_B_T_Cells$gene <- rownames(DEG_B_T_Cells)
DEG_B_T_Cells$logp <- -log10(ifelse(DEG_B_T_Cells$p_val_adj == 0, 1e-300, DEG_B_T_Cells$p_val_adj))  #using "1e-300" to avoid getting p value 0 and so that the plot is visible

Volcano_Plot_B_T_Cells <- ggplot(DEG_B_T_Cells, aes(x = avg_log2FC, y = logp)) +
  geom_point(aes(color = avg_log2FC), alpha = 0.7) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Log2 Fold Change") +
  labs(
    title = "Volcano Plots of B Cells and T Cells",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
)
print(Volcano_Plot_B_T_Cells)
ggsave("~/Desktop/SC/Figures_Assignment_1/Volcano_Plot_B_T_Cells.png", plot = Volcano_Plot_B_T_Cells, width = 10, height = 10, dpi = 300)

##DEG_Monocytes_T_Cells
monocytes <- subset(merged_with_batch, idents = c("CD14+ Monocytes(Mono CD14)", "CD16+ Monocytes (Mono CD16)"))
DEG_Monocytes_T_Cells <- FindMarkers(merged_with_batch, ident.1 = c("CD4+ T Cells(CD4)", "CD8+ T Cells(CD8)"), ident.2 = c("CD14+ Monocytes(Mono CD14)", "CD16+ Monocytes (Mono CD16)"), min.pct = 0.25, logfc.threshold = 0.25
)
DEG_Monocytes_T_Cells$gene <- rownames(DEG_Monocytes_T_Cells)
DEG_Monocytes_T_Cells$logp <- -log10(ifelse(DEG_Monocytes_T_Cells$p_val_adj == 0, 1e-300, DEG_Monocytes_T_Cells$p_val_adj))

Volcano_Plot_Monocyte_T_Cells <- ggplot(DEG_Monocytes_T_Cells, aes(x = avg_log2FC, y = logp)) +
  geom_point(aes(color = avg_log2FC), alpha = 0.7) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Log2 Fold Change") +
  labs(
    title = "Volcano Plots of T Cells and Monocytes",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
)
print(Volcano_Plot_Monocyte_T_Cells)
ggsave("~/Desktop/SC/Figures_Assignment_1/Volcano_Plot_Monocyte_T_Cells.png", plot = Volcano_Plot_Monocyte_T_Cells, width = 10, height = 10, dpi = 300)


#top5_B_T_Cells#top5_#top5_Monocytes_T_Cells
B_T_Cells5 <- DEG_B_T_Cells %>%
  arrange(p_val_adj) %>% 
  head(5) #find top 5 genes statistical analysis

Monocytes_T_Cells5 <- DEG_Monocytes_T_Cells %>%
  arrange(p_val_adj) %>% 
  head(5)

B_T_Cells5$comparison <- "B cells vs T cells"
Monocytes_T_Cells5$comparison <- "Monocytes vs T cells"

top_genes <- rbind(B_T_Cells5, Monocytes_T_Cells5)

top_genes <- top_genes %>%
  mutate(logp = ifelse(p_val_adj == 0, -log10(1e-300), -log10(p_val_adj)))  # Set logp to 300 if p_val_adj is too small

Data_Plot <- top_genes %>%
  select(gene, comparison, avg_log2FC, logp)


DEG_Comparison <- ggplot(Data_Plot, aes(x = comparison, y = gene)) + 
  geom_point(aes(size = logp, color = avg_log2FC)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    title = "DEG Analyis for comparing between B Cells vs T cells, and Monocytes vs T Cells",
    x = "Cell Type",
    y = "Gene",
    color = "Log2 Fold Change",
    size = "-Log10 Adjusted P-value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
print(DEG_Comparison)
ggsave("~/Desktop/SC/Figures_Assignment_1/DEG_Comparison.png", plot = DEG_Comparison, width = 10, height = 10, dpi = 300)

head(monocytes@meta.data)
monocytes$Group <- merged_with_batch$Group[Cells(monocytes)]
table(monocytes$Group)


#activating "Group", as comparison is between two sample group BMMC vs CD34
Idents(merged_with_batch) <- "Group"
DEG_BMMC_CD34 <- FindMarkers(object = merged_with_batch, ident.1 = "BMMC", ident.2 = "CD34", min.pct = 0.25, logfc.threshold = 0.25)
DEG_BMMC_CD34_5 <- DEG_BMMC_CD34 %>%
  arrange(p_val_adj) %>%
  head(5)
print(DEG_BMMC_CD34_5)

DEG_Monocytes_BMMC_CD34 <- FindMarkers(
  object = monocytes,
  ident.1 = "BMMC",
  ident.2 = "CD34",
  group.by = "Group",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

Top5_Monocytes_BMMC_CD34 <- DEG_Monocytes_BMMC_CD34 %>%
  arrange(p_val_adj) %>%
  head(5)
print(Top5_Monocytes_BMMC_CD34)
VlnPlot_Monocytes_Groups <- VlnPlot(monocytes, features = c("CD3D", "GATA2", "CD14"), group.by = "Group", pt.size = 0.1)
print(VlnPlot_Monocytes_Groups)
ggsave("~/Desktop/SC/Figures_Assignment_1/VlnPlot_Monocytes_Groups.png", plot = VlnPlot_Monocytes_Groups, width = 20, height = 10, dpi = 300)

#checking available databases 
available_databases <- enrichR::listEnrichrDbs()
top_genes <- rownames(DEG_BMMC_CD34[DEG_BMMC_CD34$p_val_adj < 0.05, ])
Enrichment_results <-  enrichr(top_genes, databases = c("GO_Biological_Process_2023", "KEGG_2021_Human", "Human_Gene_Atlas", "ENCODE_TF_ChIP-seq_2015", "MSigDB_Hallmark_2020"))	

View(Enrichment_results[["GO_Biological_Process_2023"]])
View(Enrichment_results[["KEGG_2021_Human"]])
View(Enrichment_results[["Human_Gene_Atlas"]])
View(Enrichment_results[["ENCODE_TF_ChIP-seq_2015"]])
View(Enrichment_results[["MSigDB_Hallmark_2020"]])

GO_Enrichr <- plotEnrich(Enrichment_results[["GO_Biological_Process_2023"]], showTerms = 30, title = "GO Terms") #pathway analysis, showing 30 most significant terms
print(GO_Enrichr)
ggsave("~/Desktop/SC/Figures_Assignment_1/GO_Enrichr.png", plot = GO_Enrichr, width = 10, height = 10, dpi = 300)
KEGG_Enrichr <- plotEnrich(Enrichment_results[["KEGG_2021_Human"]], showTerms = 30, title = "KEGG Analysis")
print(KEGG_Enrichr)
ggsave("~/Desktop/SC/Figures_Assignment_1/KEGG_Enrichr.png", plot = KEGG_Enrichr, width = 10, height = 10, dpi = 300)
HGA_Enrichr <- plotEnrich(Enrichment_results[["Human_Gene_Atlas"]], showTerms = 30, title = "Human Gene Atlas")
print(HGA_Enrichr)
ggsave("~/Desktop/SC/Figures_Assignment_1/HGA_Enrichr.png", plot = HGA_Enrichr, width = 10, height = 10, dpi = 300)
ENCODE_Enrichr <- plotEnrich(Enrichment_results[["ENCODE_TF_ChIP-seq_2015"]], showTerms = 30, title = "ENCODE Analysis")
print(ENCODE_Enrichr)
ggsave("~/Desktop/SC/Figures_Assignment_1/ENCODE_Enrichr.png", plot = ENCODE_Enrichr, width = 10, height = 10, dpi = 300)
MSigDB_Enrichr <- plotEnrich(Enrichment_results[["MSigDB_Hallmark_2020"]], showTerms = 30, title = "MSigDB Analysis")
print(MSigDB_Enrichr)
ggsave("~/Desktop/SC/Figures_Assignment_1/MSigDB_Enrichr.png", plot = MSigDB_Enrichr, width = 10, height = 10, dpi = 300)
#DEenrichRPlot(DEG_BMMC_CD34, ident.1 = "BMMC", ident.2 = "CD34", assay = "RNA", databases = "GO_Biological_Process_2023")


top_pathway_GO <- Enrichment_results[["GO_Biological_Process_2023"]] %>%
  arrange(Adjusted.P.value) %>%
  head(1)
print(top_pathway_GO)


top_pathway_KEGG <- Enrichment_results[["KEGG_2021_Human"]] %>%
  arrange(Adjusted.P.value) %>%
  head(1)
print(top_pathway_KEGG)


top_pathway_HGA <- Enrichment_results[["Human_Gene_Atlas"]] %>%
  arrange(Adjusted.P.value) %>%
  head(1)
print(top_pathway_HGA)


top_pathway_ENCODE <- Enrichment_results[["ENCODE_TF_ChIP-seq_2015"]] %>%
  arrange(Adjusted.P.value) %>%
  head(1)
print(top_pathway_ENCODE)


top_pathway_MSigDB <- Enrichment_results[["MSigDB_Hallmark_2020"]] %>%
  arrange(Adjusted.P.value) %>%
  head(1)
print(top_pathway_MSigDB)



#################################################WEEK4##########################################################
# Check current identities
Idents(merged_with_batch)
Idents(merged_with_batch) <- "seurat_clusters" #set to seurat_clusters
unique(Idents(merged_with_batch))
selected_clusters <- subset(merged_with_batch, idents = c(0, 1, 2)) #Selection the clusters for the cell types, HSC, LMPP, CLP respectively.
table(Idents(selected_clusters))

cds <- as.cell_data_set(selected_clusters, assay = "integrated") #convertinf the selected cluster to cell dataset in the integrated assay
colData(cds) #retrieving column data
fData(cds) #retrieving featur data
rownames(fData(cds))[1:10] #display first 10 rows
fData(cds)$gene_short_name <- rownames(fData(cds)) #assigning gene short names to feature data based on row names
fData(cds)
counts(cds)
recreate.partition <- c(rep(1,length(cds@colData@rownames))) #creating uniform partition vector for all cells
names(recreate.partition) <-  cds@colData@rownames #naming the partitions using cell row names
recreate.partition <- as.factor(recreate.partition) #converting partition to factor
recreate.partition
cds@clusters$UMAP$partitions <- recreate.partition #updating the UMAP with the created partition
cds@clusters$UMAP$partitions
list_cluster <- selected_clusters@active.ident
cds@clusters$UMAP$clusters <- list_cluster

cds@int_colData@listData$reducedDims$UMAP <- selected_clusters@reductions$umap@cell.embeddings #updating UMAP cell embeddings from Seurat reductions

cluster.before.trajectory <- plot_cells(cds, color_cells_by = 'seurat_clusters', label_groups_by_cluster = FALSE, group_label_size = 5) +
  scale_color_manual(values = c('red', 'orange', 'cyan')) + theme(legend.position = "right")
print(cluster.before.trajectory)

cds <- learn_graph(cds, use_partition = FALSE)
trajectory_cds <- plot_cells(cds, color_cells_by = 'seurat_clusters', label_groups_by_cluster = FALSE, label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE, group_label_size = 5)
print(trajectory_cds)
ggsave("~/Desktop/SC/Figures_Assignment_1/trajectory_cds.png", plot = trajectory_cds, width = 10, height = 10, dpi = 300)

root_cells <- colnames(cds)[colData(cds)$seurat_clusters == "1"] #Manually selecting '0'
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = root_cells)

manual_selection <- plot_cells(cds, color_cells_by = 'pseudotime', label_groups_by_cluster = FALSE, label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE,)
print(manual_selection)
ggsave("~/Desktop/SC/Figures_Assignment_1/Manual_selection.png", plot = manual_selection, width = 10, height = 10, dpi = 300)

#automatic selection
get_earliest_principal_node <- function(cds){
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex #extracting closest vertex data and convertinf to matrix using cell names
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  
  node_frequencies <- table(closest_vertex) #cslculating frequency of each vertex
  
  root_pr_node <- names(which.max(node_frequencies))    #highest frequency vertex
  root_pr_node <- paste0("Y_", root_pr_node)             # Printing all valid nodes, and here since the numbers starts with Y_
  return(root_pr_node)                                  # Checking the value of the root_pr_node, it showed 81, added Y_ to match the graph node name
  
}

root_pr_node <- get_earliest_principal_node(cds)
cds <- order_cells(cds, root_pr_nodes = root_pr_node)

automatic_selection <- plot_cells(cds, color_cells_by = 'pseudotime', label_groups_by_cluster = FALSE, label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE)
print(automatic_selection)
ggsave("~/Desktop/SC/Figures_Assignment_1/Automatic_selection.png", plot = automatic_selection, width = 10, height = 10, dpi = 300)

#CellChat
# Converting the RNA assay from Assay5 (Seurat v5) to Assay (Seurat v3)
merged_with_batch[["RNA_v3"]] <- as(object = merged_with_batch[["RNA"]], Class = "Assay") #converted v5 to v3 because i faced issues with v5, only for cellchat

# Confirming the conversion
class(merged_with_batch[["RNA_v3"]])
DefaultAssay(merged_with_batch) <- "RNA_v3"

data_matrix <- GetAssayData(merged_with_batch, assay = "RNA_v3", slot = "counts")
print(dim(data_matrix))  # show genes as rows and cells as columns
meta_data <- merged_with_batch@meta.data

#common cell barcodes
common_barcodes <- intersect(colnames(data_matrix), rownames(meta_data))

#Subsetting data_matrix and meta_data
data_matrix <- data_matrix[, common_barcodes]
meta_data <- meta_data[common_barcodes, , drop = FALSE]
rownames(meta_data) <- colnames(data_matrix)

# Adding manual annotations
manual_annotations <- c(
  "0" = "HSC(HSC)",
  "1" = "LMPP(LMPP)",
  "2" = "CLP(CLP)",
  "3" = "GMP(GMP)",
  "4" = "CMP(CMP)",
  "5" = "B Cells(B)",
  "6" = "Pre B-cell Progenitors(PreB)",
  "7" = "Plasma Cells(Plasma)",
  "8" = "CD8+ T Cells(CD8)",
  "9" = "CD4+ T Cells(CD4)",
  "10" = "NK Cells(NK)",
  "11" = "Erythrocytes(Ery)",
  "12" = "cDC(cDC)",
  "13" = "pDC(PDC)",
  "14" = "CD14+ Monocytes(Mono CD14)",
  "15" = "CD16+ Monocytes(Mono CD16)",
  "16" = "Basophils(Bas)"
)

# Applying manual labels
meta_data$manual_labels <- manual_annotations[as.character(Idents(merged_with_batch))]
meta_data$manual_labels[is.na(meta_data$manual_labels)] <- "Unknown"
Idents(merged_with_batch) <- manual_annotations[as.character(Idents(merged_with_batch))]
unique(meta_data$manual_labels)
# Add manual_labels to Seurat object's meta.data
merged_with_batch@meta.data$manual_labels <- manual_annotations[as.character(Idents(merged_with_batch))]



#CellChat object
cellchat <- createCellChat(object = data_matrix, meta = meta_data, group.by = "manual_labels")
unique(cellchat@meta[["Group"]])

CellChatDB <- CellChatDB.human #loading human intraction database in to CellChatDB,assigning it to cellchat object
cellchat@DB <- CellChatDB

#CellChat analysis
cellchat <- subsetData(cellchat) #subsetting based on predefined condition
cellchat <- identifyOverExpressedGenes(cellchat) #identifying genes that are overexpresed in the dataset
cellchat <- identifyOverExpressedInteractions(cellchat) #identifying interactions of overexpressed genes
cellchat <- computeCommunProb(cellchat) #computing probabilities
cellchat <- filterCommunication(cellchat, min.cells = 10) #fitering out communicaiton signals occuring in fewer than 10 cells. added only for this
cellchat <- computeCommunProbPathway(cellchat) #computing communication probabilities for pathways
cellchat <- aggregateNet(cellchat) #aggregating data to summarise the network interaction

cell_types_BMMC <- unique(meta_data[meta_data$Group == "BMMC", "manual_labels"])
cell_types_CD34 <- unique(meta_data[meta_data$Group == "CD34", "manual_labels"])
# Recomputing the intersection
common_cell_types <- intersect(cell_types_BMMC, cell_types_CD34)
print(common_cell_types)
BMMC_cells <- subset(merged_with_batch, subset = Group == "BMMC" & manual_labels %in% common_cell_types)
CD34_cells <- subset(merged_with_batch, subset = Group == "CD34" & manual_labels %in% common_cell_types)

cellchat_BMMC <- createCellChat(object = GetAssayData(BMMC_cells, slot = "counts"), meta = BMMC_cells@meta.data, group.by = "manual_labels")
cellchat_CD34 <- createCellChat(object = GetAssayData(CD34_cells, slot = "counts"), meta = CD34_cells@meta.data, group.by = "manual_labels")

cellchat_BMMC@DB <- CellChatDB.human # Setting the CellChat database for both objects
cellchat_CD34@DB <- CellChatDB.human

# BMMC analysis
cellchat_BMMC <- subsetData(cellchat_BMMC)
cellchat_BMMC <- identifyOverExpressedGenes(cellchat_BMMC)
cellchat_BMMC <- identifyOverExpressedInteractions(cellchat_BMMC)
cellchat_BMMC <- computeCommunProb(cellchat_BMMC, population.size = FALSE)
cellchat_BMMC <- filterCommunication(cellchat_BMMC)
cellchat_BMMC <- computeCommunProbPathway(cellchat_BMMC)
cellchat_BMMC <- aggregateNet(cellchat_BMMC)

# CD34 analysis
cellchat_CD34 <- subsetData(cellchat_CD34)
cellchat_CD34 <- identifyOverExpressedGenes(cellchat_CD34)
cellchat_CD34 <- identifyOverExpressedInteractions(cellchat_CD34)
cellchat_CD34 <- computeCommunProb(cellchat_CD34, population.size = FALSE)
cellchat_CD34 <- filterCommunication(cellchat_CD34)
cellchat_CD34 <- computeCommunProbPathway(cellchat_CD34)
cellchat_CD34 <- aggregateNet(cellchat_CD34)

pathways_BMMC <- cellchat_BMMC@netP$pathways  #listingsignaling pathways for each group
pathways_CD34 <- cellchat_CD34@netP$pathways

#common pathways between BMMC and CD34
common_pathways <- intersect(pathways_BMMC, pathways_CD34)
print(common_pathways)  # Verify the common signaling pathways

#Heatmap_cellchat
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Reds", title.name = "Number of Interactions - cellchat")
netVisual_heatmap(cellchat, measure = "weigh", color.heatmap = "Reds", title.name = "Interaction Strength - cellchat")

#Heatmap_BMMC
netVisual_heatmap(cellchat_BMMC, measure = "count", color.heatmap = "Reds", title.name = "Number of Interactions - BMMC")
netVisual_heatmap(cellchat_BMMC, measure = "weight", color.heatmap = "Purples", title.name = "Interaction Strength - BMMC")

#Heatmap_CD34
netVisual_heatmap(cellchat_CD34, measure = "count", color.heatmap = "Blues", title.name = "Number of Interactions - CD34")
netVisual_heatmap(cellchat_CD34, measure = "weight", color.heatmap = "Greens", title.name = "Interaction Strength - CD34")

selected_pathway <- "MIF"  

pathway_BMMC <- netAnalysis_contribution(cellchat_BMMC, signaling = selected_pathway)
pathway_CD34 <- netAnalysis_contribution(cellchat_CD34, signaling = selected_pathway)

netVisual_aggregate(cellchat_BMMC, signaling = "MIF", layout = "circle")
title("MIF Pathway - BMMC")

netVisual_aggregate(cellchat_CD34, signaling = "MIF", layout = "circle")
title("MIF Pathway - CD34")



saveRDS(BMMC_D1T1, file ="/Users/adwitiyaarghapb/Desktop/SC/SC_Assignment.R")

