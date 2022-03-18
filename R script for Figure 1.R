library(dplyr)
library(Seurat)  #updates from 3.2.3 to 4.0.3
library(patchwork)
library(cowplot)
library(ggplot2)
library(Signac)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(psych)
library(qgraph)
library(igraph)
library(GSVA)
library(GSEABase)
library(limma)
library(hdf5r)


#input h5 data
T1L.data <- Read10X_h5("/T1L.h5", use.names = T)
T1R.data <- Read10X_h5("/T1R.h5", use.names = T)
T2L.data <- Read10X_h5("/T2L.h5", use.names = T)
T2R.data <- Read10X_h5("/T2R.h5", use.names = T)
T3L.data <- Read10X_h5("/T3L.h5", use.names = T)
T3R.data <- Read10X_h5("/T3R.h5", use.names = T)
NT.data <- Read10X_h5("/NT.h5", use.names = T)

# Create Seurat object and quality control
T1L <- CreateSeuratObject(counts = T1L.data, project = "Thyroid_L1", min.cells = 3, min.features = 200)
T1L$stim <- "LEFT1"
T1L[["percent.mt"]] <- PercentageFeatureSet(T1L, pattern = "^MT-")
VlnPlot(T1L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T1L <- subset(T1L, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)

T2L <- CreateSeuratObject(counts = T2L.data, project = "Thyroid_L2", min.cells = 3, min.features = 200)
T2L$stim <- "LEFT2"
T2L[["percent.mt"]] <- PercentageFeatureSet(T2L, pattern = "^MT-")
VlnPlot(T2L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T2L <- subset(T2L, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)

T3L <- CreateSeuratObject(counts = T3L.data, project = "Thyroid_L3", min.cells = 3, min.features = 200)
T3L$stim <- "LEFT3"
T3L[["percent.mt"]] <- PercentageFeatureSet(T3L, pattern = "^MT-")
VlnPlot(T3L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T3L <- subset(T3L, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)

NT <- CreateSeuratObject(counts = NT.data, project = "Thyroid_L4", min.cells = 3, min.features = 200)
NT$stim <- "Normal"
NT[["percent.mt"]] <- PercentageFeatureSet(NT, pattern = "^MT-")
VlnPlot(NT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
NT <- subset(NT, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)

T1R <- CreateSeuratObject(counts = T1R.data, project = "Thyroid_R1", min.cells = 3, min.features = 200)
T1R$stim <- "RIGHT1"
T1R[["percent.mt"]] <- PercentageFeatureSet(T1R, pattern = "^MT-")
VlnPlot(T1R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T1R <- subset(T1R, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)

T2R <- CreateSeuratObject(counts = T2R.data, project = "Thyroid_R2", min.cells = 3, min.features = 200)
T2R$stim <- "RIGHT2"
T2R[["percent.mt"]] <- PercentageFeatureSet(T2R, pattern = "^MT-")
VlnPlot(T2R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T2R <- subset(T2R, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)

T2R <- CreateSeuratObject(counts = T2R.data, project = "Thyroid_R3", min.cells = 3, min.features = 200)
T2R$stim <- "T2R"
T2R[["percent.mt"]] <- PercentageFeatureSet(T2R, pattern = "^MT-")
VlnPlot(T2R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T2R <- subset(T2R, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)


##Perform integration using CCA method

thyroid.anchors <- FindIntegrationAnchors(object.list = list(T1L, T2L, T3L, NT, T1R, T2R, T3R), dims = 1:20)
thyroid.combined <- IntegrateData(anchorset = thyroid.anchors, dims = 1:20)
thyroid.combined

#Perform an integrated analysis

DefaultAssay(thyroid.combined) <- "integrated"
#NormalizeData
thyroid.combined <- NormalizeData(thyroid.combined, normalization.method = "LogNormalize", scale.factor = 10000)
#FindVariableFeatures
thyroid.combined <- FindVariableFeatures(thyroid.combined, selection.method = "vst", nfeatures = 2000)

#to check nFeature_RNA", "nCount_RNA", "percent.mt
VlnPlot(thyroid.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Run the standard workflow for visualization and clustering
thyroid.combined <- ScaleData(thyroid.combined, verbose = FALSE)
thyroid.combined <- RunPCA(thyroid.combined, npcs = 30, verbose = FALSE)


#Determine the 'dimensionality' of the dataset 
thyroid.combined <- JackStraw(thyroid.combined, num.replicate = 100)
thyroid.combined <- ScoreJackStraw(thyroid.combined, dims = 1:20)
JackStrawPlot(thyroid.combined, dims = 1:15)
ElbowPlot(thyroid.combined)


# t-SNE and Clustering
thyroid.combined <- RunUMAP(thyroid.combined, reduction = "pca", dims = 1:15)
thyroid.combined <- FindNeighbors(thyroid.combined, reduction = "pca", dims = 1:15)
thyroid.combined <- FindClusters(thyroid.combined, resolution = 0.8)
thyroid.combined <- RunTSNE(thyroid.combined, reduction = "pca", dims = 1:15)

#name the clusters
thyroid.combined <- RenameIdents(thyroid.combined, `0` = "Follicular_cells", `1` = "Pericyte", `2` = "Follicular_cells", 
                                 `3` = "Follicular_cells", `4` = "T_cells", `5` = "Follicular_cells", `6` = "Follicular_cells", `7` = "T_cells", `8` = "Myeloid_cells", `9` = "Follicular_cells", 
                                 `10` = "Endothelial_cells", `11` = "Follicular_cells", `12` = "T_cells", `13` = "Pericyte", `14` = "T_cells", `15` = "T_cells",
                                 `16` = "Endothelial_cells", `17` = "T_cells", `18` = "T_cells", `19` = "B_cells", `20` = "Fibroblast", `21` = "Endothelial_cells", `22` = "B_cells", `23` = "Endothelial_cells"
                                 , `24` = "Pericyte", `25` = "B_cells", `26` = "Mast_cells", `27` = "B_cells", `28` = "Endothelial_cells", `29` = "Endothelial_cells"
                                 , `30` = "Pericyte", `31` = "B_cells")

thyroid.combined@meta.data$celltype <- factor(thyroid.combined@meta.data$seurat_clusters, 
                                                     levels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,31), 
                                                     labels = c("Follicular_cells", "Pericyte",  "Follicular_cells", 
                                                                "Follicular_cells",  "T_cells",  "Follicular_cells",  "Follicular_cells",  "T_cells", "Myeloid_cells", "Follicular_cells", 
                                                                "Endothelial_cells", "Follicular_cells", "T_cells", "Pericyte", "T_cells",  "T_cells",
                                                                "Endothelial_cells",  "T_cells",  "T_cells",  "B_cells",  "Fibroblast",  "Endothelial_cells",  "B_cells",  "Endothelial_cells"
                                                                ,  "Pericyte",  "B_cells",  "Mast_cells",  "B_cells",  "Endothelial_cells",  "Endothelial_cells"
                                                                ,  "Pericyte", "B_cells"))


thyroid.combined@meta.data$stim <- factor(thyroid.combined@meta.data$stim, 
                                          levels = c("LEFT1", "LEFT2", "LEFT3", "LEFT4", "RIGHT1", "RIGHT2", "RIGHT3"), 
                                          labels = c("T1L", "T2L", "T3L", "NT", "T1R", "T2R", "T3R"))
View(thyroid.combined@meta.data)


# Visualization for Figure 1C and 1D
DimPlot(thyroid.combined, reduction = "umap", group.by = "stim")
DimPlot(thyroid.combined, reduction = "umap", label = TRUE)


#VlnPlot for Figure 1E
VlnPlot(thyroid.combined, features = c("TG", "HIGD1B","CD3D", "CD68", "VWF", "CD79A","PDGFRA","TPSAB1"), stack = T, flip = T)


# find markers for every cluster compared to all remaining cells, report only the positive ones
thyroid.combined.markers <- FindAllMarkers(thyroid.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(thyroid.combined.markers,"Supplementary_File_1_DEGs_among_8_known_cell_clusters.csv")

#heatmap Figure 1F
DoHeatmap(thyroid.combined, features = c("TG","KRT18","KRT19","TSHR","KRT7"
                                         ,"HIGD1B","CSRP2","CACNB2","COL25A1","RGS5"
                                         ,"CD3G","GNLY","PDCD1","CD8B","FOXP3"
                                         ,"CD86","LYZ","HLA-DRA","APOC1","S100A8"
                                         ,"VWF","ARL15","PLPP1","PTPRG","STC1"
                                         ,"CD79A","DERL3","BANK1","MS4A1","LY9"
                                         ,"PDGFRA","COL3A1","COL1A1","MGP","COL1A2"
                                         ,"KIT","CSF1","S100B"), angle = 90) + NoLegend()


##Fraction of cell types in each sample for Figure 1G
table(thyroid.combined$stim,thyroid.combined$celltype)
write.csv(table(thyroid.combined$stim,thyroid.combined$celltype),"celltypefrequency.csv")





