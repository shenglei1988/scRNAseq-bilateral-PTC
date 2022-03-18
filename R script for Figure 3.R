Subset_Cells <- subset(thyroid.combined, idents = "T_cells")
Subset_Cells <- ScaleData(Subset_Cells, verbose = FALSE)
#FindVariableFeatures
Subset_Cells <- FindVariableFeatures(Subset_Cells, selection.method = "vst", nfeatures = 2000)
#
Subset_Cells <- RunPCA(Subset_Cells, npcs = 50, verbose = FALSE)
Subset_Cells <- FindNeighbors(Subset_Cells, reduction = "pca", dims = 1:10)
Subset_Cells <- FindClusters(Subset_Cells, resolution = 0.5)
Subset_Cells <- RunUMAP(Subset_Cells, reduction = "pca", dims = 1:5)
Subset_Cells <- RunTSNE(Subset_Cells, reduction = "pca", dims = 1:5)



Subset_Cells <- subset(Subset_Cells, idents = c(0, 1, 2, 3, 4, 5, 6, 8, 9, 10))   #excluding cluster 7 like follicular cells (n=368)
Subset_Cells <- ScaleData(Subset_Cells, verbose = FALSE)
#FindVariableFeatures
Subset_Cells <- FindVariableFeatures(Subset_Cells, selection.method = "vst", nfeatures = 2000)
#
Subset_Cells <- RunPCA(Subset_Cells, npcs = 50, verbose = FALSE)
Subset_Cells <- FindNeighbors(Subset_Cells, reduction = "pca", dims = 1:10)
Subset_Cells <- FindClusters(Subset_Cells, resolution = 0.5)
Subset_Cells <- RunUMAP(Subset_Cells, reduction = "pca", dims = 1:5)
Subset_Cells <- RunTSNE(Subset_Cells, reduction = "pca", dims = 1:5)

head(Subset_Cells@meta.data)
Subset_Cells@meta.data$celltype <- factor(Subset_Cells@meta.data$seurat_clusters, 
                                                 levels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), 
                                                 labels = c("CD8 TRM.1","Treg.1","Cytotoxic CD8.1","Activated CD4.1","CD8 TRM.2","Dying T cell","Activated CD4.2","Treg.2","CD8 TRM.3","Activated CD4.3") )

#UMAP Figure 3A
DimPlot(Subset_Cells, reduction = "umap", group.by = "celltype",label = TRUE)  # 7.2*6


#Featureplot Figure 3B
FeaturePlot(Subset_Cells, features = c("CD8B", "CTLA4", "ZNF683", "ITGAE", "PRF1", "NKG7", "FOXP3", "IL2RA", "CD69", "CD40LG"), cols = c("lightgrey","red"),reduction = "umap", min.cutoff = "q5")   #20*14



#Heatmap Figure 3C
Subset_Cells$celltype <- factor(x =Subset_Cells$celltype, levels = c('CD8 TRM.1','CD8 TRM.2','CD8 TRM.3','Cytotoxic CD8.1','Treg.1','Treg.2','Activated CD4.1','Activated CD4.2','Activated CD4.3',"Dying T cell"))

DoHeatmap(Subset_Cells, features = c("CD3D", "CD8B", "CTLA4", "IFNG", "TNF", "CCL5","CCL4", "CCL3L1","CCL4L2",
                                     "CD69", "RGS1", "IL2RA", "FASLG", "CD40LG", "ICOS", "TNFRSF4", "TNFRSF9", "HLA-DRA",
                                     "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "NKG7", "GNLY", "LAMP1",
                                     "NCAM1", "NCR1", "FCGR3A", "KLRB1",
                                     "FOXP3", "ZNF683", "TOX", "TCF7", "TBX21", "ECOMES", "PRDM1", "GATA3", "LEF1", 
                                     "FOS"),angle = 90, label=F, size = 5, group.by ="celltype")


#vlnplot Figure 3D
VlnPlot(Subset_Cells, features = c("CD4","CD8B","IFNG", "FOS", "PRF1", "NKG7", "GZMA", "FOXP3", "CD40LG", "CD69", "ZNF683", "ITGAE", "PRMD1", "RGS1", "CXCR6"),stack = T, flip = T, group.by ="seurat_clusters")   

##figure 3E and 3F

table(Subset_Cells$stim,Subset_Cells$celltype)







