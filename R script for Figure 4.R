Subset_Cells <- subset(thyroid.combined, idents = "Myeloid_cells")
Subset_Cells <- ScaleData(Subset_Cells, verbose = FALSE)
#FindVariableFeatures
Subset_Cells <- FindVariableFeatures(Subset_Cells, selection.method = "vst", nfeatures = 2000)
#
Subset_Cells <- RunPCA(Subset_Cells, npcs = 50, verbose = FALSE)
Subset_Cells <- FindNeighbors(Subset_Cells, reduction = "pca", dims = 1:10)
Subset_Cells <- FindClusters(Subset_Cells, resolution = 0.5)
Subset_Cells <- RunUMAP(Subset_Cells, reduction = "pca", dims = 1:5)
Subset_Cells <- RunTSNE(Subset_Cells, reduction = "pca", dims = 1:5)


Subset_Cells <- subset(Subset_Cells, idents = c(0, 1, 2, 3, 4, 5, 6))   #excluding cluster 7 like follicular cells
Subset_Cells <- ScaleData(Subset_Cells, verbose = FALSE)
#FindVariableFeatures
Subset_Cells <- FindVariableFeatures(Subset_Cells, selection.method = "vst", nfeatures = 2000)
#
Subset_Cells <- RunPCA(Subset_Cells, npcs = 50, verbose = FALSE)
Subset_Cells <- FindNeighbors(Subset_Cells, reduction = "pca", dims = 1:10)
Subset_Cells <- FindClusters(Subset_Cells, resolution = 0.5)
Subset_Cells <- RunUMAP(Subset_Cells, reduction = "pca", dims = 1:5)
Subset_Cells <- RunTSNE(Subset_Cells, reduction = "pca", dims = 1:5)

#UMAP Supplementary Figure S5A and S5B
DimPlot(Subset_Cells, reduction = "umap")
DimPlot(Subset_Cells, reduction = "umap", group.by = "stim")

Subset_Cells <- RenameIdents(Subset_Cells, `0` = "M1-like", `1` = "M2-like", `2` = "Marker-low cells", 
                             `3` = "CD14+",`4` = "CD14+",`5` = "CD14+/CD16+ monocyte",`6` = "CD14+/CD16+ monocyte")


#Name clusters 
head(Subset_Cells@meta.data)
Subset_Cells@meta.data$celltype <- factor(Subset_Cells@meta.data$seurat_clusters, 
                                          levels = c(0, 1, 2, 3, 4, 5, 6), 
                                          labels = c("M1-like", "M2-like", "Marker-low cells", 
                                                     "CD14+","CD14+", "CD14+/CD16+ monocyte","CD14+/CD16+ monocyte") )


# Figure 4A and 4B
DimPlot(Subset_Cells, reduction = "umap", group.by = "stim",cols = c("red", "red", "red", "green", "red", "red", "red"))
DimPlot(Subset_Cells, reduction = "umap", cols = c("red", "green", "orange", "blue", "purple"))

#Figure 4C
VlnPlot(Subset_Cells, features = c("TNF","IL1B","CTSD","FN1","S100A8","S100A9","CD14","FCGR3A"), idents=c('M1-like','M2-like','CD14+','CD14+/CD16+ monocyte'),stack = T, flip = T)   #5*6

##headmap Figure 4D
mi <- subset(Subset_Cells, idents = c('M1-like',"M2-like"))
DoHeatmap(Subset_Cells, features = c("IL1B","TNF","HLA-DPB1","CD86","IL2RA",
                                     "CCL18","CCL23","CTSD","FN1","GAS7","HMOX1","PPARG","LIPA","MS4A4A"),angle = 90, label=F, size = 5,group.colors = c("red", "green"), cells = colnames(mi))

##data for Figure 4E
table(Subset_Cells$stim,Subset_Cells$celltype)


##UMAP for Figure 4F
DimPlot(Subset_Cells, reduction = "umap", split.by = "stim", cells.highlight = CellsByIdentities(Subset_Cells, idents=c('M1-like','M2-like')), cols.highlight =c("green", "red"))




###Monocle for Figure 4I to 4M
install.packages("devtools")
devtools::install_github("cole-trapnell-lab/monocle-release@develop")

library(monocle)


view(Subset_Cells@meta.data )
sce=Subset_Cells


##
sample_ann <- sce@meta.data 
#sample_ann$celltype <- sample_ann@active.ident  ##sample_ann$celltype <- Idents(sce)
head(sample_ann)
dim(sample_ann)
##
gene_ann <- data.frame(gene_short_name = row.names(sce@assays$RNA),row.names = row.names(sce@assays$RNA))
##
head(gene_ann)
dim(gene_ann)

ct=as.data.frame(sce@assays$RNA@counts)
ct[1:4,1:4]
dim(ct)

pd <- new('AnnotatedDataFrame', data = sample_ann) 
fd <- new('AnnotatedDataFrame', data = gene_ann)

sc_cds <- newCellDataSet(as.matrix(ct),
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())


sc_cds
sc_cds <- detectGenes(sc_cds, min_expr = 1)
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 10,]   


cds <- sc_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds


disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds)
plot_pc_variance_explained(cds,return_all = F)


cds <- reduceDimension(cds, max_components = 2, num_dim = 6, reduction_method = 'tSNE',verbose = T)
cds <- clusterCells(cds, num_clusters = 6)  
plot_cell_clusters(cds, 1, 2)



table(pData(cds)$Cluster)
colnames(pData(cds))
table(pData(cds)$seurat_clusters)
table(pData(cds)$Cluster,pData(cds)$seurat_clusters)
table(pData(cds)$Cluster,pData(cds)$celltype)
plot_cell_clusters(cds, 1, 2)


save(cds, file = 'input_cds.Rdata')
#加载数据
load('input_cds.Rdata')


pData(cds)$Cluster=pData(cds)$seurat_clusters 
table(pData(cds)$Cluster)
view(pData(cds))

Sys.time()
diff_test_res <- differentialGeneTest(cds,fullModelFormulaStr = "~Cluster")
Sys.time()

head(diff_test_res)

#select genes that are significant at a FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes = sig_genes[order(sig_genes$pval),]
head(sig_genes[,c("gene_short_name","pval","qval")])

cg=as.character(head(sig_genes$gene_short_name))
plot_genes_jitter(cds[cg,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow = 3,
                  ncol = NULL)

cg2=as.character(tail(sig_genes$gene_short_name))
plot_genes_jitter(cds[cg2,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow = 3,
                  ncol = NULL)


plot_genes_jitter(cds[cg,],
                  grouping = "State",
                  color_by = "State",
                  nrow = 3,
                  ncol = NULL)


ordering_genes <- row.names(subset(diff_test_res, qval< 0.01))
ordering_genes
cds <- setOrderingFilter(cds,ordering_genes)
plot_ordering_genes(cds)

Sys.time()
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
Sys.time()



cds <- orderCells(cds)


plot_cell_trajectory(cds, color_by = "Cluster")
ggsave('monocle_plot_cell_trajectory_for_seurat')

#
phe=pData(cds)
boxplot(phe$Pseudotime,phe$Cluster) 



my_cds_subset=cds
head(pData(my_cds_subset))

Sys.time()
my_pseudotime_de <- differentialGeneTest(my_cds_subset,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 1)   

Sys.time()
head(my_pseudotime_de)
save(my_cds_subset,my_pseudotime_de, file = 'output_of_monocle.Rdata')

load(file = 'output_of_monocle.Rdata')
cds=my_cds_subset
colnames(pData(cds))
view(pData(cds))


library(ggsci)
p1=plot_cell_trajectory(cds,color_by = "Cluster")    
p1

p6=plot_cell_trajectory(cds,color_by = "Cluster") + facet_wrap(~Cluster, nrow = 4)
p6
ggsave('trajectory_by_cluster.pdf')

plot_cell_trajectory(cds,color_by = "celltype")


p2=plot_cell_trajectory(cds,color_by = "Pseudotime") 
p2
ggsave('trajectory_by_Pseudotime.pdf')

p3=plot_cell_trajectory(cds,color_by = "State") 
p3
ggsave('trajectory_by_State.pdf')

p4=plot_cell_trajectory(cds,color_by = "stim") 
p4

p5=plot_cell_trajectory(cds,color_by = "stim") + facet_wrap(~stim, nrow = 1)
p5



library(patchwork)
p1 + p2/p3                 

phe=pData(cds)
head(phe)
table(phe$State,phe$Cluster)
table(phe$stim,phe$Cluster)

write.csv(table(phe$State,phe$Cluster),"StatevsCluster.csv")
write.csv(table(phe$stim,phe$Cluster),"stimvsCluster.csv")
write.csv(table(phe$stim,phe$State),"stimvsState.csv")


library(dplyr)
my_pseudotime_de %>% arrange(qval) %>% head()

#save the top 6 genes
my_pseudotime_de %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_gene
my_pseudotime_gene = my_pseudotime_gene[,1]
my_pseudotime_gene
plot_gene_in_pseudotime(my_cds_subset[my_pseudotime_gene,])+ scale_color_npg()
ggsave('monocle_top6_pseudotime_by_state.pdf')

plot_genes_jitter(my_cds_subset[my_pseudotime_gene,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow = 3,
                  ncol = NULL)+scale_color_nejm()
ggsave('monocle_top6_pseudotime_by_cluster.pdf')


#cluster the top 50 genes that vary as a function of pseudotime
my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster[,1]
gene_to_cluster
colnames(pData(my_cds_subset))
dim(pData(my_cds_subset))

table(pData(my_cds_subset)$Cluster,pData(my_cds_subset)$State)
ac=pData (my_cds_subset)[,c('Cluster','State','Pseudotime')]   #Cluster替换celltype
head(ac)
dim(ac)
#这个热图绘制的并不是纯粹的细胞基因表达量矩阵，而是被pseudotime好了的100例
my_pseudotime_cluster <- plot_pseudotime_heatmap (my_cds_subset[gene_to_cluster,],
                                                  #num_clusters =2,
                                                  #add_annotation_col = ac,
                                                  show_rownames = TRUE,
                                                  return_heatmap = TRUE)

my_pseudotime_cluster



pdf('monocle_top50_heatmap')
print(my_pseudotime_cluster)

my_pseudotime_de %>% arrange(qval) %>% head(100) %>% select(gene_short_name) -> my_pseudotime_gene
my_pseudotime_gene = my_pseudotime_gene[,1]
my_pseudotime_gene

library(pheatmap)

n=t(scale(t( counts[my_pseudotime_gene,])))
n[n>2] = 2
n[n<-2] = -2
n[1:4,1:4]

typeof(t)

pheatmap (n,show_colnames =F, show_rownames =F)
ac=phe[,c(10,16,17)] 
head(ac)
rownames(ac)=colnames(n)
dim(n)

pheatmap (n,show_colnames =F, 
          show_rownames =F,
          annotation_col = ac)

od=order(ac$Pseudotime)

pheatmap (n[,od],show_colnames =F, 
          show_rownames =F, cluster_cols = F, 
          annotation_col = ac[od,])











