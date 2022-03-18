
#Subset cells 
Subset_Cells <- subset(thyroid.combined, idents = "Follicular_cells")

#scaledate
Subset_Cells <- ScaleData(Subset_Cells, verbose = FALSE)
#FindVariableFeatures
Subset_Cells <- FindVariableFeatures(Subset_Cells, selection.method = "vst", nfeatures = 2000)
#
Subset_Cells <- RunPCA(Subset_Cells, npcs = 50, verbose = FALSE)
Subset_Cells <- FindNeighbors(Subset_Cells, reduction = "pca", dims = 1:10)
Subset_Cells <- FindClusters(Subset_Cells, resolution = 0.5)
Subset_Cells <- RunUMAP(Subset_Cells, reduction = "pca", dims = 1:5)
Subset_Cells <- RunTSNE(Subset_Cells, reduction = "pca", dims = 1:5)


#UMAP for Figure 2A and Figure 2B
DimPlot(Subset_Cells, reduction = "umap",group.by = "stim", cols = c("red", "red", "red", "blue", "red","red", "red"))
DimPlot(Subset_Cells, reduction = "umap")


#cell type frequency for Figure 2C
table(Subset_Cells$stim,Subset_Cells$seurat_clusters)



#DEGs in follicular cells between tumor and non tumor (Supplementary file 2)
b.interferon.response <- FindMarkers(Subset_Cells, group.by = "stim", ident.1 = c("T1L","T2L","T3L","T1R","T2R","T3R"), 
                                     ident.2 = c("NT"), verbose = TRUE)
write.csv(b.interferon.response,"Follicularcells_TumorvsNonTumor.csv")


#KEGG pathways enriched in malignant follicular cells in Figure 2D
sig_dge.celltype <- subset(b.interferon.response, p_val_adj<0.01&abs(avg_log2FC)>1)

genelist <- bitr(row.names(sig_dge.celltype), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
p1 <- barplot(ekegg, showCategory=50)
p2 <- dotplot(ekegg, showCategory=50, orderBy = "x")
p1
p2

#KEGG_allTumorvsNontumor1 after selecting pathways of interest as shown Figure 2D
ekegg1 <- read.csv("KEGG_allTumorvsNontumor1.csv", row.names = 1)
ggplot(ekegg1,aes(x = GeneRatio,y = Description),showCategory=50)+
  geom_point(aes(color = p.adjust,
                 size = Count))+
  scale_color_gradient(low = "red", high = "blue")+
  xlab("GeneRatio")+
  theme_bw()


#Enhandced Vocalno plot in Figure 2E

#install EnhancedVolcano
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)

volcano<-subset(b.interferon.response,select = c(p_val_adj,avg_log2FC))
threshold<-as.factor(abs(volcano$avg_log2FC)>1&volcano$p_val_adj<0.01)


EnhancedVolcano(volcano,lab = rownames(volcano),x = "avg_log2FC",y = "p_val_adj",
                #selectLab = c("CXCL10","GDF15","CXCL8","CXCL9","CXCL3","CXCL1","IL32","CCL4L2","CCL20","CCL3","CCL3L1","LTB","IL1A","TNFRSF12A","IL18","CCL4","CCL5"),
                selectLab = c("CXCL10","GDF15","CXCL8","CXCL9","CXCL3","IL32","CCL4","CCL4L2","CCL20","CCL3","IL1A","CCL5","TNFRSF12A","CXCL1","IL18"),
                xlim = c(-100, 400),
                ylim = c(0, 120),
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                border = "full",
                #legendPosition = "bottom",
                borderWidth = 1.5,
                borderColour = "blue",
                drawConnectors = TRUE,
                widthConnectors = 1,
                #colConnectors = c('black','black','orange','orange'),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                boxedLabels = TRUE)





####GSVA for ploting 2H-2J
library(GSVA)
library(GSEABase)
library(limma)

#
keggSet = getGmt("/c2.cp.kegg.v7.4.symbols.gmt")   #(browse 186 gene sets)

#subset
Subset_Cells <- subset(thyroid.combined, idents = "Follicular_cells")
Subset_Cells <- subset(Subset_Cells, subset= (stim == "T1L" | stim == "T1R") ) 
table(Subset_Cells$stim)

expr <- as.data.frame(Subset_Cells@assays$RNA@data)
write.csv(expr,"expr.csv")


keggEs = gsva(expr=as.matrix(expr), gset.idx.list = keggSet, kcdf="Gaussian", parallel.sz=1, method="gsva")  #适用于标准化后的log2FPKM等数据
write.csv(keggEs,"keggEs.csv")


exprSet <- keggEs

meta <- Subset_Cells@meta.data[,c("stim")]  
group <- factor(meta,levels = c("T1L","T1R"),ordered = F)  

design <- model.matrix(~group)

colnames(design) <- levels(group)

fit <- lmFit(exprSet,design)
fit2 <- eBayes(fit)

allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf)
write.csv(allDiff,"allDiffKEGG_T1LvsT1R.csv")

#n=12
up <- c("KEGG_GLYCOSAMINOGLYCAN_DEGRADATION","KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM","KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM","KEGG_OXIDATIVE_PHOSPHORYLATION","KEGG_CELL_CYCLE","KEGG_APOPTOSIS","KEGG_ERBB_SIGNALING_PATHWAY","KEGG_CELL_ADHESION_MOLECULES_CAMS","KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","KEGG_THYROID_CANCER","KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY","KEGG_MAPK_SIGNALING_PATHWAY")
#n=11
down <- c("KEGG_HEDGEHOG_SIGNALING_PATHWAY","KEGG_ECM_RECEPTOR_INTERACTION","KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM","KEGG_LINOLEIC_ACID_METABOLISM","KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION","KEGG_STEROID_HORMONE_BIOSYNTHESIS","KEGG_RETINOL_METABOLISM","KEGG_GLYCOSAMINOGLYCAN_BIOSYNTHESIS_HEPARAN_SULFATE","KEGG_PPAR_SIGNALING_PATHWAY","KEGG_PRIMARY_IMMUNODEFICIENCY",
          "KEGG_PHOSPHATIDYLINOSITOL_SIGNALING_SYSTEM")
TEST <- c(up,down)

top10 = keggEs[rownames(keggEs) %in% TEST,]

dim(keggEs)
view(keggEs)
top10_1 = top10[,1:626]
top10_2 = top10[,627:1014]
top10_1$mean = rowMeans(top10_1)
top10_2$mean = rowMeans(top10_2)
data = data.frame(Normal = top10_1$mean,Treat = top10_2$mean)
rownames(data) = rownames(top10)

library(pheatmap)
annotation_col = data.frame(group = c(rep("T1L",1),rep("T1R",1)))
rownames(annotation_col) = colnames(data)

view(data)
data$group <- c(1,2,2,2,1,1,2,1,2,2,2,1,1,2,2,1,1,2,2,1,1,1,1)

data1 <- subset(data,data$group==1)
data1 <- data1[order(data1$Treat, decreasing = T),]

data2 <- subset(data,data$group==2)
data2 <- data2[order(data2$Normal, decreasing = T),]

data3 <- bind_rows(data1,data2)

data4 <- data3[,-3]


pheatmap(data4,
         annotation_col = annotation_col,
         cluster_cols = F,
         cluster_row = FALSE,
         show_colnames = F)

##similar method to plot Figure 2I and 2J




####InferCNV 


sce=thyroid.combined


table(Idents(sce))
table(sce@meta.data$seurat_clusters) 
table(sce@meta.data$orig.ident) 

dat=GetAssayData(sce,
                 slot='counts',assay='RNA')
dat[1:4,1:4]
dim(dat)
groupinfo=data.frame(v1=colnames(dat),
                     v2= Idents(sce) )
head(groupinfo)

BiocManager::install("AnnoProbe")


library(broom)
library(AnnoProbe)
geneInfor=annoGene(rownames(dat),
                   "SYMBOL",'human')
colnames(geneInfor)
head(geneInfor)
geneInfor=geneInfor[with(geneInfor,
                         order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
dim(geneInfor)


dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
head(groupinfo)
dat[1:4,1:4]
table(groupinfo$v2)
head(groupinfo)
dim(groupinfo)


expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)

groupFiles='groupFiles.txt'

groupinfo[,1]=gsub("-",".",groupinfo[,1])  ##to replace “-” with “.”
write.table(groupinfo,file = groupFiles,sep = '\t',
            quote = F,col.names = F,row.names = F)

head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',
            quote = F,col.names = F,row.names = F)


options(stringsAsFactors = F)
library(Seurat)
library(gplots)
library(ggplot2)

expFile='expFile.txt' 
groupFiles='groupFiles.txt'  
geneFile='geneFile.txt'

library(edgeR)
library(rjags)
library(infercnv)
packageVersion("Signac")
install.packages("edgeR")



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR", force=TRUE)
BiocManager::install("infercnv", force=TRUE)


infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=c("T_cells","B_cells") )

# cutoff=0.1 works well for 10x Genomics to generate the Figure 2K
infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir= "infercnv_output",  # dir is auto-created for storing outputs
                              cluster_by_groups=TRUE,   # cluster 
                              #cluster_references = F,
                              HMM=FALSE,
                              denoise=TRUE,
                              hclust_method="ward.D2", plot_steps=F)

# for subsequent analysis and generating Figure 2L and 2M and Supplementary Figure S4
infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir= "infercnv_output",  # dir is auto-created for storing outputs
                              cluster_by_groups=F,   # cluster 若为FALSE，and cluster_references = F，no error occured
                              cluster_references = F,
                              HMM=FALSE,
                              denoise=TRUE,
                              hclust_method="ward.D2", plot_steps=F)



library(phylogram)
library(gridExtra)
library(grid)
library(dendextend)
library(ggthemes)
library(miscTools)

#import inferCNV dendrogram
infercnv.dend <- read.dendrogram(file = "infercnv_output/infercnv.observations_dendrogram.txt")


#cut tree
infercnv.labels <- cutree(infercnv.dend, k = 6, order_clusters_as_data =FALSE)
table(infercnv.labels)

#color labels
the_bars <- as.data.frame(tableau_color_pal("Tableau 20")(20)[infercnv.labels])
colnames(the_bars) <- "inferCNV_tree"

the_bars$inferCNV_tree <- as.character(the_bars$inferCNV_tree)

infercnv.dend %>% set("labels", rep("",nobs(infercnv.dend))) %>% plot(main="inferCNV dendrogram") %>% 
  colored_bars(colors = as.data.frame(the_bars), dend = infercnv.dend, sort_by_labels_order = FALSE, add = T, y_scale=100 , y_shift = 0)  

infercnv.labels = as.data.frame(infercnv.labels)
groupFiles = 'groupFiles.txt'
meta=read.table(groupFiles, sep = '\t')
infercnv.labels$V1=rownames(infercnv.labels)
meta=merge(meta,infercnv.labels,by = 'V1')
table(meta[,2:3])

if( ! file.exists("cnv_scores.csv")){
  tmp=read.table("inferCNV_output/infercnv.references.txt", header = T)
  down=mean(rowMeans(tmp)) - 2 * mean(apply (tmp,1,sd))
  up=mean(rowMeans(tmp)) + 2 * mean(apply (tmp,1,sd))
  oneCopy=up-down
  oneCopy
  a1=down-2*oneCopy
  a2=down-1*oneCopy
  down;up
  a3=up+1*oneCopy
  a4=up+2*oneCopy
  
  cnv_table <- read.table("inferCNV_output/infercnv.observations.txt", header = T)
  
  cnv_score_table <- as.matrix(cnv_table)
  
  cnv_score_mat <- as.matrix(cnv_table)
  
  #scoring
  cnv_score_table[cnv_score_mat >0 & cnv_score_mat < a2] <- "A"  #complete loss. 2pts
  cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B"  #loss of one copy. 1pts
  cnv_score_table[cnv_score_mat >= down & cnv_score_mat < up] <- "C"  #Neutral . opts
  cnv_score_table[cnv_score_mat >= up & cnv_score_mat <= a3] <- "D"  #addition of one copy. 1pts
  cnv_score_table[cnv_score_mat > a3 & cnv_score_mat <= a4] <- "E"  #addition of two copies. 2pts
  cnv_score_table[cnv_score_mat > a4] <- "F"  #addition of more than two copies. 2pts
  
  #check
  table(cnv_score_table[,1])
  #replace with score
  cnv_score_table_pts <- cnv_table
  rm(cnv_score_mat)
  #
  cnv_score_table_pts[cnv_score_table == "A"] <-2
  cnv_score_table_pts[cnv_score_table == "B"] <-1
  cnv_score_table_pts[cnv_score_table == "C"] <-0
  cnv_score_table_pts[cnv_score_table == "D"] <-1
  cnv_score_table_pts[cnv_score_table == "E"] <-2
  cnv_score_table_pts[cnv_score_table == "F"] <-2
  
  #scores are stored in "cnv_score_table_pts", use colSums to add up scores for each cell and store as 
  cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
  colnames(cell_scores_CNV) <- "cnv_score"
  head(cell_scores_CNV)
  write.csv(x=cell_scores_CNV, file="cnv_scores.csv")
  
  
}


cell_scores_CNV = read.csv('cnv_scores.csv', row.names = 1)
head(cell_scores_CNV)


sce=thyroid.combined
phe=sce@meta.data
head(phe)
phe$celltype=Idents(sce)
head(rownames(phe))
head(rownames(cell_scores_CNV))

#rowname(phe)=paste0('X',rownames(phe))
rownames(phe) = gsub("-",".",rownames(phe))

head(rownames(phe))
head(rownames(cell_scores_CNV))

head(rownames(phe))
phe=phe[rownames(phe) %in% rownames(cell_scores_CNV),]
identical(rownames(phe), rownames(cell_scores_CNV))

dim(phe)

infercnv.labels <- cutree(infercnv.dend, k=6, order_clusters_as_data = FALSE)
phe$inferCNV = infercnv.labels[match(rownames(phe),names(infercnv.labels))]

phe$cnv_scores = cell_scores_CNV[rownames(phe),]

table(phe$celltype, phe$inferCNV)
head(rownames(phe))
dim(phe)

library(ggpubr)

p1=ggboxplot(phe,'celltype','cnv_scores', fill ="celltype")
p1 = p1 + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
p2=ggboxplot(phe,'inferCNV','cnv_scores',fill="inferCNV")

phe1=phe[phe$celltype=="Follicular_cells",]
head(phe1)
p3=ggboxplot(phe1,'stim','cnv_scores',fill="stim",order = c("T1L", "T1R", "T2L", "T2R", "T3L", "T3R", "NT"))
p3 = p3 + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))


library(patchwork)

####Figure 2L
p1
p2

####Figure 2M
p3

ggsave(filename = 'anno_CNVscore.pdf')
table(phe$celltype, phe$inferCNV)

write.csv(table(phe$celltype, phe$inferCNV),"celltype_inferCNV.csv")



phe2=phe[phe$celltype=="Pericyte",]
head(phe2)

phe3=phe[phe$celltype=="Endothelial_cells",]
head(phe3)

phe4=phe[phe$celltype=="Myeloid_cells",]
head(phe4)

phe5=phe[phe$celltype=="Fibroblast",]
head(phe5)

phe6=phe[phe$celltype=="Mast_cells",]
head(phe6)

p4=ggboxplot(phe2,'stim','cnv_scores',fill="stim",order = c("T1L", "T1R", "T2L", "T2R", "T3L", "T3R", "NT"))
p4 = p4 + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

p5=ggboxplot(phe3,'stim','cnv_scores',fill="stim",order = c("T1L", "T1R", "T2L", "T2R", "T3L", "T3R", "NT"))
p5 = p5 + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

p6=ggboxplot(phe4,'stim','cnv_scores',fill="stim",order = c("T1L", "T1R", "T2L", "T2R", "T3L", "T3R", "NT"))
p6 = p6 + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

p7=ggboxplot(phe5,'stim','cnv_scores',fill="stim",order = c("T1L", "T1R", "T2L", "T2R", "T3L", "T3R", "NT"))
p7 = p7 + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

p9=ggboxplot(phe7,'stim','cnv_scores',fill="stim",order = c("T1L", "T1R", "T2L", "T2R", "T3L", "T3R", "NT"))
p8 = p8 + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

##Supplementary Figure S4
p4+p5+p6+p7+p8













