###CellphoneDB Figure 5

table(thyroid.combined@meta.data$stim)

#T1L
LEFT1_count <-thyroid.combined@assays$RNA@data[,1:3720]
LEFT1_meta <-thyroid.combined@meta.data[1:3720,]
write.table(as.matrix(LEFT1_count), 'LEFT1_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(LEFT1_meta), LEFT1_meta[,'seurat_clusters', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, 'LEFT1_meta.txt', sep='\t', quote=F, row.names=F)

#T2L
LEFT2_count <-thyroid.combined@assays$RNA@data[,3721:7618]
LEFT2_meta <-thyroid.combined@meta.data[3721:7618,]
write.table(as.matrix(LEFT2_count), 'LEFT2_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(LEFT2_meta), LEFT2_meta[,'seurat_clusters', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, 'LEFT2_meta.txt', sep='\t', quote=F, row.names=F)

#T3L
LEFT3_count <-thyroid.combined@assays$RNA@data[,7619:12066]
LEFT3_meta <-thyroid.combined@meta.data[7619:12066,]
write.table(as.matrix(LEFT3_count), 'LEFT3_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(LEFT3_meta), LEFT3_meta[,'seurat_clusters', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, 'LEFT3_meta.txt', sep='\t', quote=F, row.names=F)


#NT
LEFT4_count <-thyroid.combined@assays$RNA@data[,12067:16905]
LEFT4_meta <-thyroid.combined@meta.data[12067:16905,]
write.table(as.matrix(LEFT4_count), 'LEFT4_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(LEFT4_meta), LEFT4_meta[,'seurat_clusters', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, 'LEFT4_meta.txt', sep='\t', quote=F, row.names=F)

#T1R
RIGHT1_count <-thyroid.combined@assays$RNA@data[,12906:20538]
RIGHT1_meta <-thyroid.combined@meta.data[12906:20538,]
write.table(as.matrix(RIGHT1_count), 'RIGHT1_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(RIGHT1_meta), RIGHT1_meta[,'seurat_clusters', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, 'RIGHT1_meta.txt', sep='\t', quote=F, row.names=F)


#T2R
RIGHT2_count <-thyroid.combined@assays$RNA@data[,20539:23826]
RIGHT2_meta <-thyroid.combined@meta.data[20539:23826,]
write.table(as.matrix(RIGHT2_count), 'RIGHT2_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(RIGHT2_meta), RIGHT2_meta[,'seurat_clusters', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, 'RIGHT2_meta.txt', sep='\t', quote=F, row.names=F)


#T3R
RIGHT3_count <-thyroid.combined@assays$RNA@data[,23827:29561]
RIGHT3_meta <-thyroid.combined@meta.data[23827:29561,]
write.table(as.matrix(RIGHT3_count), 'RIGHT3_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(RIGHT3_meta), RIGHT3_meta[,'seurat_clusters', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, 'RIGHT3_meta.txt', sep='\t', quote=F, row.names=F)

#TUMOR
TUMOR_count <-thyroid.combined@assays$RNA@data[,c(1:12066,12906:29561)]
TUMOR_meta <-thyroid.combined@meta.data[c(1:12066,12906:29561),]
write.table(as.matrix(TUMOR_count), 'TUMOR_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(TUMOR_meta), TUMOR_meta[,'seurat_clusters', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, 'TUMOR_meta.txt', sep='\t', quote=F, row.names=F)



#Run following script in Pycharm with python 3.7
#cellphonedb method statistical_analysis  TUMOR_meta.txt  TUMOR_count.txt      --counts-data=gene_name 
#cellphonedb plot dot_plot 
#cellphonedb plot heatmap_plot TUMOR_meta.txt 












