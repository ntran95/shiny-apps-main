library(Seurat)
library(dplyr)

setwd("./bgmp/shiny-apps-main/smrtseq_vs_10X_scRNAseq/data/seurat_obj_input/")
seurat_obj <- readRDS("scaled.filtered_adj_fpkm_1828_smartseq_integ.rds")

############# downsample by 75% per cluster identity ############# 
set.seed(1)
sampled_seurat_obj <- subset(seurat_obj, cells = sample(Cells(seurat_obj), size = round(length(colnames(seurat_obj)) * .75)))
table(Idents(sampled_seurat_obj))
table(Idents(seurat_obj))

print(object.size(sampled_seurat_obj), units = "MB")

saveRDS(sampled_seurat_obj, file = "../subsampled.75.scaled.filtered_adj_fpkm_1828_smartseq_integ.rds")

############# downsample by 50% per cluster identity ############# 
set.seed(1)
sampled_seurat_obj <- subset(seurat_obj, cells = sample(Cells(seurat_obj), size = round(length(colnames(seurat_obj)) * .50)))
table(Idents(sampled_seurat_obj))
table(Idents(seurat_obj))

print(object.size(sampled_seurat_obj), units = "MB")

saveRDS(sampled_seurat_obj, file = "../subsampled.50.scaled.filtered_adj_fpkm_1828_smartseq_integ.rds")


############# downsample by 30% per cluster identity ############# 
set.seed(1)
sampled_seurat_obj <- subset(seurat_obj, cells = sample(Cells(seurat_obj), size = round(length(colnames(seurat_obj)) * .30)))
table(Idents(sampled_seurat_obj))
table(Idents(seurat_obj))
DimPlot(sampled_seurat_obj)

print(object.size(sampled_seurat_obj), units = "MB")

saveRDS(sampled_seurat_obj, file = "../subsampled.30.scaled.filtered_adj_fpkm_1828_smartseq_integ.rds")



downsample_list <- c(.75, .50, .30)
sampled_seurat_obj <- list()
set.seed(1)
for (i in 1:length(downsample_list)){
  print(downsample_list[i])
  print(sampled_seurat_obj[i])
  #set.seed(1)
  #sampled_seurat_obj[[i]] <- subset(seurat_obj, cells = sample(Cells(seurat_obj), size = round(length(colnames(seurat_obj)) * as.numeric(downsample_list[i]))))
  #table(Idents(sampled_seurat_obj))
  #table(Idents(seurat_obj))
  
  #print(object.size(sampled_seurat_obj), units = "MB")
}
sampled_seurat_obj