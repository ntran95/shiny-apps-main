seurat_obj <- readRDS("./data/filtered_adj_fpkm_1828_smartseq_integ.RDS")


files <- list.files("./data", pattern = ".RDS", full.names = TRUE)
file_list <- list()

print("Loading Seurat objects...")
for (i in 1:length(files)) {
  file_list[[i]] <- readRDS(files[i])
  DefaultAssay(file_list[[i]]) <- "RNA"
  file_list[[i]] <- ScaleData(file_list[[i]])
}

saveRDS(file_list[[1]], "./data/seurat_obj_input/scaled.filtered_adj_fpkm_1828_smartseq_integ.RDS")
