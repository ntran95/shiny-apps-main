library(devtools)
#dev_mode(on=T)
#devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(hrbrthemes)
library(tidyr)

if (TRUE) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  script_name <- "modify_seuratObj"
  
  date <-Sys.Date()
  
  figurePath <- function(filename, format){paste0(script_name, "_figures/", filename)}
  
  if (exists("seurat_obj") == "FALSE"){
    seurat_obj <- readRDS("./TRIMMED_SeurObj_all_LL_cells_regen_v1.2_.RDS")
  }
  
  load(file = "./saved_env_2020-07-10.RData")
  
}

files <- list.files(".", pattern = "TRIMMED", full.names = TRUE)
file_list <- list()

#names(file_list) <- as.character(c("all she-pos. cells"))

print("Loading Seurat objects...")
for (i in 1:length(files)) {
  file_list[[i]] <- readRDS(files[i])
  DefaultAssay(file_list[[i]]) <- "RNA"
  #create new column in meta.data 
  if ("cell.type.ident" %in% colnames(file_list[[i]])){
    file_list[[i]]@meta.data$cell.type.ident <- plyr::revalue(
      file_list[[i]]@meta.data$cell.type.ident, c("early-HCs" = "young-HCs"))
  file_list[[i]]@meta.data$cell.type.ident.by.data.set <- paste(file_list[[i]]@meta.data$cell.type.ident,
                                                                file_list[[i]]@meta.data$data.set,sep="_")
  }
  if ("cell_type" %in% colnames(file_list[[i]]@meta.data)){
    print(files[i])
    file_list[[i]]@meta.data$cell.type.ident.by.data.set <- paste(file_list[[i]]@meta.data$cell_type,
                                                                  file_list[[i]]@meta.data$data.set,sep="_")
    if ("early-HCs" %in% file_list[[i]]$cell_type){
      file_list[[i]]@meta.data$cell_type <- plyr::revalue(
        file_list[[i]]@meta.data$cell_type, c("early-HCs" = "young-HCs"))
    }
  }
  saveRDS(file_list[[i]], file = files[i])
}

