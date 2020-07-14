library(devtools)
#dev_mode(on=T)
#devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(hrbrthemes)
library(tidyr)

save_envir <- "FALSE"
load_envir <- "FALSE"
if (TRUE) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  script_name <- "modify_seuratObj"
  
  date <-Sys.Date()
  
  figurePath <- function(filename, format){paste0(script_name, "_figures/", filename)}
  
  if (exists("seurat_obj") == "FALSE"){
    seurat_obj <- readRDS("./TRIMMED_SeurObj_all_LL_cells_regen_v1.2_.RDS")
  }
  if (load_envir ==TRUE){
  load(list.files(path = ".", pattern = "saved_env"))
  }
  if (save_envir == TRUE && file.exists(list.files(path = ".", pattern = ".RData")) == "TRUE"){
    #if (file.exists(list.files(path = ".", pattern = ".RData")) == "TRUE"){
      print("This file exists already:")
      print(list.files(path = ".", pattern = ".RData"))
      system("rm *.RData")
      print("saving updated RData env")
      save.image(file = paste0("saved_env_", date, ".RData"))
    }else{
      print("saving updated RData env")
      save.image(file = paste0("saved_env_", date, ".RData"))
    }
  }

files <- list.files(".", pattern = "TRIMMED", full.names = TRUE)
file_list <- list()

cell.type <- c("mature-HCs","central-cells","young-HCs","amp-SCs","HC-prog" ,"mantle-cells","AP-cells",  
               "DV-cells","Inm","blood", "spi1b-pos","krt17-pos","twist3-pos","tm4sf4-pos")
treatments <- c("homeo" ,"0min" , "30min", "1hr", "3hr","5hr", "10hr")

readSeuratObj <- TRUE
modifySeuratObj <-TRUE
print("Loading Seurat objects...")
for (i in 1:length(files)) {
  if (readSeuratObj == TRUE){
  file_list[[i]] <- readRDS(files[i])
  DefaultAssay(file_list[[i]]) <- "RNA"
  }
  if (modifySeuratObj ==TRUE){
  #create new column in meta.data 
    #applied to all-she-pos and neuromast analysis
  if ("cell.type.ident" %in% colnames(file_list[[i]])){
    file_list[[i]]@meta.data$cell.type.ident <- plyr::revalue(
      file_list[[i]]@meta.data$cell.type.ident, c("early-HCs" = "young-HCs"))
  file_list[[i]]@meta.data$cell.type.ident.by.data.set <- factor(paste(file_list[[i]]@meta.data$cell.type.ident,
                                                                file_list[[i]]@meta.data$data.set,sep="_"))
  #reorder cell.type.ident.by.data.set
  Idents(file_list[[i]]) <- "cell.type.ident.by.data.set"
  my_levels <- as.vector(t(outer(cell.type, treatments, paste, sep="_"))) 
  file_list[[i]]@active.ident <- factor(file_list[[i]]@active.ident, levels= my_levels)

  }
  if ("cell_type" %in% colnames(file_list[[i]]@meta.data)){
    print(files[i])
    file_list[[i]]@meta.data$cell.type.ident.by.data.set <- factor(paste(file_list[[i]]@meta.data$cell_type,
                                                                  file_list[[i]]@meta.data$data.set,sep="_"))
    if ("early-HCs" %in% file_list[[i]]$cell_type){
      file_list[[i]]@meta.data$cell_type <- plyr::revalue(
        file_list[[i]]@meta.data$cell_type, c("early-HCs" = "young-HCs"))
    }
  }
  file_list[[i]]$cell.type.ident.by.data.set <- factor(file_list[[i]]$cell.type.ident.by.data.set)
    
  #saveRDS(file_list[[i]], file = files[i])
  }
}

temp<- file_list[[1]]
Idents(temp) <- "cell.type.ident.by.data.set"
cell.type <- c("mature-HCs","central-cells","young-HCs","amp-SCs","HC-prog" ,"mantle-cells","AP-cells",  
               "DV-cells","Inm","blood", "spi1b-pos","krt17-pos","twist3-pos","tm4sf4-pos")
treatments <- c("homeo" ,"0min" , "30min", "1hr", "3hr","5hr", "10hr")
my_levels <- as.vector(t(outer(cell.type, treatments, paste, sep="_"))) 
temp@active.ident <- factor(temp@active.ident, levels= my_levels)

features <- c("atoh1a", "her4.1", "hes2.2", "dld", "sox4a*1", "myclb", "gadd45gb.1",
              "insm1a", "wnt2", "sost", "sfrp1a", "pcna", "mki67", "isl1", "slc1a3a", "glula", "lfng", "cbln20", "ebf3a",
              "znf185", "si:ch211-229d2.5", "si:ch73-261i21.5", "spaca4l", "foxp4", "crip1")
DotPlot(temp, features = features)
Idents(file_list[[4]])
