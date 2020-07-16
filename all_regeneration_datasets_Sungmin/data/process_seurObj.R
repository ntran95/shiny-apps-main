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
  
  }

files <- list.files(".", pattern = "TRIMMED", full.names = TRUE)
file_list <- list()

cell.type <- c("mature-HCs","young-HCs","HC-prog" ,"central-cells", "DV-cells","AP-cells",  
               "amp-SCs","mantle-cells", "Inm","blood", "spi1b-pos","krt17-pos","twist3-pos","tm4sf4-pos")
treatments <- c("homeo" ,"0min" , "30min", "1hr", "3hr","5hr", "10hr")

readSeuratObj <- TRUE
modifySeuratObj <-FALSE

for (i in 1:length(files)) {
  if (readSeuratObj){
  print("Loading Seurat objects...")
  file_list[[i]] <- readRDS(files[i])
  DefaultAssay(file_list[[i]]) <- "RNA"
  }
  if (modifySeuratObj){
  #create new column in meta.data 
    #applied to all-she-pos and neuromast analysis
  if ("cell.type.ident" %in% colnames(file_list[[i]])){
    file_list[[i]]@meta.data$cell.type.ident <- plyr::revalue(
      file_list[[i]]@meta.data$cell.type.ident, c("early-HCs" = "young-HCs"))
  file_list[[i]]@meta.data$cell.type.ident.by.data.set <- factor(paste(file_list[[i]]@meta.data$cell.type.ident,
                                                          file_list[[i]]@meta.data$data.set,sep="_"))
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
    
  saveRDS(file_list[[i]], file = files[i])
  }
}

  #reorder cell.type.ident.by.data.set in all-she-pos analysis & neuromast analysis
for (i in 6:5){
  print(file_list[[i]])
  Idents(file_list[[i]]) <- "cell.type.ident.by.data.set"
  my_levels <- as.vector(t(outer(cell.type, treatments, paste, sep="_"))) 
  file_list[[i]]@active.ident <- factor(file_list[[i]]@active.ident, levels= my_levels)
  file_list[[i]]$cell.type.ident.by.data.set <- factor(file_list[[i]]$cell.type.ident.by.data.set, 
                                                       levels= my_levels)
  file_list[[i]]@active.ident <- droplevels(file_list[[i]]@active.ident)
  file_list[[i]]$cell.type.ident.by.data.set <- droplevels(file_list[[i]]$cell.type.ident.by.data.set)
  Idents(file_list[[i]]) <- "cell.type.ident"
  saveRDS(file_list[[i]], file = files[i])
  
  
}

# ================== test new heatmap ======================
features <- c("atoh1a", "her4.1", "hes2.2", "dld", "sox4a*1", "myclb", "gadd45gb.1",
              "insm1a", "wnt2", "sost", "sfrp1a", "pcna", "mki67", "isl1", "slc1a3a", "glula", "lfng", "cbln20", "ebf3a",
              "znf185", "si:ch211-229d2.5", "si:ch73-261i21.5", "spaca4l", "foxp4", "crip1")

seurat_Obj <- file_list[[6]]

dotplot <- DotPlot(seurat_Obj, features = features,
                   group.by = "cell.type.ident.by.data.set")

dotplot$data$groupIdent <- gsub("(.+?)(\\_.*)", "\\1",dotplot$data$id)
dotplot$data$groupIdent <- factor(dotplot$data$groupIdent,levels=cell.type)

g <- ggplot(dotplot$data, aes(id, features.plot,fill= avg.exp.scaled, width = 1, height = 1)) + 
  geom_tile() +
  scale_fill_distiller(
    palette = "RdYlBu") +
  theme_ipsum()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5,size = 13),
        axis.title.y.right = element_text(size=13),panel.spacing = unit(.35, "lines")) + facet_grid( ~ groupIdent, scales='free_x')



g <- g + labs(title = paste("Selected analysis:",
                            as.character(input$Analysis)), subtitle = "", caption = "") +
  theme(plot.title = element_text(face = "plain", size = 14))





