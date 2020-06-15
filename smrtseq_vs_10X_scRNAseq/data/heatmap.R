library(Seurat)
library(ggplot2)

# ======================================================== function ===================================
obj_integrated@meta.data$data.set <- plyr::revalue(
  obj_integrated@meta.data$data.set, c("homeo-2047" = "homeo-10X-2047",
                                       "homeo-2410-7" = "homeo-10X-2410-7",
                                       "homeo-2410-8" = "homeo-10X-2410-8"))

saveRDS(obj_integrated, "./scaled.filtered_adj_fpkm_1828_smartseq_integ.RDS")

#ids <- as.list(levels(obj_integrated$data.set))

setwd("/home/ntran2/bgmp/shiny-apps-main/smrtseq_vs_10X_scRNAseq/data")

obj_integrated <- readRDS("./scaled.filtered_adj_fpkm_1828_smartseq_integ.RDS")

smartseq_tomatch <- c("1hr-smrtseq", "homeo-smrtseq")

tenX_tomatch <- c("homeo-10X-isl1", "homeo-10X-2047", "homeo-10X-2410-7", "homeo-10X-2410-8")

split_heatmap <- function(seurat_obj, method, tomatch){
  split_obj <- subset(seurat_obj, subset = seq.method == method)
  meta <- split_obj@meta.data
  
  adj.data.set <- as.vector(split_obj@meta.data$data.set)
  
  for (i in 1:length(tomatch)) {
    print(tomatch[[i]])
    meta <- meta%>% mutate(adj.data.set =case_when(str_detect(data.set, 
            paste(tomatch[[i]])) ~ tomatch[[i]],
            TRUE ~ as.vector(split_obj@meta.data$data.set)))
  }
  split_obj@meta.data$adj.data.set <- meta$adj.data.set
  
  return(split_obj)
}

smartseq <- split_heatmap(obj_integrated, method = "smartseq2", tomatch = smartseq_tomatch)

 s <- DoHeatmap(smartseq, features = selected, group.by = "adj.data.set")
 
tenX <- split_heatmap(obj_integrated, method = "10X", tomatch = tenX_tomatch)

t <- DoHeatmap(tenX, features = selected, group.by = "adj.data.set")
 
s + t
