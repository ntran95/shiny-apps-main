
multiGrep2 <- function(toMatch, toSearch, ...) {
  toMatch <- ifelse(grepl("*", toMatch),
    gsub("\\*","\\\\*", toMatch), toMatch <- toMatch)
  
  toMatch <- paste(toMatch, collapse = "|")
  inCommon <- grep(toMatch, toSearch, value = FALSE)
  return(inCommon)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

getLenInput <- function(input) {
  selected <- unlist(strsplit(input, " "))
  
  ifelse(selected %in% com_name,
    selected <- gene_df[com_name %in% selected, 3],
    
    ifelse(selected %in% ens_id,
      selected <- gene_df[ens_id %in% selected, 3],"")
  )
  len <- length(selected)
  return(len)
}

files <- list.files("./data", pattern = "TRIMMED", full.names = TRUE)
file_list <- list()

for (i in 1:length(files)) {
  file_list[[i]] <- readRDS(files[i])
  DefaultAssay(file_list[[i]]) <- "RNA"
}

# !! items to check/change for project (START) !!
file_list <- file_list[c(6,5,1:4)]

# seurat_obj <- file_list[[1]]
print(object.size(file_list), units = "MB")

names(file_list) <- as.character(c(
  "all she-pos. cells", "neuromast cells","AP cells",
  "central cells", "HC progenitors", "mantle cells"))

avg_mtx <- readRDS(paste0("./data/mtx_CLR_nrml_scld_tmpts_",
  "in_cell_type_all_LL_cells_regen_anchored_seurat3_v1.2_.RDS"))
trt_colors <- c("green3", "gold", "darkorange",
  "deeppink", "mediumorchid1", "deepskyblue", "blue")

smpl_genes_sm <- paste0("atoh1a her4.1")
smpl_genes_lg <- paste0("atoh1a her4.1 hes2.2 dld sox4a*1 myclb gadd45gb.1",
" insm1a wnt2 sost sfrp1a pcna mki67 isl1 slc1a3a glula lfng cbln20 ebf3a",
" znf185 si:ch211-229d2.5 si:ch73-261i21.5 spaca4l foxp4 crip1")

app_title <- "Neuromast Regeneration scRNA-seq"

gene_df <- read.table("./data/Danio_Features_unique_Ens91_v2.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# ! items to check/change for project (END) !

ens_id <- gene_df$Gene.stable.ID
com_name <- gene_df$Gene.name.uniq

