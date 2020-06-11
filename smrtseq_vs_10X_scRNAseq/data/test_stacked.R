library(Seurat)
library(ggplot2)

# ======================================== ==============================================


obj_integrated <- readRDS("/home/ntran2/bgmp/shiny-apps-main/smrtseq_vs_10X_scRNAseq/data/subsampled.30.scaled.filtered_adj_fpkm_1828_smartseq_integ.rds")

ids <- as.list(levels(obj_integrated$data.set))


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


obj_trt_list <- list()[1:length(ids)]
for (i in 1:length(ids)) {
  print(ids[[i]])
  obj_trt_list[[i]] <- obj_integrated[,obj_integrated[["data.set"]] == ids[[i]]]
}
#})


stacked_violin_plot <- function(goi, obj_trt_list){
  for (i in 1:length(ids)) {
    vln_obj <- VlnPlot(
      obj_trt_list[[i]], features = goi, pt.size = 0) +
      xlab("") + ylab(ids[i]) + ggtitle("") +
      theme(legend.position = "none", axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = rel(1), angle = 0),
            axis.text.y = element_text(size = rel(1)),
            plot.margin = unit(c(-0.75, 0.5, -0.75, 0.5), "cm"))
    trt_plot_list[[i]] <- vln_obj
  }
  
  trt_plot_list[[length(trt_plot_list)]]<- trt_plot_list[[length(trt_plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  # change the y-axis tick to only max value
  ymaxs <- purrr::map_dbl(trt_plot_list, extract_max)
  trt_plot_list <- purrr::map2(trt_plot_list, ymaxs, function(x, y) x +
                                 scale_y_continuous(breaks = c(y)) + expand_limits(y = y))
  grid_obj <- cowplot::plot_grid(plotlist = trt_plot_list,
                                 nrow = length(ids), ncol = 1, axis = "l", align = "hv") +
    theme(plot.margin = margin(1.0, 1.0, 1.0, 1.0, unit = "in")) +
    ggtitle(goi) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  

  return(grid_obj)
  
}

selected <- c("atoh1a","wnt2")
grid_obj <- list()[1:length(selected)]
for (i in 1:length(goi)){
  print(goi[[i]])
  grid_obj[[i]] <- stacked_violin_plot(goi = selected[[i]], obj_trt_list = obj_trt_list)
  
}
names(grid_obj) <- selected

final_grid <- cowplot::plot_grid(plotlist = grid_obj, nrow = length(grid_obj))
final_grid
