setwd("bgmp/shiny-apps-main/smrtseq_vs_10X_scRNAseq/")
obj_integrated <- readRDS("./data/seurat_obj_input/scaled.filtered_adj_fpkm_1828_smartseq_integ.rds")

# ==== Stacked violin plot
ids <- c("1hr-smrtseq", "homeo-10X-isl1", "homeo-2047", "homeo-2410-7", "homeo-2410-8", "homeo-smrtseq")
ids <- as.list(levels(obj_integrated$data.set))
ids <- as.list(levels(obj_integrated$cell.type.ident))
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
  obj_trt_list[[i]] <- obj_integrated[,obj_integrated[["cell.type.ident"]] == ids[[i]]]
}


goi <- c("atoh1a")
trt_plot_list <- list()[1:length(ids)]
names(trt_plot_list) <- ids
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
#png(figurePath(paste0(goi,"_stacked_vln.png")),
 #   width = 16, height = 8, units = "in", res = 300)
print(grid_obj)
#dev.off()


# =================== other way ========================

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0.2, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, group.by = "cell.type.ident" )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0.2, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

selected <- c("atoh1a", "her4.1")
StackedVlnPlot(obj = obj_integrated, features = selected)
