rm(list = ls()) # remove all object in the current enviroment
gc()

# You can restart the R session to avoid any R package conflicts

#load required library
library(tidyverse)

#please set the working directory
#setwd("/your_path/CRC_Data_Fig6")

mydir <- c("data/GSE17536/") #required data files for the GSE17536 dataset
cell_files <- list.files(path = mydir, pattern = "(.*)rda$")
cell_Names <- tools::file_path_sans_ext(list.files(path = mydir, pattern = "(.*)rda$"))

i <- 1
for(cell in cell_files) {
  load(paste0(mydir,"/", cell)) # if you are using windows OS, please use this code "load(paste0(mydir,"\\", cell))"
  assign(cell_Names[i],  pcss)
  i <- i+1
}
rm(pcss)

cell.names <- c("Endothelial", "Epithelial", "GC Cells", "Memory B Cells", "Naive Cells", "Plasma Cells", "CAF-S4", "detox-iCAF", "ecm-myCAF", "IL-iCAF", "TGFB-myCAF", "wound-myCAF", "CD1C-Dendritic Cells", "Granulocytes", "Monocytes", "TAMs", "CD4-CXCL13", "CD4-Resting Cells", "CD4-TH17", "CD4-TREGs", "CD8-GZMA", "CD8-GZMB", "CD8-GZMK", "CD8-MAIT Cells", "CD8-TRM", "ILCs", "NK Cells")
cells <- c("Endothelial", "Epithelial", "GERMINAL_CENTER", "MEMORY_B_CELLS", "NAIVE_CELLS", "PLASMA_CELLS", "CAF_S4", "detox_iCAF", "ecm_myCAF", "IL_iCAF", "TGFB_myCAF", "wound_myCAF", "CD1C_DENDRITIC_CELLS", "GRANULOCYTES", "MONOCYTES", "TUMOR_ASSOSIATED_MACROPHAGES", "CD4_CXCL13", "CD4_Resting_cells", "CD4_TH17", "CD4_Tregs", "CD8_GZMA", "CD8_GZMB", "CD8_GZMK", "CD8_MAIT_Cells", "CD8_TRM", "ILC", "NK_cells")

colors <- c("#dcc134","#008000","#37c8ab","#ff5555")#cms1, cms2, cms3, cms4
colors.toplot <- c(colors, "grey")
names(colors.toplot) <- c(paste0("CMS", 1:4), "unlabeled")

# create plots
cplots <- list()

j <- 1
for (set in cells){
  eset.tmp <- get(set)
  cplots[[j]] <- eset.tmp %>% 
    ggplot(aes(x = PCSS1, y = PCSS2, color = cms_label_SSP)) +ggtitle(cell.names[j])+
    geom_point(aes(alpha = (cms_label_SSP == "not labeled"))) +
    scale_color_manual(values = colors.toplot, name = "CMS subtype") +
    scale_alpha_manual(values = c("TRUE" = 0.5, "FALSE" = 1), guide = F) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold"),axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),legend.direction = "horizontal",
          legend.position = c(0, 1),
          legend.justification=c(0,1),
          legend.background = element_rect(colour="black"))
  j <- j+1
}

# Merge plots in a single figure
library(ggpubr)

path_plots <- c("Figures/GSE17536/")# create a directory to save the plots
if(!file.exists(path_plots)) dir.create(path_plots, recursive=TRUE, showWarnings = FALSE)

figure <- ggarrange(plotlist=cplots, ncol = 6, nrow = 5, common.legend = TRUE,hjust = 0,vjust = 1,legend = "top")
figure <- annotate_figure(figure, 
                           top = text_grob("PCSS1 vs PCSS2 for each cells types", color = "blue", face = "bold", size = 14)
)


f_name <- paste0(path_plots,"PCSS1v2_CellTypes_GSE17536.png")
saveplot::save_png(figure,f_name, width=25, height=20, ppi = 600)

