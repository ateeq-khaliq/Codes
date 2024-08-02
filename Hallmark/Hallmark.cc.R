suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dittoSeq))
suppressPackageStartupMessages(library(ggplot2))

library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(GSEABase)

print("reading pdac")
setwd("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/hallmark/minuscc")

H.gsea <- getGmt("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/correlation/h.all.v2023.1.Hs.symbols.gmt", collectionType=BroadCollection(category="h"), geneIdType=SymbolIdentifier())

pdac <- readRDS("/Users/akhaliq/Desktop/data_for analysis/st_data/Pdac_updated.rds")

Idents(pdac) <- pdac$cc_ischia_10

a <- subset(pdac, idents=c("CC1","CC2","CC3","CC4","CC5","CC7","CC10")) 

print(a)
rm(pdac)

# CC6,9 & 8

#a <- readRDS("Pdac_allres_most_updated.rds")
print("Starting Es for Minus CCs")

DefaultAssay(a) <- "Spatial"
a <- RenameAssays(object = a, Spatial = 'RNA')
Es <- enrichIt(obj = a, 
               gene.sets = H.gsea, 
               groups = 1000, cores = 8)

saveRDS(Es,"ES_minusCCs.rds")
print("Es_minusCCs Done !!! ")


print("Starting for all PDAC!")

print("Reading pdac")
setwd("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/hallmark/allcc")

H.gsea <- getGmt("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/correlation/h.all.v2023.1.Hs.symbols.gmt", collectionType=BroadCollection(category="h"), geneIdType=SymbolIdentifier())

pdac <- readRDS("/Users/akhaliq/Desktop/data_for analysis/st_data/Pdac_updated.rds")

Idents(pdac) <- pdac$cc_ischia_10

#a <- subset(pdac, idents=c("CC1","CC2","CC3","CC4","CC5","CC7","CC10")) 
a=pdac
print(a)
rm(pdac)

# CC6,9 & 8

#a <- readRDS("Pdac_allres_most_updated.rds")
print("Starting Es for All CCs")

DefaultAssay(a) <- "Spatial"

a <- RenameAssays(object = a, Spatial = 'RNA')
Es <- enrichIt(obj = a, gene.sets = H.gsea, groups = 1000, cores = 8)

saveRDS(Es,"ES_allcc.rds")
print("Es_allcc Done !!! ")


###
Es <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/hallmark/allcc/ES_allcc.rds")
library(stringr)

colnames(Es) <- str_remove(colnames(Es), "^HALLMARK_")

tibble_data <- Es %>% 
  rownames_to_column(var = "condition") %>%
  pivot_longer(cols = -condition, names_to = "source", values_to = "score") %>%
  arrange(source, condition) %>%
  select(source, condition, score)

a = pdac

a[['mgsigdb']] <-tibble_data %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

DefaultAssay(a) <- "mgsigdb"
a <- ScaleData(a)

a@assays$mgsigdb@data <- a@assays$mgsigdb@scale.data
Idents(a) <- a$cc_ischia_10


data <- a

rm(a)

df <- t(as.matrix(data@assays$mgsigdb@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))
#`summarise()` has grouped output by 'cluster'. You can override using the `.groups` argument.
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

####  
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

pdf("Hallmark_allcc1.pdf",width = 25, height = 10)
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()
####

# Load required libraries
library(pheatmap)

# Define the number of colors in the palette
palette_length <- 100

# Generate a color palette using colorRampPalette
my_color <- colorRampPalette(c("blue", "white", "red"))(palette_length)

# Create legend labels to match the breaks
legend_labels <- c(
  paste("<= ", round(my_breaks[1:ceiling(length(my_breaks)/2)], 2)),
  paste(">= ", round(my_breaks[ceiling(length(my_breaks)/2) + 1:length(my_breaks)], 2))
)

# Set up a PDF graphics device for high-quality output
pdf("hallmark_minuscc.pdf", width = 20, height = 10)

# Create a heatmap plot using pheatmap
pheatmap(
  top_acts_mat,
  border_color = NA,
  color = my_color,
  breaks = my_breaks,
  cellwidth = 20,
  cellheight = 20,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "Heatmap For Hallmark Genesets",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  show_colnames = TRUE,
  show_rownames = TRUE,
  angle_col = 90
)

# End the PDF graphics device
dev.off()


