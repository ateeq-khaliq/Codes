save_markers_to_excel_pages_by_cluster <- function(marker_txt_file_path, marker_excel_path) {
  
library(dplyr)
library(writexl)
marker_Genes <- read.delim(marker_txt_file_path)
cluster_count <- table(marker_Genes$cluster)
for (i in 0:(length(cluster_count)+1)) {
  print(i)
  j <- paste0("cluster_",as.character(i)) 
  marker_Genes_filtered <- as.data.frame(marker_Genes %>% select(everything()) %>% filter(p_val <= 0.05 & avg_log2FC>0.5 & cluster==i) %>% arrange(desc(avg_log2FC))) # 
  assign(j, marker_Genes_filtered)
}

data_list <- list()

for (i in 0:(length(cluster_count)+1)) {
  nData <- paste0("cluster_",i)
  vData <- get(paste0("cluster_",i))
  data_list[[nData]] <- vData
}

system.time(write_xlsx(data_list, marker_excel_path))
print(paste0("Excel files saved to : ",marker_excel_path))
}