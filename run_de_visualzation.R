## Visualize differentially expressed genes

# Receive variables from the user (Script Only) -----------------------------------------
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  stop("please provide 5 args: <path_to_obj> <path_to_de> <path_to_plots> <K> <type_of_de>")
  
} else {
  path_to_obj <- args[1] # Obj
  path_to_de <- args[2] # DE 
  path_to_plots <- args[3]
  K <- args[4] # K value
  LFC_type <- args[5] # the type of the calculation of the LogFoldChange (vsnull, k, le) more info in the "FastTopic packege")
}
# Loading libraries  ------------------------------------------------------
load_libraries <- function() { 
  library(Matrix)
  library(fastTopics)
  library(ggplot2)
  library(cowplot)
  library(RColorBrewer)
  library(hexbin)
  library(viridis)
  library(Seurat)
  library(SeuratDisk)
  library(SeuratData)
  library(magrittr)
  library(dplyr)
  library(gridExtra)
  library(ggplot2)
  source("utils.R")
  library(grid)
  library(viridisLite)


}

# Global variables  --------------------------------------------------------
font <- "Arial"
additional_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                       "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                       "#cab2d6", "#6a3d9a")
# Optional Part (see doc) -------------------------------------------------
# This is temp processes. While running the flow few rows (0 rows) of the
#obj@meta.data we removed as part of the flow, this rows (0 rows)
#where NOT remove from the Counts Matrix and therefore need to be removed
remove_zeros_from_counts_mat <- function(obj)
{
  counts <- obj@assays$RNA@counts
  row_sum <- rowSums(counts)
  nonzero_rows <- row_sum != 0
  counts <- counts[nonzero_rows,]
  counts <- t(counts)
  return (counts)
}
# VolcanoPlot -------------------------------------------------------------


#' Generate Volcano Plots for Differential Expression Analysis
#'
#' This function generates a set of volcano plots for differential expression analysis using the DESeq2 package.
#' Each volcano plot represents the results of differential expression analysis for a specific condition or group.
#'
#' @details The function calculates and visualizes differential expression (DE) results for multiple conditions or groups.
#' It iterates through each condition and generates a volcano plot showing log-fold change (LFC) versus statistical significance (-log10 p-value) for each gene.
#' 
#' @return None
#' @export
#' 
#' @examples
#' # Generate volcano plots for differential expression analysis
#' generate_VolcanoPlot()
#' 
#' @import 
#' 
generate_VolcanoPlot <- function(de, path_to_plot, counts, K=10, LFC_type="le") {
  print("Generating Volcano Plots...")
  
  de_VolcanoPlot <- list()
  for (i in 1:K) {
    de_VolcanoPlot[[i]] <- volcano_plot(de,
                                        k = i,
                                        labels = colnames(counts)) 
  }
  combined_plot <- arrangeGrob(grobs = de_VolcanoPlot,
                               ncol = 2,
                               nrow = 5)
  ggsave(combined_plot, 
         file = paste0(path_to_plot,
                       "VolcanoPlot_DE_LFC_", 
                       LFC_type, 
                       ".pdf"),
         width = 12,
         height = 20)
  print("Done with Volcano Plots")
}



# Find Markers ------------------------------------------------------------
significant_genes <- list()
genes <- rownames(de$lfsr)

for (i in 1:K) {
 markers <- intersect(genes[de$lfsr[, i] < 0.001], genes[de$est[, i] > 0])
 significant_genes <- append(significant_genes, list(markers))
}

# open the markers file
markers <- read.csv(path_to_markers, header = TRUE, row.names = 1)
markers$representative_topic <- c(1, 2,3,4,4,5,5,5)
#markers



# Main Flow Run -----------------------------------------------------------
main <- function(){
  rm(list = ls())
  
  load_libraries()
  ## local run - REMOVE before running in the cluster
  path_to_obj <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/objects/obj_and_topics.h5seurat"
  path_to_de <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/objects/de_vsnull"
  path_to_plot <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/plots/"
  path_to_markers <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/cell_marker_24_8_23.csv"
  K <- 10
  LFC_type <- "VS_NULL"
  COMPLEX_CASE <- T
  
  
  # genrate new directory for the plots
  folder_name <- paste0("Plots_for_DE_K_", K, "_LFC_Type_", LFC_type)
  dir.create(file.path(path_to_plot, folder_name))
             #path_to_plot, folder_name)
  
  path_to_plot <- paste0(path_to_plot, folder_name, "/")
  #path_to_plot + folder_name
  # load objects 
  print("Loading the objects - START")
  obj <- LoadH5Seurat(path_to_obj)
  de <-readRDS(path_to_de)
  print("Loading the objects - FINISH")
  
  
  # Optional see documentation 
  counts <- remove_zeros_from_counts_mat(obj)
  
  generate_VolcanoPlot(de = de, 
                       path_to_plot = path_to_plot, 
                       counts=counts, 
                       K=as.integer(K), 
                       LFC_type =LFC_type)
  
  generate_VlnPlot(obj = obj,
                   path_to_plot = path_to_plot,
                   LFC_type = LFC_type,
                   K=as.integer(K))
  
  generate_vln_plot_complex(obj = obj,
                            path_to_plot = path_to_plot, 
                            LFC_type = LFC_type,
                            K=as.integer(K))
  
  print(paste0("The Run is Finished. The Plots are saved at: ",path_to_plot))
  
}

main()

##### sandbox  ####





load_libraries()

# from matrix de$postmean, get the top 100 genes (Accoridng to the val in col5)
# and save them in a list
total <- 100
topic <- 3
pattern <- "^AL"

top_100 <- de$postmean[order(de$postmean[,topic], decreasing = TRUE),][1:total,]
top_100_topic3 <- de$postmean[order(de$postmean[,3], decreasing = TRUE),][1:total,]
top_100_topic5 <- de$postmean[order(de$postmean[,3], decreasing = TRUE),][1:total,]
start_with_AC_topic_3 <- rownames(top_100_topic3)[grepl(pattern, rownames(top_100_topic3))]
start_with_AC_topic_5 <- rownames(top_100_topic5)[grepl(pattern, rownames(top_100_topic5))]
# removing all the col but k5
inter <- intersect(start_with_AC_topic_3, start_with_AC_topic_5)
union <- union(start_with_AC_topic_3, start_with_AC_topic_5)
length(inter)
length(union)

top_100 <- data.frame(top_100[,topic])
colnames(top_100) <- c("k5")

# count how many col in top_100 staring with "AL"

start_with_AC_ <-
# list of all the row names in top_100 staring with "AC"




prec_start_with_AL <-(sum(grepl("^AL", rownames(top_100))))/total
prec_start_with_LI <-sum(grepl("^LI", rownames(top_100)))/total
prec_start_with_ST <-sum(grepl("^ST", rownames(top_100)))/total
prec_start_with_TT <-sum(grepl("^TT", rownames(top_100)))/total
prec_start_with_LIN <-sum(grepl("^LIN", rownames(top_100)))/total
prec_start_with_AP <-sum(grepl("^AP", rownames(top_100)))/total
prec_start_with_BEX <-sum(grepl("^BEX", rownames(top_100)))/total

sum(prec_start_with_AP, prec_start_with_LIN, prec_start_with_AL, prec_start_with_AC, prec_start_with_LI, prec_start_with_ST, prec_start_with_TT)/2

K <- 10
significant_genes <- list()
genes <- rownames(de$lfsr)
markers <- intersect(genes[de$lfsr[, 5] < 0.001], genes[de$postmean[, 5] > 3])
k_5 <- fit$F[,markers]

# significant_genes <- append(significant_genes, list(markers))



# (Still working on it) Find Markers of Clusters  -----------------------------------------------

obj_markers <- FindAllMarkers(object = obj,
                              only.pos = TRUE,
                              min.pct=.25)


obj_markers_with_cell_type <- obj_markers %>% 
  #group_by(cluster) %>% 
  mutate(celltype = case_when(
  cluster %in% c(0, 1) ~ "Mature_oligodendrocytes",
  cluster == 2 ~ "Astrocytes",
  cluster %in% c(3, 5, 7, 8, 9) ~ "Neurons",
  cluster == 4 ~ "Microglia",
  cluster == 6 ~ "OPC",
  cluster == 10 ~ "vascular",
  cluster == 11 ~ "Astrocytes"#,
  #TRUE ~ NA_character_
))


obj_markers_with_cell_type1 <- obj_markers_with_cell_type %>% group_by(cluster) %>% slice_max(n = 100, order_by = avg_log2FC)
obj_markers_with_cell_type_top10 <- obj_markers_with_cell_type1 %>%
  group_by(cluster) %>% 
  slice_max(n = 10, order_by = avg_log2FC)





################################################################################
## visualize DE
################################################################################

## Topic Proportion - Matrix L (using UMAP or t-SNE)
de_le <-readRDS("de_le")
de_vsnull <- readRDS("de_vsnull")

setwd("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/3rd_year_project2/plots")
for (i in 1:K) {
  de_le_VlnPlot <- volcano_plot(de_vsnull,
                                k = K,
                                labels = colnames(counts)) +
    labs(title = paste0("DE with LFC=vsnull: k = ", i))

  ggsave(paste0("DE_with_LFC_VS_null_k_", i, ".pdf"), plot = de_le_VlnPlot)
}

de_le_VlnPlot
head(de_le$lfsr,2)
de_null_VlnPlot <- volcano_plot(de_vsnull,
                       k = K,
                       labels = genes$symbol
                       )
