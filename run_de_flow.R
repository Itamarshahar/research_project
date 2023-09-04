################################################################################
## visualize DE
################################################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("please provide 3 args: <path_to_obj> <path_to_fit> <path_to_plots>")
  
} else {
  path_to_de <- args[1] # DE obj
  K <- args[2] # K value
  path_to_plots <- args[3]
}

################################################################################
## Global variables
################################################################################
font <- "Arial"
K <- 10

additional_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                       "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                       "#cab2d6", "#6a3d9a")

################################################################################
## Global imports and libraries
################################################################################
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



## tmp 
path_to_obj <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/3rd_year_project2/five_prec.h5seurat"
obj <- LoadH5Seurat(path_to_obj) 

path_to_fit <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/3rd_year_project2/fitted_k_10.rds"
fit <- readRDS(path_to_fit)

################################################################################
## load objects 
################################################################################
print("Loading the objects - START")
obj <- LoadH5Seurat(path_to_obj)
fit <- readRDS(path_to_fit)
print("Loading the objects - DONE")


counts <- obj@assays$RNA@counts
row_sum <- rowSums(counts)
nonzero_rows <- row_sum != 0
counts <- counts[nonzero_rows,]
counts <- t(counts)
#setwd("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/3rd_year_project2/plots")

## Topic Proportion - Matrix L (using UMAP or t-SNE)
setwd("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/plots")

path_to_de <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/databases/de_le"
de <-readRDS(path_to_de)
de_VolcanoPlot <- list()

for (i in 1:K) {
  de_VolcanoPlot[[i]]<- volcano_plot(de,
                                k = K,
                                labels = colnames(counts)) +
    labs(title = paste0("DE with LFC vs-null Mode: k = ", i))
}

combined_plot <- grid.arrange(grobs = de_le_VolcanoPlot,
                              ncol= 2,
                              nrow= 5,
                              #widths = , 
                              #heights= , 
                                )
ggsave(combined_plot, file = "VolcanoPlot_DE_with_LFC_VS_null_k_all1111.pdf")

################################################################################  
## Volcano Plot
################################################################################  

de_VlnPlot <- list()

for (i in 1:K) {
  de_VlnPlot[[i]]<- VlnPlot(obj,
                            features = counts,
                            #labels = colnames(counts), 
                            ) +
    labs(title = paste0("DE with LFC vs-null Mode: k = ", i))
}
?VlnPlot
combined_plot <- grid.arrange(grobs = de_VlnPlot,
                              ncol= 2,
                              nrow= 5,
                              #widths = , 
                              #heights= , 
)
combined_plot
ggsave(combined_plot, file = "VlnPlots_DE_with_LFC_VS_null.pdf")

vl

