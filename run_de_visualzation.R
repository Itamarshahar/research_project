################################################################################
## visualize DE
################################################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("please provide 3 args: <path_to_obj> <path_to_fit> <path_to_plots>")
  
} else {
  path_to_obj <- args[1] # Obj
  path_to_fit <- args[2] # Fit 
  path_to_de <- args[3] # DE 
  K <- args[4] # K value
  path_to_plots <- args[5]
}

################################################################################
## local run - REMOVE before running in the cluster
################################################################################
path_to_obj <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/databases/obj_and_topics.h5seurat"
path_to_fit <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/3rd_year_project2/fitted_k_10.rds"
path_to_de <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/databases/de_le"
path_to_plot <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/plots"
path_to_markers <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/cell_marker_24_8_23.csv"
K <- 10

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
library(viridisLite)
################################################################################
## Global variables
################################################################################
font <- "Arial"


additional_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                       "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                       "#cab2d6", "#6a3d9a")
################################################################################
## load objects 
################################################################################
print("Loading the objects - START")
obj <- LoadH5Seurat(path_to_obj)
fit <- readRDS(path_to_fit)
de <-readRDS(path_to_de)

print("Loading the objects - DONE")

################################################################################
## This is temp processes. While running the flow few rows (0 rows) of the 
#obj@meta.data we removed as part of the flow, this rows (0 rows) 
#where NOT remove from the Counts Matrix and therefore need to be removed 
################################################################################
counts <- obj@assays$RNA@counts
row_sum <- rowSums(counts)
nonzero_rows <- row_sum != 0
counts <- counts[nonzero_rows,]
counts <- t(counts)

################################################################################
## VolcanoPlot
################################################################################
setwd(path_to_plot)

de_VolcanoPlot <- list()
for (i in 1:K) {
  de_VolcanoPlot[[i]]<- volcano_plot(de,
                                     k = K,
                                     labels = colnames(counts)) + labs(title = paste0("DE with LFC vs-null Mode: k = ", i))
}

combined_plot <- grid.arrange(grobs = de_VolcanoPlot,
                              ncol= 2,
                              nrow= 5,
                              #widths = , 
                              #heights= , 
)
ggsave(combined_plot, file = paste0(path_to_plot, "\\" , "VolcanoPlot_DE_with_LFC_VS_null.pdf"))

################################################################################  
## VlnPlot grouping.by clusters 
################################################################################  

de_VlnPlot <- list()

for (i in 1:K) {
  de_VlnPlot[[i]] <- VlnPlot(obj,
                           features = paste0("k", i),
                           pt.size=0)+ NoLegend() + labs(x = NULL) 
}

      
combined_plot <- arrangeGrob(grobs = de_VlnPlot,
                                      ncol = 2,
                                      nrow = 5,
                             )
ggsave(combined_plot, 
       file = paste0(path_to_plot, "\\" ,"VlnPlots.pdf"),
       width =15,
       height=10)
 
################################################################################           
## VlnPlot grouping.by SampleID 
################################################################################           
de_VlnPlot <- list()

for (i in 1:K) {
  de_VlnPlot[[i]] <- VlnPlot(obj,
                             features = paste0("k", i),
                             pt.size=0,
                             group.by = "SampleID")+ NoLegend() + labs(x = NULL) 
}

?VlnPlot
combined_plot <- arrangeGrob(grobs = de_VlnPlot,
                             ncol = 2,
                             nrow = 5,
)
ggsave(combined_plot, 
       file = paste0(path_to_plot, "\\" ,"VlnPlotsSamplesID.pdf"),
       width =15,
       height=10)


################################################################################           
## VlnPlot grouping.by Gender 
################################################################################           
de_VlnPlot <- list()

for (i in 1:K) {
  de_VlnPlot[[i]] <- VlnPlot(obj,
                             features = paste0("k", i),
                             pt.size=0,
                             group.by = "Gender")+ NoLegend() + labs(x = NULL) 
}

?VlnPlot
combined_plot <- arrangeGrob(grobs = de_VlnPlot,
                             ncol = 2,
                             nrow = 5,
)
ggsave(combined_plot, 
       file = paste0(path_to_plot, "\\" ,"VlnPlotsGender.pdf"),
       width =15,
       height=10)

################################################################################           
## VlnPlot grouping.by Diagnosis 
################################################################################           
de_VlnPlot <- list()

for (i in 1:K) {
  de_VlnPlot[[i]] <- VlnPlot(obj,
                             features = paste0("k", i),
                             pt.size=0,
                             group.by = "Diagnosis")+ NoLegend() + labs(x = NULL) 
}


combined_plot <- arrangeGrob(grobs = de_VlnPlot,
                             ncol = 2,
                             nrow = 5,
)
ggsave(combined_plot, 
       file = paste0(path_to_plot, "\\" ,"VlnPlotsAccordingTo", "Diagnosis.pdf"),
       width =15,
       height=10)












################################################################################           
## 
################################################################################ 
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
