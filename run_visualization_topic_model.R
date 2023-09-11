
################################################################################
## This script suppose to run the flow of the topics.
################################################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("please provide 3 args: <path_to_obj> <path_to_fit> <path_to_plots>")
  
} else {
path_to_obj <- args[1] # Seurat obj
path_to_fit <- args[2] # RDS obj
path_to_plots <- args[3]
}

path_to_obj <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/objects/five_prec.h5seurat"
path_to_fit <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/objects/fitted_topic_model_k_15.rds"
path_to_plots <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/plots"
################################################################################
##  libreries
################################################################################

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "http://cran.us.r-project.org")
}

library(remotes)
if (!requireNamespace("fastTopics", quietly = TRUE)) {
  remotes::install_github("stephenslab/fastTopics")
}

install.packages('remotes')

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
library(gridExtra)
library(dplyr)
library(gridExtra)
source("utils.R")

################################################################################
## # Load the object
################################################################################
print("Loading the objects - START")
obj <- LoadH5Seurat(path_to_obj)
fit <- readRDS(path_to_fit)
print("Loading the objects - DONE")

################################################################################
## Globals stuff
################################################################################
font <- "Arial"
K <- dim(fit$L)[2]

additional_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                       "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                       "#cab2d6", "#6a3d9a")
print(paste0("The number of Topics are:", K))



################################################################################
## Adding the topics as column to the obj@meta.data
################################################################################

obj@meta.data <- cbind(obj@meta.data, fit$L)

setwd("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/databases")
SaveH5Seurat(obj,
             filename= paste0("obj_and_topics_", K),
             overwrite = TRUE)
################################################################################
## Create folder for all the plots
################################################################################

dir.create(file.path(path_to_plots, paste0("Plots_for_k=", K)))
path_to_plots <- paste0(path_to_plots, "/Plots_for_k=", K)

################################################################################
## visualize the plans (topics) for each cell
################################################################################
print("visualize the plans (topics) for each cell")

structure_plot(fit,
               topics = 1:K,
               colors = additional_colors,
               gap = 25,
               grouping = obj$celltype,
               )

ggsave(paste0("structure_plot_celltype", ".pdf"),
       limitsize = FALSE,
       path=path_to_plots,
)
structure_plot(fit,
               topics = 1:K,
               colors = additional_colors,
               gap = 25,
               grouping = obj$seurat_clusters,
)
ggsave(paste0("structure_plot_clusters", ".pdf"),
       limitsize = FALSE,
       path=path_to_plots,
)
################################################################################
## FeaturePlot - umap
################################################################################
custom_feature_names <- paste0("k", 1:K)
g<-arrangeGrob(grobs=FeaturePlot(object = obj,
            features =custom_feature_names,
            label = TRUE,
            combine = F,
            ))
ggsave("FeaturePlot_topics_distribution_over_clusters.pdf",
       plot = g,
       width = 35,
       height = 25,
       limitsize = FALSE,
       path=path_to_plots,
)
g<-arrangeGrob(grobs=FeaturePlot(object = obj,
                                  features =custom_feature_names,
                                  label = TRUE,
                                  combine = F,
                                  split.by = "celltype",
))
ggsave("FeaturePlot_topics_distribution_over_celltype.pdf",
       plot = g,
       width = 35,
       height = 25,
       limitsize = FALSE,
       path=path_to_plots,
)

################################################################################
## DotPlot
################################################################################

g<- dotplot_topics(obj  = obj,
               topic_columns=custom_feature_names,
               group.by = "celltype",
               alpha.threshold=0.01,
               order=F,
              )  + labs(y = "Topics", x = "Cell Type", title = "Cells Type and the Average Expression of Over Different Topics ")

ggsave("DotPlot_topics_vs_celltype.pdf",
       limitsize = FALSE,
       path=path_to_plots,
)








