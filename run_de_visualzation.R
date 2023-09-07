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


# VlnPlot -----------------------------------------------------------------

##' Generate Violin Plot(s) for Grouped Data
#'
#' This function generates one or more violin plots for grouped data using the Seurat package.
#'
#' @param obj A Seurat object containing the data for plotting.
#' @param path_to_plot A character string specifying the path where the generated plots will be saved.
#' @param group_by_parm A character string specifying the variable by which to group the data for plotting.
#' @param K An integer specifying the number of violin plots to generate.
#' @param LFC_type A character string specifying the log-fold change (LFC) type for plotting.
#' @param complex_case_text A character string (optional) for indicating a complex case in plot titles.
#'
#' @details The function generates violin plots for specified groups or clusters in the input Seurat object (`obj`).
#' Each group specified in `group_by_parm` will be represented as a separate panel in the output plot.
#' 
#' @return None
#' @export
#' 
#' @examples
#' # Generate violin plots for Seurat clusters
#' plot_vln(seurat_obj, "path/to/plots", "seurat_clusters", K = 10, LFC_type = "le")
#' 
#' # Generate violin plots for other variables
#' plot_vln(seurat_obj, "path/to/plots", "Diagnosis", K = 10, LFC_type = "le")
#' 
#' @seealso
#' \code{\link{generate_vln_plot}}, \code{\link{generate_vln_plot_complex}}
#' 
#' @import Seurat
#' 
plot_vln <- function(obj,
                     path_to_plot,
                     group_by_parm,
                     K=10 ,
                     LFC_type="le",
                     complex_case_text=NA ) 
  {
  
  de_VlnPlot <- list()
  K <- as.integer(K)
  for (i in 1:K) {
    de_VlnPlot[[i]] <- VlnPlot(obj,
                               features = paste0("k", i),
                               pt.size = 0,
                               group.by = group_by_parm) +
      NoLegend() + labs(x = NULL) + theme(axis.text.x = element_text(size = 8))
  }
  combined_plot <- arrangeGrob(grobs = de_VlnPlot,
                               ncol = 2,
                               nrow = 5,
  )
  if (!is.na(complex_case_text)) {
    group_by_parm <- paste0(group_by_parm, "_", complex_case_text)
  }
  ggsave(combined_plot, 
         file = paste0(path_to_plot ,"VlnPlot-", "LFC_TYPE_",LFC_type,group_by_parm, ".pdf"),
         width = 15,
         height = 10)
}

#' Generate Multiple Violin Plots
#'
#' This function generates a set of violin plots for different grouping variables.
#'
#' @param obj A Seurat object containing the data for plotting.
#' @param path_to_plot A character string specifying the path where the generated plots will be saved.
#' @param K An integer specifying the number of violin plots to generate.
#' @param LFC_type A character string specifying the log-fold change (LFC) type for plotting.
#'
#' @details The function calls `plot_vln` with different grouping variables to generate multiple violin plots.
#' 
#' @return None
#' @export
#' 
#' @examples
#' # Generate multiple violin plots
#' generate_vln_plot(seurat_obj, "path/to/plots", K = 10, LFC_type = "le")
#' 
#' @seealso
#' \code{\link{plot_vln}}, \code{\link{generate_vln_plot_complex}}
#' 
#' @import Seurat
#' 
generate_vln_plot <- function(obj,path_to_plot, K=10, LFC_type="le") {
  print("Generating Violin Plots...")
  plot_vln(obj, 
           path_to_plot, 
           "seurat_clusters",  
           K=10, 
           LFC_type="le")  # Violin plot grouping by Seurat clusters
  plot_vln(obj,
           path_to_plot,  
           "Diagnosis",  
           K=10, 
           LFC_type="le")       # Violin plot grouping by Diagnosis
  plot_vln(obj, 
           path_to_plot, 
           "celltype",  
           K=10, 
           LFC_type="le")        # Violin plot grouping by cell type
  plot_vln(obj,
           path_to_plot,
           "SampleID",
           K=10,
           LFC_type="le")        # Violin plot grouping by SampleID
  plot_vln(obj, 
           path_to_plot,
           "Age",  
           K=10,
           LFC_type="le") # Violin plot grouping by Age

  
  print("Done with Violin Plots...")
  
}
#' Generate Violin Plots for Complex Cases
#'
#' This function generates a set of violin plots for complex cases in differential expression analysis.
#'
#' @param obj A Seurat object containing the data for plotting.
#' @param path_to_plot A character string specifying the path where the generated plots will be saved.
#' @param K An integer specifying the number of violin plots to generate.
#' @param LFC_type A character string specifying the log-fold change (LFC) type for plotting.
#'
#' @details The function generates violin plots for complex cases by grouping the data based on the "Diagnosis" variable.
#' It creates separate plots for each complex case specified in the "Diagnosis" variable and saves them with appropriate titles.
#' 
#' @return None
#' @export
#' 
#' @examples
#' # Generate violin plots for complex cases
#' generate_vln_plot_complex(seurat_obj, "path/to/plots", K = 10, LFC_type = "le")
#' 
#' @seealso
#' \code{\link{plot_vln}}, \code{\link{generate_vln_plot}}
#' 
#' @import Seurat
#' 

generate_vln_plot_complex <- function(obj,path_to_plot, K=10, LFC_type="le") {
  print("Generating Violin Plots Complex Case ...")
  
  plot_vln(subset(x = obj, subset = Diagnosis == "AD"),
           path_to_plot, 
           "celltype",  
           K=10,
           LFC_type="le", 
           complex_case_text = "AD") # Violin plot to only AD grouping by celltype
  
  plot_vln(subset(x = obj, subset = Diagnosis == "MCI"),
           path_to_plot,
           "celltype",
           K=10,
           LFC_type="le", 
           complex_case_text = "MCI")
  
  plot_vln(subset(x = obj, subset = Diagnosis == "Young CTRL"), 
           path_to_plot, 
           "celltype",  
           K=10,
           LFC_type="le", 
           complex_case_text = "YoungCTRL") 

  plot_vln(subset(x = obj, subset = Diagnosis == "HA"),
           path_to_plot,
           "celltype", 
           K=10,
           LFC_type="le", 
           complex_case_text = "HA") 
  
  plot_vln(subset(x = obj, subset = Diagnosis == "SuperAgers"),
           path_to_plot, 
           "celltype", 
           K=10,
           LFC_type="le", 
           complex_case_text = "SuperAgers") 
  
  plot_vln(subset(x = obj, subset = Diagnosis == c("SuperAgers", "HA", "Young CTRL")),
           path_to_plot, 
           "celltype", 
           K=10,
           LFC_type="le", 
           complex_case_text = "SuperAgers&HA&YoungCTRL")
  
  plot_vln(subset(x = obj, subset = Diagnosis == c("SuperAgers", "HA")),
           path_to_plot, 
           "celltype", 
           K=10,
           LFC_type="le", 
           complex_case_text = "SuperAgers&HA")
  
  plot_vln(subset(x = obj, subset = Diagnosis == c("AD", "MCI")),
           path_to_plot, 
           "celltype", 
           K=10,
           LFC_type="le", 
           complex_case_text = "AD&MCI")
  print("Done with Violin Plots...")
  
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
  
  generate_vln_plot(obj = obj,
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



