
################################################################################
## This script suppose to run the flow of the topics.
################################################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("please provide 3 args: <path_to_obj> <path_to_fit> <path_to_plots>", args)
  
} else {
path_to_obj <- args[1] # Seurat obj
path_to_fit <- args[2] # RDS obj
path_to_plots <- args[3]
}

load_libraries <- function() {
  print("loading libraries\n")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "http://cran.us.r-project.org")
  }

  library(remotes)
  if (!requireNamespace("fastTopics", quietly = TRUE)) {
    remotes::install_github("stephenslab/fastTopics")
  }
  library(RColorBrewer)
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
  library(grid)
  source("utils.R")
  print("loading libraries finished", "\n")

}

#dir.create(file.path(path_to_plots, paste0("Plots_for_k=", K)))
#path_to_plots <- paste0(path_to_plots, "/Plots_for_k=", K)

#
generate_structure_plot<- function (obj = obj,
                        path_to_plots = path_to_plots,
                        additional_colors = additional_colors,
                        K = K,
                        fit = fit
)
{
  structure_plot(fit,
                 topics = 1:K,
                 colors = additional_colors,
                 gap = 25,
                 grouping = obj$celltype,
                 )

  ggsave(paste0("structure_plot_celltype", ".pdf"),
         limitsize = FALSE,
         path=path_to_plots,
         width = 15,
         height = 10
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
        width = 15,
         height = 10
  )
  structure_plot(fit,
                 topics = 1:K,
                 colors = additional_colors,
                 gap = 25,
                 grouping = obj$Diagnosis,
                 )

  ggsave(paste0("structure_plot_diagnosis", ".pdf"),
         limitsize = FALSE,
         path=path_to_plots,
         width = 15,
         height = 10
  )
}
generate_umap <- function(obj,
                          path_to_plots,
                          custom_feature_names,
                          font = "Helvetica") {

  title_clusters <- "Cell expresison in different topics splitted by clusters"
  g <- arrangeGrob(grobs = FeaturePlot(object = obj,
                                         features = custom_feature_names,
                                         label = TRUE,
                                         combine = F,
                                         label.size = 8,
    ),
   top = textGrob(title_clusters,gp=gpar(fontsize=40,font=font))
  )
  ggsave(filename = "FeaturePlot_topics_distribution_over_clusters.pdf",
           plot = g,
           width = 35,
           height = 25,
           limitsize = FALSE,
           path = path_to_plots,
  )
  title_celltype <- "Cell expresison in different topics splitted by cell type"
  g <- arrangeGrob(grobs = FeaturePlot(object = obj,
                                         features = custom_feature_names,
                                         label = TRUE,
                                         combine = F,
                                         label.size = 8,
                                        split.by = "celltype",

    ),
   top = textGrob(title_celltype,gp=gpar(fontsize=40,font=font))
  )
ggsave(filename = "FeaturePlot_topics_distribution_over_celltype.pdf",
         plot = g,
         width = 35,
         height = 25,
         limitsize = FALSE,
         path = path_to_plots,
)

}

generate_DotPlot <- function (obj,custom_feature_names ) {
  g<- dotplot_topics(obj  = obj,
               topic_columns=custom_feature_names,
               group.by = "celltype",
               alpha.threshold=0.01,
               order=F,
              )  + labs(y = "Topics", x = "Cell Type", title = "Cells Type and the Average Expression of Over Different Topics ")

ggsave("DotPlot_topics_vs_celltype.pdf",
       plot = g,
       limitsize = FALSE,
       path=path_to_plots,
        width = 15,
         height = 10

)
}
#'
#' This function generates one or more violin plots for grouped data using the Seurat package.
#'
#' @param obj A Seurat object containing the data for plotting.
#' @param path_to_plot A character string specifying the path where the generated plots will be saved.
#' @param group_by_param A character string specifying the variable by which to group the data for plotting.
#' @param K An integer specifying the number of violin plots to generate.
#' @param LFC_type A character string specifying the log-fold change (LFC) type for plotting.
#' @param complex_case_text A character string (optional) for indiprinting a complex case in plot titles.
#'
#' @details The function generates violin plots for specified groups or clusters in the input Seurat object (`obj`).
#' Each group specified in `group_by_param` will be represented as a separate panel in the output plot.
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
plot_Vln <- function(obj,
                     path_to_plots,
                     group_by_param,
                     K=10,
                     complex_case_text=NA )
  {
  de_VlnPlot <- list()
  #K <- as.integer(K)
  for (i in 1:K) {
    de_VlnPlot[[i]] <- VlnPlot(obj,
                               features = paste0("k", i),
                               pt.size = 0,
                               group.by = group_by_param) +
      NoLegend() + labs(x = NULL) + theme(axis.text.x = element_text(size = 8))
  }
  combined_plot <- arrangeGrob(grobs = de_VlnPlot,
                               ncol = 1,
                               nrow = K,
  )
  if (!is.na(complex_case_text)) {
    group_by_param <- paste0(group_by_param, "_", complex_case_text)
  }
  ggsave(combined_plot,
         file = paste0(path_to_plots ,"/VlnPlot-",group_by_param, ".pdf"),
         width = 15,
         height = K*2.5)
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
generate_VlnPlot <- function(obj, path_to_plots, K=10) {
  print(K)
  #print("Generating Violin Plots...", "K=", as.integer(K))
  plot_Vln(obj,
           path_to_plots,
           "seurat_clusters",
           K)  # Violin plot grouping by Seurat clusters
  plot_Vln(obj,
           path_to_plots,
           "Diagnosis",
           K)       # Violin plot grouping by Diagnosis
  plot_Vln(obj,
           path_to_plots,
           "celltype",
           K)        # Violin plot grouping by cell type
  plot_Vln(obj,
           path_to_plots,
           "SampleID",
           K)        # Violin plot grouping by SampleID
  plot_Vln(obj,
           path_to_plots,
           "Age",
           K) # Violin plot grouping by Age


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
generate_vln_plot_complex <- function(obj,path_to_plots, K=10) {
  print("Generating Violin Plots Complex Case ...")

  plot_Vln(subset(x = obj, subset = Diagnosis == "AD"),
           path_to_plots,
           "celltype",
           K=K,
           complex_case_text = "AD") # Violin plot to only AD grouping by celltype

  plot_Vln(subset(x = obj, subset = Diagnosis == "MCI"),
           path_to_plots,
           "celltype",
           K=K,
           complex_case_text = "MCI")

  plot_Vln(subset(x = obj, subset = Diagnosis == "Young CTRL"),
           path_to_plots,
           "celltype",
           K=K,
           complex_case_text = "YoungCTRL")

  plot_Vln(subset(x = obj, subset = Diagnosis == "HA"),
           path_to_plots,
           "celltype",
           K=K,
           complex_case_text = "HA")

  plot_Vln(subset(x = obj, subset = Diagnosis == "SuperAgers"),
           path_to_plots,
           "celltype",
           K=K,
           complex_case_text = "SuperAgers")

  plot_Vln(subset(x = obj, subset = Diagnosis == c("SuperAgers", "HA", "Young CTRL")),
           path_to_plots,
           group_by_param = "celltype",
           K=K,
           complex_case_text = "SuperAgers&HA&YoungCTRL")

  plot_Vln(subset(x = obj, subset = Diagnosis == c("SuperAgers", "HA")),
           path_to_plots,
           "celltype",
           K=K,
           complex_case_text = "SuperAgers&HA")

  plot_Vln(subset(x = obj, subset = Diagnosis == c("AD", "MCI")),
           path_to_plots,
           "celltype",
           K=K,
           complex_case_text = "AD&MCI")
  print("Done with Violin Plots...")

}
get_qualitative_colors <- function(n) {
  qual_col_pals = brewer.pal.info[brewer.pal.info$printegory == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  return (col_vector[1:n])
}
load_objects <- function(path_to_fit, path_to_obj){
   # load objects
  print("Loading the objects - START")
  obj <- LoadH5Seurat(path_to_obj)
  fit  <-readRDS(path_to_fit)
  print("Loading the objects - FINISH")

  return (list(obj = obj, fit = fit))
}
run_main_flow <- function (obj, fit) {
  load_libraries()
  font <- "Helvetica"


  obj <- combine_topics_and_meta_data(obj = obj, fit = fit)
  K <- dim(fit$L)[2]
  all_K <- paste0("k", 1:K)

  additional_colors <- get_qualitative_colors(K)

  print(paste0("The number of Topics are: ", K), "\n")
  generate_structure_plot(obj = obj,
                        path_to_plots = path_to_plots,
                        additional_colors = additional_colors,
                        K = K,
                        fit = fit)
  print("Generating UMAP...\n")
  # generate_umap(obj = obj,
  #               path_to_plots= path_to_plots,
  #               custom_feature_names=all_K,
  # )
  print("Finish UMAP...\n")

  generate_VlnPlot(obj = obj,
                   path_to_plots = path_to_plots,
                   K=K)

  generate_vln_plot_complex(obj = obj,
                            path_to_plots = path_to_plots,
                            K=K)
  generate_DotPlot (obj = obj,
                    custom_feature_names=all_K )


  print("visualize the plans (topics) for each cell")


  print(paste0("The Run is Finished. The Plots are saved at: ",path_to_plots))

}


main <- function() {
  #load_libraries()## local run - REMOVE before running in the cluster
  path_to_fit <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/objects/fitted_topic_model_k_15.rds"
  path_to_obj <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/objects/five_prec.h5seurat"
  path_to_plots <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/k15/plots"
  load_libraries()
  loads <- load_objects(path_to_fit, path_to_obj)
  obj <- loads$obj
  fit <- loads$fit
  run_main_flow(obj, fit)
}

main()
