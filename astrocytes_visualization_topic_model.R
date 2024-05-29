################################################################################
## This script suppose to run the flow of the topics.
################################################################################


load_libraries <- function() {
  library(logger)
  log_info("loading libraries")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "http://cran.us.r-project.org")
  }

  library(remotes)
  if (!requireNamespace("fastTopics", quietly = TRUE)) {
    remotes::install_github("stephenslab/fastTopics")
  }
  source("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/utils/constants.R")
  source("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/utils/utils.R")
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
  library(logger)

  #library(SeuratData)
  library(magrittr)
  library(gridExtra)
  library(dplyr)
  library(grid)
  source("./utils/utils.R")
  library(ComplexHeatmap)
  library(circlize)
  log_info("loading libraries finished")

}


#' Generate Structure Plot for Topics
#'
#' This function generates structure plots to visualize the relationship between topics and different groupings
#' such as Seurat clusters, and diagnosis.
#'
#' @param obj A Seurat object containing the data for plotting.
#' @param path_to_plots A character string specifying the path where the generated plots will be saved.
#' @param additional_colors A vector of color values for topics to be used in the plot.
#' @param K An integer specifying the number of topics to include in the plot.
#' @param fit The fitted topic model.
#'
#' @details The function creates structure plots for each of the specified groupings (Seurat clusters, and diagnosis)
#' to visualize how topics are distributed across these groupings. The plots show the relationship between topics (1 to K)
#' and the specified grouping variable.
#'
#' @return None
#' @export
#'
#' @examples
#' # Generate structure plots for topics vs. cell types, clusters, and diagnosis
#' generate_structure_plot(seurat_obj, path_to_plots = "plots", additional_colors = c("blue", "red"), K = 10)
#'
#' @seealso
#' \code{\link{structure_plot}}
#'
#' @import Seurat
#'
generate_structure_plot <- function(obj = obj,
                                    path_to_plots = path_to_plots,
                                    additional_colors = additional_colors,
                                    K = K,
                                    fit = fit
)
{
  log_info("Generating structure_plot...")

  structure_plot(fit,
                 topics = 1:K,
                 colors = additional_colors,
                 gap = 25,
                 grouping = obj$seurat_clusters,
  )
  ggsave(paste0("structure_plot_clusters", ".pdf"),
         limitsize = FALSE,
         path = path_to_plots,
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
         path = path_to_plots,
         width = 15,
         height = 10
  )
  structure_plot(fit,
                 topics = 1:K,
                 colors = additional_colors,
                 gap = 25,
                 grouping = obj$SampleID,
  )

  ggsave(paste0("structure_plot_", SAMPLE_ID, ".pdf"),
         limitsize = FALSE,
         path = path_to_plots,
         width = 15,
         height = 10
  )
  log_info("Finish structure_plot... ")
}

#' Generate UMAP-based Feature Plots
#'
#' This function generates UMAP-based feature plots to visualize the expression of custom features
#' across different topics and split by clusters or cell types.
#'
#' @param obj A Seurat object containing the data for plotting.
#' @param path_to_plots A character string specifying the path where the generated plots will be saved.
#' @param custom_feature_names A character vector specifying the names of custom features to be visualized.
#' @param font A character string specifying the font to be used for plot titles.
#'
#' @details The function generates two types of feature plots based on UMAP coordinates:
#' - Feature plot for different topics split by clusters.
#' - Feature plot for different topics split by cell types.
#' Each feature plot displays the expression of custom features and labels for the specified features.
#'
#' @return None
#' @export
#'
#' @examples
#' # Generate UMAP-based feature plots for custom features
#' generate_umap(seurat_obj, path_to_plots = "plots", custom_feature_names = c("Gene1", "Gene2"), font = "Arial")
#'
#' @seealso
#' \code{\link{FeaturePlot}}, \code{\link{arrangeGrob}}
#'
#' @import Seurat
#'
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
                   top = textGrob(title_clusters, gp = gpar(fontsize = 40, font = font))
  )
  ggsave(filename = "FeaturePlot_topics_distribution_over_clusters.pdf",
         plot = g,
         width = 35,
         height = 25,
         limitsize = FALSE,
         path = path_to_plots,
  )

  title_SampleID <- "Cell expresison in different topics splitted by SampleID"
  g <- arrangeGrob(grobs = FeaturePlot(object = obj,
                                       features = custom_feature_names,
                                       label = TRUE,
                                       combine = F,
                                       label.size = 8,
                                       #split.by = SAMPLE_ID,
                                       split.by = DIAGNOSIS_SAMPLE,

  ),
                   # top = textGrob(title_celltype, gp = gpar(fontsize = 40, font = font))
  )
  ggsave(filename = "FeaturePlot_topics_distribution_over_SampleID.pdf",
         plot = g,
         width = 35,
         height = 25,
         limitsize = FALSE,
         path = path_to_plots,
  )

  title_Diagnosis <- "Cell expresison in different topics splitted by Diagnosis"
  g <- arrangeGrob(grobs = FeaturePlot(object = obj,
                                       features = custom_feature_names,
                                       label = TRUE,
                                       combine = F,
                                       label.size = 8,
                                       split.by = "Diagnosis",

  ),
                   # top = textGrob(title_celltype, gp = gpar(fontsize = 40, font = font))
  )
  ggsave(filename = "FeaturePlot_topics_distribution_over_Diagnosis.pdf",
         plot = g,
         width = 35,
         height = 25,
         limitsize = FALSE,
         path = path_to_plots,
  )

}

#' Generate Dot Plot for Topics vs. Cell Type
#'
#' This function generates a dot plot to visualize the relationship between topics and cell types.
#'
#' @param obj A Seurat object containing the data for plotting.
#' @param custom_feature_names A character vector specifying the names of topics or features to plot.
#'
#' @details The function creates a dot plot to display the average expression of different topics
#' across various cell types in the Seurat object. Each dot represents the average expression of a topic
#' within a specific cell type.
#'
#' @return None
#' @export
#'
#' @examples
#' # Generate a dot plot for topics vs. cell types
#' generate_DotPlot(seurat_obj, custom_feature_names = c("k1", "k2", "k3"))
#'
#' @seealso
#' \code{\link{dotplot_topics}}
#'
#' @import Seurat
#'
generate_DotPlot <- function(obj, custom_feature_names, path_to_plots) {
  log_info("Generating DotPlots...")

  g_Age <- dotplot_topics(obj = obj,
                          topic_columns = custom_feature_names,
                          group.by = AGE,
                          alpha.threshold = 0.01,
                          order = F,
  ) + labs(y = "Topics", x = AGE, title = "Age and the Average Expression of Over Different Topics ") & RotatedAxis()

  g_Diagnosis <- dotplot_topics(obj = obj,
                                topic_columns = custom_feature_names,
                                group.by = DIAGNOSIS,
                                alpha.threshold = 0.01,
                                order = F,
  ) + labs(y = "Topics", x = DIAGNOSIS, title = "Diagnosis and the Average Expression of Over Different Topics ") & RotatedAxis()


  g_SampleID <- dotplot_topics(obj = obj,
                               topic_columns = custom_feature_names,
                               # group.by = "Diagnosis_SampleID",
                               group.by = DIAGNOSIS_SAMPLE,
                               #group.by = SAMPLE_ID,
                               alpha.threshold = 0.01,
                               order = F,
                               #mean_threshold = 0.001
  ) &
    labs(y = "Topics", x = SAMPLE_ID, title = "SampleID and the Average Expression of Over Different Topics") &
    RotatedAxis()

  g_SampleID_2 <- dotplot_topics(obj = obj,
                                 topic_columns = custom_feature_names,
                                 group.by = DIAGNOSIS_SAMPLE,
                                 # group.by = "Diagnosis_SampleID",
                                 #group.by = SAMPLE_ID,
                                 alpha.threshold = 0.001,
                                 order = F,
                                 mean_threshold = 0.001
  ) &
    labs(y = "Topics", x = SAMPLE_ID, title = "SampleID and the Average Expression of Over Different Topics") &
    RotatedAxis()

  g_seurat_clusters <- dotplot_topics(obj = obj,
                                      topic_columns = custom_feature_names,
                                      group.by = SEURAT_CLUSTERS,
                                      alpha.threshold = 0.01,
                                      order = F,
  ) + labs(y = "Topics", x = SEURAT_CLUSTERS, title = "seurat_clusters and the Average Expression of Over Different Topics ") & RotatedAxis()
  plot_grid(g_Age, g_Diagnosis, g_SampleID, g_SampleID_2, g_seurat_clusters, nrow = 5, ncol = 1)
  ggsave("DotPlot_topics_vs_different_features.pdf",
         #plot = g_all,
         limitsize = FALSE,
         path = path_to_plots,
         width = 15,
         height = 30

  )

}

#' Generate Violin Plots
#'
#' This function generates a set of violin plots for grouped data using the Seurat package.
#'
#' @param obj A Seurat object containing the data for plotting.
#' @param path_to_plots A character string specifying the path where the generated plots will be saved.
#' @param group_by_param A character string specifying the variable by which to group the data for plotting.
#' @param K An integer specifying the number of violin plots to generate.
#' @param complex_case_text A character string (optional) for incorporating a complex case in plot titles.
#'
#' @details The function generates violin plots for specified groups or clusters in the input Seurat object (`obj`).
#' Each group specified in `group_by_param` will be represented as a separate panel in the output plot.
#'
#' @return None
#' @export
#'
#' @examples
#' # Generate violin plots for Seurat clusters
#' plot_Vln(seurat_obj, "path/to/plots", SEURAT_CLUSTERS, K = 10, complex_case_text = "ComplexCase")
#'
#' # Generate violin plots for other variables
#' plot_Vln(seurat_obj, "path/to/plots", DIAGNOSIS, K = 10)
#'
#' @seealso
#' \code{\link{generate_VlnPlot}}, \code{\link{generate_vln_plot_complex}}
#'
#' @import Seurat
#'
plot_Vln <- function(obj,
                     path_to_plots,
                     group_by_param,
                     K = 10,
                     complex_case_text = NA,
                     do_group_by = TRUE)
{
  de_VlnPlot <- list()
  for (i in 1:K) {
    de_VlnPlot[[i]] <- VlnPlot(obj,
                               features = generate_names_for_fit_row_and_cols(K, generate_single_name = T, i),
                               # features = paste0("X.10.", "k", i),
                               pt.size = 0,
                               group.by = group_by_param) +
      NoLegend() +
      labs(x = NULL) +
      theme(axis.text.x = element_text(size = 8))
  }
  combined_plot <- arrangeGrob(grobs = de_VlnPlot,
                               ncol = 1,
                               nrow = K,
  )
  if (!is.na(complex_case_text)) {
    group_by_param <- paste0(group_by_param, "_", complex_case_text)
  }
  ggsave(combined_plot,
         # file = paste0(path_to_plots, "/VlnPlot-", group_by_param, ".pdf"),
         file = paste0(path_to_plots, "/", "VlnPlot-", group_by_param, ".pdf"),
         width = 15,
         height = K * 2.5)
}

generate_VlnPlot <- function(obj, path_to_plots, K = 10) {
  plot_Vln(obj,
           path_to_plots,
           group_by_param = SEURAT_CLUSTERS,
           K)  # Violin plot grouping by Seurat clusters
  plot_Vln(obj,
           path_to_plots,
           DIAGNOSIS,
           K)       # Violin plot grouping by Diagnosis
  plot_Vln(obj,
           path_to_plots,
           DIAGNOSIS_SAMPLE,
           K)        # Violin plot grouping by SampleID
  plot_Vln(obj,
           path_to_plots,
           AGE,
           K) # Violin plot grouping by Age
  log_info("Done with Violin Plots...")

}

#' Generate Multiple Violin Plots
#'
#' This function generates a set of violin plots for different grouping variables using the Seurat package.
#'
#' @param obj A Seurat object containing the data for plotting.
#' @param path_to_plots A character string specifying the path where the generated plots will be saved.
#' @param K An integer specifying the number of violin plots to generate.
#'
#' @details The function calls the `plot_Vln` function for different grouping variables to generate multiple violin plots.
#'
#' @return None
#' @export
#'
#' @examples
#' # Generate multiple violin plots
#' generate_VlnPlot(seurat_obj, "path/to/plots", K = 10)
#'
#' @seealso
#' \code{\link{plot_Vln}}, \code{\link{generate_vln_plot_complex}}
#'
#' @import Seurat
#'
generate_vln_plot_complex <- function(obj, path_to_plots, K = 10) {

  log_info("Generating Violin Plots Complex Case ...")
  ggsave(filename = glue::glue(path_to_plots, "/", "VlnPlot-sick_vs_healthy.pdf"),
         plot = VlnPlot(object = obj,
                        features = generate_names_for_fit_row_and_cols(K),
                        pt.size = 0,
                        group.by = "sick_vs_healthy"),
         width = 20,
         height = K * 1.6)

  if (old_plot) {
  plot_Vln(subset(x = obj, subset = Diagnosis == "AD"),
           path_to_plots,
           group_by_param = "cell.type",
           K = K,
           complex_case_text = "AD") # Violin plot to only AD grouping by celltype

  plot_Vln(subset(x = obj, subset = Diagnosis == "MCI"),
           path_to_plots,
           group_by_param = "cell.type",
           K = K,
           complex_case_text = "MCI")

  plot_Vln(subset(x = obj, subset = Diagnosis == "Young CTRL"),
           path_to_plots,
           group_by_param = "cell.type",
           K = K,
           complex_case_text = "YoungCTRL")

  plot_Vln(subset(x = obj, subset = Diagnosis == "HA"),
           path_to_plots,
           group_by_param = "cell.type",
           K = K,
           complex_case_text = "HA")

  plot_Vln(subset(x = obj, subset = Diagnosis == "SuperAgers"),
           path_to_plots,
           group_by_param = "cell.type",
           K = K,
           complex_case_text = "SuperAgers")

  plot_Vln(subset(x = obj, subset = Diagnosis == c("SuperAgers", "HA", "Young CTRL")),
           path_to_plots,
           group_by_param = "cell.type",
           K = K,
           complex_case_text = "SuperAgers&HA&YoungCTRL")

  plot_Vln(subset(x = obj, subset = Diagnosis == c("SuperAgers", "HA")),
           path_to_plots,
           group_by_param = "cell.type",
           K = K,
           complex_case_text = "SuperAgers&HA")

  plot_Vln(subset(x = obj, subset = Diagnosis == c("AD", "MCI")),
           path_to_plots,
           group_by_param = "cell.type",
           K = K,
           complex_case_text = "AD&MCI")
  }
  log_info("Done with Violin Plots...")
  #
}

get_qualitative_colors <- function(n) {
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  return(col_vector[1:n])
}

load_objects <- function(path_to_fit, path_to_obj, rename_fit_cols_rows = TRUE, generate_columns_for_complex_vln = TRUE) {
  # load objects
  log_info("Loading the objects - START")
  obj <- readRDS(path_to_obj)
  fit <- readRDS(path_to_fit)
  if (rename_fit_cols_rows) {
    all_K <- generate_names_for_fit_row_and_cols(K = as.integer(dim(fit$L)[2]))
    fit <- rename_fit_row_and_cols(fit, all_K)
  }
  log_info("Loading the objects - FINISH")
  if (generate_columns_for_complex_vln)
  {
    obj@meta.data$sick_vs_healthy <- ifelse(obj@meta.data$Diagnosis %in% c("AD", "MCI"), yes = "AD_MCI",
                                            no = ifelse(obj@meta.data$Diagnosis %in% c("SuperAger", "HA"), yes = "SuperAger_HA", no = "YOUNG_CTRL"))
  }
  return(list(obj = obj, fit = fit))
}

run_main_flow <- function(obj, fit, path_to_plots) {
  load_libraries()
  obj <- combine_topics_and_meta_data(obj = obj, fit = fit)
  K <- as.integer(dim(fit$L)[2])
  all_K <- generate_names_for_fit_row_and_cols(K = K)
  additional_colors <- get_qualitative_colors(K)
  generate_DotPlot(obj = obj,
                   custom_feature_names = all_K,
                   path_to_plots = path_to_plots)
  generate_structure_plot(obj = obj,
                          path_to_plots = path_to_plots,
                          additional_colors = additional_colors,
                          K = K,
                          fit = fit)


  log_info("Generating UMAP...")
  generate_umap(obj = obj,
                path_to_plots = path_to_plots,
                custom_feature_names = all_K,
  )
  log_info("Finish UMAP...")
  log_info("Generating VlnPlot...")

  generate_VlnPlot(obj = obj,
                   path_to_plots = path_to_plots,
                   K = K)

  generate_vln_plot_complex(obj = obj,
                            path_to_plots = path_to_plots,
                            K = K)
  log_info("Finish complex VlnPlot")
  log_info(paste0("The Run is Finished. The Plots are saved at: ", path_to_plots))

}

main <- function() {
  source("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/astrocytes_visualization_topic_model.R")
  load_libraries()
  path_to_obj <- "/Volumes/habib-lab/shmuel.cohen/astrocytes/objects/filltered_astrocytes_v1.rds"

  for (K in 14:17) {
    log_info(paste0("Running for K = ", K))
    path_to_fit <- glue::glue("/Volumes/habib-lab/shmuel.cohen/astrocytes/objects/fitted_topic_model_k_", K, ".rds")
    path_to_plots <- glue::glue("/Volumes/habib-lab/shmuel.cohen/astrocytes/plots/k", K)
    loads <- load_objects(path_to_fit, path_to_obj)
    run_main_flow(loads$obj, loads$fit, path_to_plots)
    log_info(paste0("Finish for K = ", K))
  }
}
