###############################################################################
###############################################################################
# Run the DE visualization flow:


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
  library(grid)
  library(viridisLite)
  library(ggpointdensity)
  library(reshape2)
  library(glue)
  source("utils.R")
  library(ComplexHeatmap)
  library(circlize)
  library(combinat)
  library(glue)
  library(stringr)
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
  return(counts)
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
generate_VolcanoPlot <- function(de, path_to_plots, counts, K = 10, LFC_type = "VS_NULL") {
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
         file = paste0(path_to_plots,
                       "VolcanoPlot_DE_LFC_",
                       LFC_type,
                       ".pdf"),
         width = 12,
         height = 20)
  print("Done with Volcano Plots")
}

generate_ScatterPlot <- function(fit, de, variable_name = "topic", LFC_type = "VS_NULL", path_to_plots, K = 10
) {

  ###############################################################################
  # scheme of de_enriched_matrix
  #       k1  k2  k3  k4  k5  k6  k7  k8  k9 k10  f_value lfsr  lfc
  # geneA
  # geneB
  # geneC
  #
  #       k
  # geneA <the_val_of_geneA_in_k=1>
  # geneA <the_val_of_geneA_in_k=2>
  # geneA <the_val_of_geneA_in_k=3>
  ###############################################################################

  merge_by <- c("gene_name", variable_name)

  mat_f <- as.data.frame(fit$F)
  #combine mat_f and de$postmean
  mat_f$gene_name <- as.factor(rownames(mat_f))
  melted_mat_f <- melt(mat_f,
                       id.vars = "gene_name",
                       variable.name = variable_name,
                       value.name = "k_values",
  )
  melted_mat_f <- melted_mat_f[order(melted_mat_f$gene_name),]

  gene_score <- as.data.frame(de$F)
  gene_score$gene_name <- as.factor(rownames(gene_score))
  melted_gene_score <- melt(mat_f,
                            id.vars = "gene_name",
                            variable.name = variable_name,
                            value.name = "gene_score",
  )
  melted_gene_score <- melted_gene_score[order(melted_gene_score$gene_name),]

  lfc <- as.data.frame(de$postmean)
  lfc$gene_name <- as.factor(rownames(lfc))

  melted_lfc <- melt(mat_f,
                     id.vars = "gene_name",
                     variable.name = variable_name,
                     value.name = "lfc",
  )
  melted_lfc <- melted_lfc[order(melted_lfc$gene_name),]
  lfsr <- as.data.frame(de$lfsr)
  lfsr$gene_name <- as.factor(rownames(lfsr))
  melted_lfsr <- melt(mat_f,
                      id.vars = "gene_name",
                      variable.name = variable_name,
                      value.name = "lfsr",
  )

  melted_lfsr <- melted_lfsr[order(melted_lfsr$gene_name),]


  combined <- merge(melted_mat_f, melted_gene_score, by = merge_by)
  combined <- merge(combined, melted_lfsr, by = merge_by)
  combined <- merge(combined, melted_lfc, by = merge_by)


  de_enriched <- combined[order(combined$gene_name, combined$k),]
  rownames(combined) <- seq(1, nrow(combined))

  generate_ScatterPlot_DE(de_enriched = de_enriched,
                          path_to_plots = path_to_plots,
                          K = K,
                          file_name = glue(LFC_type, "_ScatterPlot.pdf")
  )
}


load_objects_fit_obj_de <- function(path_to_fit, path_to_de, path_to_obj = NA) {
  # load objects
  print("Loading the objects - START")

  if (!is.na(path_to_obj)) {
    obj <- LoadH5Seurat(path_to_obj)
    de <- readRDS(path_to_de)
    fit <- readRDS(path_to_fit)
    print("Loading the objects - FINISH")

    return(list(fit = fit, de = de, obj = obj))
  }
  de <- readRDS(path_to_de)
  fit <- readRDS(path_to_fit)
  print("Loading the objects - FINISH")

  return(list(fit = fit, de = de,))
}

# Main Flow Run -----------------------------------------------------------
run_de_visualzation <- function(K = 10,
                                LFC_type = "vs_null",
                                remove_zeros_from_counts_matrix = T) {
  load_libraries()
  print(glue("Runnig the DE visualization flow with K={K} and LFC_type={LFC_type}"))
  ## local run - REMOVE before running in the cluster
  plots <- "plots"
  k_folder <- glue("k", K)
  objects <- "objects"
  obj_name <- "obj_and_topics"
  seurat_postfix <- ".h5seurat"
  fitted_obj_name <- glue("fitted_k_{K}.rds")
  left_path <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/SuperAgers_n=20/5_percent"
  de_obj_name <- "de_vsnull"
  obj_name <- "obj_and_topics"

  path_to_plots <- glue("{left_path}/{k_folder}/{plots}/")
  path_to_de <- glue("{left_path}/{k_folder}/{objects}/{de_obj_name}")
  path_to_fit <- glue("{left_path}/{k_folder}/{objects}/{fitted_obj_name}")
  path_to_obj <- glue("{left_path}/{k_folder}/{objects}/{obj_name}{seurat_postfix}")

  loaded <- load_objects_fit_obj_de(path_to_fit, path_to_de, path_to_obj)
  if (remove_zeros_from_counts_matrix == TRUE) {
    counts <- remove_zeros_from_counts_mat(loaded$obj)
  } else {
    counts <- loaded$obj@assays$RNA@counts
  }

  generate_VolcanoPlot(de = loaded$de,
                       path_to_plots = path_to_plots,
                       counts = counts,
                       K = as.integer(K),
                       LFC_type = LFC_type)
  generate_ScatterPlot(fit = loaded$fit,
                       de = loaded$de,
                       path_to_plots = path_to_plots,
                       K = as.integer(K),
                       LFC_type = LFC_type)

  #print(glue("The Run is Finished. The Plots are saved at: ", path_to_plots))

}


###############################################################################
###############################################################################

# Run the topic evaluations flow:


generate_all_permutations <- function(lst) {
  # Initialize an empty vector to store the permutations
  all_permutations <- character(0)

  # Generate unique permutations of length 1
  for (i in 1:length(lst)) {
    all_permutations <- c(all_permutations, paste(lst[i], lst[i], sep = "_"))
  }

  # Generate unique permutations of length 2
  for (i in 1:length(lst)) {
    for (j in i:length(lst)) {
      if (i != j) {
        all_permutations <- c(all_permutations, paste(lst[i], lst[j], sep = "_"))
      }
    }
  }

  return(all_permutations)
}

extract_k <- function(path_to_fit) {
  fit <- readRDS(path_to_fit)
  return(list(fit = fit, k = dim(fit$F)[2]))
}

generate_fits_list <- function(paths) {
  fits_list <- list()  # Initialize the fits_list

  for (path in paths) {
    obj <- extract_k(path)
    fit <- obj$fit
    k <- obj$k
    fits_list[[as.character(k)]] <- fit  # Use as.character to ensure k is a character key
  }
  return(fits_list)
}

helper_topic_evaluation <- function(fits_list, path_to_plots, type = "cells", correlation_method = "pearson") {
  print(glue("Running topic evaluation flow for {type} and {correlation_method}"))
  col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  col_fun(seq(-20, 20))

  fits_list <- generate_fits_list(fits_list)
  all_permutations <- generate_all_permutations(names(fits_list))
  #print(all_permutations)

  for (per in all_permutations) {
    #print(per)
    split_parts <- unlist(strsplit(per, "_"))
    k_left <- split_parts[1]
    k_right <- split_parts[2]
    fit_k_left <- fits_list[[k_left]]
    fit_k_right <- fits_list[[k_right]]
    if (type == "cells") {
      correlation <- claculate_correlation(fit_k_left$L, fit_k_right$L, method = correlation_method)
      file_name <- glue('{correlation_method}_corrlation_between_{type}_')
    }
    else {
      correlation <- claculate_correlation(fit_k_left$F, fit_k_right$F, method = correlation_method)
      file_name <- glue('{correlation_method}_corrlation_between_genes_')
    }

    pdf(glue("{path_to_plots}{file_name}k={k_left}_with_k={k_right}.pdf"))
    draw(Heatmap(correlation,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 col = col_fun,
                 column_title = glue("The Corralation Between K={k_left} with K={k_right} Topics"),
                 column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                 name = "Correlation",
                 rect_gp = gpar(col = "white", lwd = 2),
                 column_names_rot = 45,
                 cell_fun = HeatmapHelper_add_values_to_display(correlation = correlation,
                                                                OnlyPositive = TRUE)
    )
    )
    dev.off()
    load_libraries()
  }
  #print("Done with the heatmaps")

}

run_topic_evaluation <- function(fit_list = NA, local_run = TRUE) {
  print("Loading libraries")
  load_libraries()
  if (local_run) {
    obj_folder <- "objects"
    fit_name_k_10 <- glue("fitted_k_", 10, ".rds")
    fit_name_k_15 <- glue("fitted_k_", 15, ".rds")
    fit_name_k_12 <- glue("fitted_k_", 12, ".rds")
    fit_name_k_20 <- glue("fitted_k_", 20, ".rds")
    initial_object_name <- "5_precent.h5seurat"
    plots <- "plots"
    left_path <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/SuperAgers_n=20/5_percent"

    path_to_plots <- glue("{left_path}/{plots}/")
    #path_to_fit_k_10 <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/SuperAgers_n=20/5_percent/k10/objects/fitted_k_10.rds"
    path_to_fit_k_10 <- glue("{left_path}/k10/{obj_folder}/{fit_name_k_10}")
    path_to_fit_k_15 <- glue("{left_path}/k15/{obj_folder}/{fit_name_k_15}")
    path_to_fit_k_12 <- glue("{left_path}/k12/{obj_folder}/{fit_name_k_12}")
    path_to_fit_k_20 <- glue("{left_path}/k20/{obj_folder}/{fit_name_k_20}")

  }

  #path_to_fit_k_15 <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/SuperAgers_n=20/5_percent/k15/objects/fitted_k_15.rds"
  helper_topic_evaluation(c(path_to_fit_k_10, path_to_fit_k_15, path_to_fit_k_12, path_to_fit_k_20), path_to_plots, type = "genes") #, correlation_method = "kendall")
  helper_topic_evaluation(c(path_to_fit_k_10, path_to_fit_k_15, path_to_fit_k_12, path_to_fit_k_20), path_to_plots, type = "cells") #, correlation_method = "kendall")
}

###############################################################################
###############################################################################
# Run the topic visualization flow:
###############################################################################
###############################################################################
################################################################################
## This script suppose to run the flow of the topics.
################################################################################
# args = commandArgs(trailingOnly = TRUE)
#
# if (length(args) != 3) {
#   stop("please provide 3 args: <path_to_obj> <path_to_fit> <path_to_plots>", args)
#
# } else {
#   path_to_obj <- args[1] # Seurat obj
#   path_to_fit <- args[2] # RDS obj
#   path_to_plots <- args[3]
# }


#' Generate Structure Plot for Topics
#'
#' This function generates structure plots to visualize the relationship between topics and different groupings
#' such as cell types, Seurat clusters, and diagnosis.
#'
#' @param obj A Seurat object containing the data for plotting.
#' @param path_to_plots A character string specifying the path where the generated plots will be saved.
#' @param additional_colors A vector of color values for topics to be used in the plot.
#' @param K An integer specifying the number of topics to include in the plot.
#' @param fit The fitted topic model.
#'
#' @details The function creates structure plots for each of the specified groupings (cell types, Seurat clusters, and diagnosis)
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
  structure_plot(fit,
                 topics = 1:K,
                 colors = additional_colors,
                 gap = 25,
                 grouping = obj$celltype,
  )

  ggsave(paste0("structure_plot_celltype", ".pdf"),
         limitsize = FALSE,
         path = path_to_plots,
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

  title_clusters <- str_to_title("Cell expresison in different topics splitted by clusters", locale = "en")
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
  title_celltype <- str_to_title("Cell expresison in different topics splitted by cell type", locale = "en")
  g <- arrangeGrob(grobs = FeaturePlot(object = obj,
                                       features = custom_feature_names,
                                       label = TRUE,
                                       combine = F,
                                       label.size = 8,
                                       split.by = "celltype",

  ),
                   top = textGrob(title_celltype, gp = gpar(fontsize = 40, font = font))
  )
  ggsave(filename = "FeaturePlot_topics_distribution_over_celltype.pdf",
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
  g_cell_type <- dotplot_topics(obj = obj,
                                topic_columns = custom_feature_names,
                                group.by = "celltype",
                                alpha.threshold = 0.01,
                                order = F,
  ) &
    RotatedAxis() &
    labs(y = "Topics", x = "Cell Type", title = "Cells Type and the Average Expression of Over Different Topics")

  g_Age <- dotplot_topics(obj = obj,
                          topic_columns = custom_feature_names,
                          group.by = "Age",
                          alpha.threshold = 0.01,
                          order = F,
  ) + labs(y = "Topics", x = "Age", title = "Age and the Average Expression of Over Different Topics ") & RotatedAxis()

  g_Diagnosis <- dotplot_topics(obj = obj,
                                topic_columns = custom_feature_names,
                                group.by = "Diagnosis",
                                alpha.threshold = 0.01,
                                order = F,
  ) + labs(y = "Topics", x = "Diagnosis", title = "Diagnosis and the Average Expression of Over Different Topics ") & RotatedAxis()


  g_SampleID <- dotplot_topics(obj = obj,
                               topic_columns = custom_feature_names,
                               group.by = "SampleID",
                               alpha.threshold = 0.01,
                               order = F,
  ) &
    labs(y = "Topics", x = "SampleID", title = "SampleID and the Average Expression of Over Different Topics ") &
    RotatedAxis()

  g_seurat_clusters <- dotplot_topics(obj = obj,
                                      topic_columns = custom_feature_names,
                                      group.by = "seurat_clusters",
                                      alpha.threshold = 0.01,
                                      order = F,
  ) + labs(y = "Topics", x = "seurat_clusters", title = "seurat_clusters and the Average Expression of Over Different Topics ") & RotatedAxis()
  plot_grid(g_cell_type, g_Age, g_Diagnosis, g_SampleID, g_seurat_clusters, nrow = 5, ncol = 1)
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
#' plot_Vln(seurat_obj, "path/to/plots", "seurat_clusters", K = 10, complex_case_text = "ComplexCase")
#'
#' # Generate violin plots for other variables
#' plot_Vln(seurat_obj, "path/to/plots", "Diagnosis", K = 10)
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
                     split_by = NA)
{
  de_VlnPlot <- list()
  if (!is.na(split_by)) {
    for (i in 1:K) {
      de_VlnPlot[[i]] <- VlnPlot(obj,
                                 features = paste0("k", i),
                                 pt.size = 0,
                                 group.by = group_by_param,
                                 split.by = "SampleID") +
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
           file = glue(path_to_plots, "/VlnPlot-", group_by_param, "/split_by_SampleID.pdf"),
           width = 15,
           height = K * 2.4)
  }
  for (i in 1:K) {
    de_VlnPlot[[i]] <- VlnPlot(obj,
                               features = paste0("k", i),
                               pt.size = 0,
                               group.by = group_by_param,
                               #                           split.by = "SampleID"
    ) +
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
         file = paste0(path_to_plots, "/VlnPlot-", group_by_param, ".pdf"),
         width = 15,
         height = K * 2.4)
}

generate_VlnPlot <- function(obj, path_to_plots, K = 10, split_by = NA) {
  plot_Vln(obj,
           path_to_plots,
           "seurat_clusters",
           K,
           split_by = split_by
  )  # Violin plot grouping by Seurat clusters
  plot_Vln(obj,
           path_to_plots,
           "Diagnosis",
           K,
           split_by = split_by
  )       # Violin plot grouping by Diagnosis
  plot_Vln(obj,
           path_to_plots,
           "celltype",
           K,
           split_by = split_by
  )        # Violin plot grouping by cell type
  plot_Vln(obj,
           path_to_plots,
           "SampleID",
           K,
           split_by = split_by
  )        # Violin plot grouping by SampleID
  plot_Vln(obj,
           path_to_plots,
           "Age",
           K,
           split_by = split_by
  ) # Violin plot grouping by Age


  print("Done with Violin Plots...")

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
  print("Generating Violin Plots Complex Case ...")

  plot_Vln(subset(x = obj, subset = Diagnosis == "AD"),
           path_to_plots,
           "celltype",
           K = K,
           complex_case_text = "AD") # Violin plot to only AD grouping by celltype

  plot_Vln(subset(x = obj, subset = Diagnosis == "MCI"),
           path_to_plots,
           "celltype",
           K = K,
           complex_case_text = "MCI")

  plot_Vln(subset(x = obj, subset = Diagnosis == "Young CTRL"),
           path_to_plots,
           "celltype",
           K = K,
           complex_case_text = "YoungCTRL")

  plot_Vln(subset(x = obj, subset = Diagnosis == "HA"),
           path_to_plots,
           "celltype",
           K = K,
           complex_case_text = "HA")

  plot_Vln(subset(x = obj, subset = Diagnosis == "SuperAgers"),
           path_to_plots,
           "celltype",
           K = K,
           complex_case_text = "SuperAgers")

  plot_Vln(subset(x = obj, subset = Diagnosis == c("SuperAgers", "HA", "Young CTRL")),
           path_to_plots,
           group_by_param = "celltype",
           K = K,
           complex_case_text = "SuperAgers&HA&YoungCTRL")

  plot_Vln(subset(x = obj, subset = Diagnosis == c("SuperAgers", "HA")),
           path_to_plots,
           "celltype",
           K = K,
           complex_case_text = "SuperAgers&HA")

  plot_Vln(subset(x = obj, subset = Diagnosis == c("AD", "MCI")),
           path_to_plots,
           "celltype",
           K = K,
           complex_case_text = "AD&MCI")
  print("Done with Violin Plots...")

}

get_qualitative_colors <- function(n) {
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  return(col_vector[1:n])
}

load_objects_fit_and_obj <- function(path_to_fit, path_to_obj) {
  # load objects
  print("Loading the objects - START")
  obj <- LoadH5Seurat(path_to_obj)
  fit <- readRDS(path_to_fit)
  print("Loading the objects - FINISH")

  return(list(obj = obj, fit = fit))
}

helper_visualization_topic_model <- function(obj, fit, path_to_plots) {
  load_libraries()
  font <- "Conqueror Sans Medium"


  obj <- combine_topics_and_meta_data(obj = obj, fit = fit)
  K <- as.integer(dim(fit$L)[2])
  all_K <- paste0("k", 1:K)

  additional_colors <- get_qualitative_colors(K)
  generate_DotPlot(obj = obj,
                   custom_feature_names = all_K,
                   path_to_plots = path_to_plots)
  print("Generating structure_plot...")
  print(additional_colors)
  generate_structure_plot(obj = obj,
                          path_to_plots = path_to_plots,
                          additional_colors = additional_colors,
                          K = K,
                          fit = fit)
  print("Finish structure_plot... ")

  print("Generating UMAP...")
  generate_umap(obj = obj,
                path_to_plots = path_to_plots,
                custom_feature_names = all_K,
  )
  print("Finish UMAP...")
  print("Generating VlnPlot...")

  generate_VlnPlot(obj = obj,
                   path_to_plots = path_to_plots,
                   K = K,
                   split_by = "celltype")

  generate_vln_plot_complex(obj = obj,
                            path_to_plots = path_to_plots,
                            K = K)

  print("Finish VlnPlot...")
  #print(paste0("The Run is Finished. The Plots are saved at: ", path_to_plots))

}


run_visualization_topic_model <- function(K = 10) {
  load_libraries()

  k_folder <- glue("k", K)
  obj_folder <- "objects"
  fit_name <- glue("fitted_k_", K, ".rds")
  initial_object_name <- "5_precent.h5seurat"
  plots <- "plots"
  left_path <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/SuperAgers_n=20/5_percent"

  path_to_fit <- glue("{left_path}/{k_folder}/{obj_folder}/{fit_name}")
  path_to_obj <- glue("{left_path}/{obj_folder}/{initial_object_name}")
  path_to_plots <- glue("{left_path}/{k_folder}/{plots}/")

  loads <- load_objects_fit_and_obj(path_to_fit, path_to_obj)
  obj <- loads$obj
  fit <- loads$fit

  helper_visualization_topic_model(obj, fit, path_to_plots)
}


main <- function() {
  rm(list = ls())
  start_time <- Sys.time()


  print("Starting to run the topic visualization flow ")


  run_topic_evaluation()
  run_visualization_topic_model(10)
  run_visualization_topic_model(12)
  run_visualization_topic_model(15)
  run_visualization_topic_model(20)

  source("deprecated-run_de_visualzation.R")
  run_de_visualzation(K = 10,
                      LFC_type = "vs_null",
                      remove_zeros_from_counts_matrix = TRUE)
  end_time <- Sys.time()
    elapsed_time <- end_time - start_time

  print(glue("The whole procces is finished took {elapsed_time}"))

}

main()






