# Run the topic evaluations flow:
load_libraries <- function() {
  library(Matrix) #-c r r-matrix conda-forge r-ggplot2 conda-forge r-cowplot conda-forge r-seuratdisk conda-forge r-glue conda-forge r-magrittr
  library(fastTopics)
  library(ggplot2)
  library(cowplot)
  library(RColorBrewer)
  library(hexbin)
  library(viridis)
  library(Seurat)
  library(SeuratDisk)
  # library(SeuratData)
  library(magrittr)
  library(dplyr)
  library(gridExtra)
  library(ggplot2)
  library(grid)
  library(viridisLite)
  # library(ggpointdensity)
  library(reshape2)
  library(glue)
  source("utils.R")
  library(ComplexHeatmap)
  library(circlize)
  # library(combinat)
  #library(glue)
  library(stringr)
}


generate_all_permutations <- function(lst) {
  # Initialize an empty vector to store the permutations
  all_permutations <- character(0)
  
  # Generate all permutations 
  for (i in 1:length(lst)) {
    for (j in i:length(lst)) {
      all_permutations <- c(all_permutations, paste(lst[i], lst[j], sep = "_"))
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
    fits_list[[as.character(obj$k)]] <- obj$fit  # Use as.character to ensure k is a character key
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
                 column_title_gp = gpar(fontsize = 10, fontface = "bold"),
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

}



run_topic_evaluation <- function(fit_list, path_to_plots, correlation_method="pearson") {
  print("Loading libraries")
  load_libraries()
  helper_topic_evaluation(fit_list , glue(path_to_plots, "gene_correlation/"), type = "genes", correlation_method = correlation_method) #, correlation_method = "kendall")
  helper_topic_evaluation(fit_list, glue(path_to_plots, "cell_correlation/"), type = "cells" , correlation_method = correlation_method)
}



