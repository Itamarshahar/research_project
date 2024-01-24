# rm(list = ls
rm(list = ls())

load_libraries <- function() {
  library(ComplexHeatmap)
  library(circlize)
  source("utils.R")
  library(glue)
  library(combinat)
}

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

run_main_flow <- function(fits_list, path_to_plots, type = "cells", correlation_method = "pearson") {
  print("Running main flow")
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
      correlation <- claculate_correlation(fit_k_left$L, fit_k_right$L,method=correlation_method)
      file_name <- glue('{correlation_method}_corrlation_between_{type}_')
    }
    else {
      correlation <- claculate_correlation(fit_k_left$F, fit_k_right$F, method=correlation_method)
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

run_topic_evaluation <- function(fit_list=NA, local_run=TRUE) {
  print("Loading libraries")
  load_libraries()
  if (local_run) {
    obj_folder <- "objects"
  fit_name_k_10 <- glue("fitted_k_", 10, ".rds")
  fit_name_k_15 <- glue("fitted_k_", 15, ".rds")
  fit_name_k_12 <- glue("fitted_k_", 12, ".rds")
  fit_name_k_20<- glue("fitted_k_", 20, ".rds")
  initial_object_name <- "five_prec.h5seurat"
  plots <- "plots"
  left_path <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/SuperAgers_n=20/5_percent"

  path_to_plots <- glue("{left_path}/{plots}/")
  #path_to_fit_k_10 <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/SuperAgers_n=20/5_percent/k10/objects/fitted_k_10.rds"
  path_to_fit_k_10 <-  glue("{left_path}/k10/{obj_folder}/{fit_name_k_10}")
  path_to_fit_k_15 <-  glue("{left_path}/k15/{obj_folder}/{fit_name_k_15}")
  path_to_fit_k_12 <-  glue("{left_path}/k12/{obj_folder}/{fit_name_k_12}")
  path_to_fit_k_20 <-  glue("{left_path}/k20/{obj_folder}/{fit_name_k_20}")

  }

  #path_to_fit_k_15 <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/SuperAgers_n=20/5_percent/k15/objects/fitted_k_15.rds"
  run_main_flow(c(path_to_fit_k_10, path_to_fit_k_15,path_to_fit_k_12, path_to_fit_k_20), path_to_plots, type = "genes")#, correlation_method = "kendall")
  run_main_flow(c(path_to_fit_k_10, path_to_fit_k_15,path_to_fit_k_12, path_to_fit_k_20), path_to_plots, type = "cells")#, correlation_method = "kendall")
}





