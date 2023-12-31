# function to calculate the correlation between topic modeling from different data bases

# 6 -  this is work!
######################################################################################################
### Predict the L matrix based on the F matrix and compute the correlation with the cells matrix. ###
###################################################################################################### 
#obj <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/filtered_microglia.rds")


correlation_with_cortex <- function(obj, cortex_fit_15, fits_list, path_to_plots, type = "cells", correlation_method = "pearson", predicted="NA"){
  #extract the counts matrix 
  obj_count <- obj@assays$RNA@counts
  obj_count <- t(obj_count)
  col_sums <- colSums(obj_count)
  nonzero_cols <- col_sums != 0
  obj_count <- obj_count[, nonzero_cols]
  
  #intersection to the commonly gene
  intersection_gene <- (intersect(colnames(obj_count), rownames(cortex_fit_15$F)))
  cortex_fit_15$F <- cortex_fit_15$F[rownames(cortex_fit_15$F) %in% intersection_gene,]
  obj_count <- obj_count[,colnames(obj_count) %in% intersection_gene]
  
  #predict the L matrix by the F of 500 and count of 18
  if (all(!is.na(predicted))) {
    predicted <- predict(cortex_fit_15, obj_count, numiter = 100)
  }
  
  # print the heatmap of correlations
  helper_topic_evaluation(fits_list=fit_files_paths, path_to_plots="/Users/shmuel/microglia/plots/gene_correlation/vs_15topics_of_500/", type = "cells", correlation_method = "pearson", L500=predicted)
  
}



# print the heatmap of correlations
helper_topic_evaluation <- function(fits_list, path_to_plots, type = "cells", correlation_method = "pearson", L500="NA") {
  extract_k <- function(path_to_fit) {
    fit <- readRDS(path_to_fit)
    cat(dim(fit$F)[2], "||")
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
  generate_all_permutations <- function(lst) {
    # Initialize an empty vector to store the permutations
    all_permutations <- character(0)
    
    # Generate all permutations 
    for (i in 1:length(lst)) {
      all_permutations <- c(all_permutations, paste(lst[i]))
    }
    return(all_permutations)
  }
  print(glue("Running topic evaluation flow for {type} and {correlation_method}"))
  col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  col_fun(seq(-20, 20))
  
  fits_list <- generate_fits_list(fits_list)
  all_permutations <- generate_all_permutations(names(fits_list))

  k_right <- "L500"
  fit_k_right <- L500
  for (per in all_permutations) {
    cat("per = ", per, ".")
    split_parts <- unlist(strsplit(per," "))
    k_left <- split_parts[1]
    fit_k_left <- fits_list[[k_left]]
    if (type == "cells") {
      #cat(dim(fit_k_left$L), "/n")
      correlation <- claculate_correlation(fit_k_left$L, fit_k_right, method = correlation_method)
      file_name <- glue('{correlation_method}_corrlation_between_{type}_')
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

