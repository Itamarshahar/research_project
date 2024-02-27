
################################################################################
######calculate the correlation between the diagnosis to topics ##########
################################################################################
calculate_correlation_by_diagosis.percell <- function(obj, fitted, path_to_plots, correlation_method="pearson") {
  
  new <- dplyr::select(obj@meta.data, "Diagnosis" )
  original_row_names <- row.names(new)
  dummy_matrix <- fastDummies::dummy_cols(new, select_columns = "Diagnosis", remove_first_dummy = FALSE)
  row.names(dummy_matrix) <- original_row_names
  dummy_matrix <- dummy_matrix[,-1]
  
  correlation <- claculate_correlation(fitted$L, dummy_matrix, method = correlation_method)
  order_row = c("Diagnosis_Young CTRL", "Diagnosis_SuperAgers", "Diagnosis_HA" ,"Diagnosis_MCI", "Diagnosis_AD")
  correlation <- correlation[order_row,]
  
  rownames(correlation) <- c("Young", "SuperAgers", "HA", "MCI", "AD")
  col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  col_fun(seq(-20, 20))
  k <- ncol(correlation)
  file_name <- glue("{path_to_plots}{correlation_method}_Corralation_Between_Diagnosis_in_{k}_Topics_percell.pdf")
  pdf(file_name)
  
  draw(Heatmap(correlation,
          cluster_rows = T,
          cluster_columns = T,
          col = col_fun,
          column_title = glue("The Corralation Between Topic and Diagnosis"),
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
          name = "Correlation",
          rect_gp = gpar(col = "white", lwd = 2),
          column_names_rot = 45,
          cell_fun = HeatmapHelper_add_values_to_display(correlation = correlation,
                                                      OnlyPositive = TRUE)
  ))
  
  dev.off()
  
}

# Itamar fix it 
calculate_correlation_by_diagosis.persampleID <- function(obj, fitted, path_to_plots, correlation_method="pearson") {
  
  new <- dplyr::select(obj@meta.data, "SampleID_Diagnosis" )
  original_row_names <- row.names(new)
  dummy_matrix <- fastDummies::dummy_cols(new, select_columns = "SampleID_Diagnosis", remove_first_dummy = FALSE)
  row.names(dummy_matrix) <- original_row_names
  dummy_matrix <- dummy_matrix[,-1]
  
  correlation <- claculate_correlation(fitted$L, dummy_matrix, method = correlation_method)

  # Rename rows
  col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  col_fun(seq(-20, 20))
  k <- ncol(correlation)
  file_name <- glue("{path_to_plots}_Corralation_Diagnosis_in_{k}_Topics_SampleID.pdf")
  pdf(file_name)
  
  draw(Heatmap(correlation,
               cluster_rows = T,
               cluster_columns = T,
               col = col_fun,
               column_title = glue("The Corralation Between Topic and Diagnosis"),
               column_title_gp = gpar(fontsize = 10, fontface = "bold"),
               name = "Correlation",
               rect_gp = gpar(col = "white", lwd = 2),
               column_names_rot = 45,
               cell_fun = HeatmapHelper_add_values_to_display(correlation = correlation,
                                                              OnlyPositive = TRUE)
  ))
  
  dev.off()
  
}

#calculate the correlation between the diagnosis to topics 
calculate_correlation_by_diagosis.persample_good_maybe <- function(obj, fitted, path_to_plots, correlation_method="pearson") {
  
  new <- dplyr::select(obj@meta.data, "SampleID" )
  
  mean.of.topic.persample <- do.call(rbind, sapply(unique(obj@meta.data$SampleID), function(sampleID){
    cells.of.sample <- rownames(obj@meta.data)[obj@meta.data$SampleID == sampleID]
    colMeans(fitted$L[cells.of.sample,])
  }, USE.NAMES = T, simplify = F))
  
  dummy_matrix <- fastDummies::dummy_cols(obj@meta.data[,c("Diagnosis", "SampleID")] %>% unique() , select_columns = "Diagnosis", remove_first_dummy = FALSE)
  
  correlation <- claculate_correlation(mean.of.topic.persample, dummy_matrix, method = correlation_method)
  #order_row = c("Diagnosis_Young CTRL", "Diagnosis_SuperAgers", "Diagnosis_HA" ,"Diagnosis_MCI", "Diagnosis_AD")
  #correlation <- correlation[order_row,]
  
  #rownames(correlation) <- c("Young CTRL", "SuperAgers", "HA", "MCI", "AD")
  col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  col_fun(seq(-20, 20))
  k <- ncol(correlation)
  file_name <- glue("{path_to_plots}{correlation_method}_Corralation_Between_Diagnosis_in_{k}_Topics_persample.pdf")
  pdf(file_name)
  
  draw(Heatmap(correlation,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               col = col_fun,
               column_title = glue("The Corralation Between Diagnosis in {k} Topics"),
               column_title_gp = gpar(fontsize = 10, fontface = "bold"),
               name = "Correlation",
               rect_gp = gpar(col = "white", lwd = 2),
               column_names_rot = 45,
               cell_fun = HeatmapHelper_add_values_to_display(correlation = correlation,
                                                              OnlyPositive = TRUE)
  ))
  
  dev.off()
  
}



#path_to_obj <-"/Users/shmuel/microglia/objects/filtered_microglia_.rds" 
path_to_plots <- "/Volumes/habib-lab/shmuel.cohen/microglia/plots/correlation/"
#obj <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/filtered_microglia.rds")
calculate_correlation_by_diagosis.persample(obj = obj,
                                  fitted = hippocampus_15,
                                  path_to_plots = path_to_plots)
calculate_correlation_by_diagosis.percell(obj = obj,
                                            fitted = hippocampus_15,
                                            path_to_plots = path_to_plots)

calculate_correlation_by_diagosis.persampleID(obj = obj,
                                              fitted = hippocampus_15,
                                              path_to_plots = path_to_plots)



calculate_correlation_by_diagosis.persampleID(obj = obj,
                                              fitted = hippocampus_15,
                                              path_to_plots = path_to_plots)
