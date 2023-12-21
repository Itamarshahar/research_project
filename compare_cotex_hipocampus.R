#compare between the 14 topics of hipocampus to 15 topics of cortex



calculate_corralation_by_diagosis <- function(obj, fitted, path_to_plots, correlation_method="pearson") {
  new <- select(obj@meta.data, "Diagnosis")
  original_row_names <- row.names(new)
  dummy_matrix <- dummy_cols(new, select_columns = "Diagnosis", remove_first_dummy = FALSE)
  row.names(dummy_matrix) <- original_row_names
  dummy_matrix <- dummy_matrix[,-1]
  
  correlation <- claculate_correlation(fitted$L, dummy_matrix, method = correlation_method)
  order_row = c("Diagnosis_Young CTRL", "Diagnosis_SuperAgers", "Diagnosis_HA" ,"Diagnosis_MCI", "Diagnosis_AD")
  correlation <- correlation[order_row,]
  
  rownames(correlation) <- c("Young CTRL", "SuperAgers", "HA", "MCI", "AD")
  
  k <- ncol(correlation)
  file_name <- glue("{path_to_plots}{correlation_method}_Corralation_Between_Diagnosis_in_{k}_Topics.pdf")
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




######## 
# continue here Shmuel 
####### 

setwd('/Users/shmuel/microglia/plots')

path_to_hippocampus_14 <- "/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_14.rds"
path_to_obj <-"/Users/shmuel/microglia/objects/filtered_microglia_.rds" 
path_to_plots <- "/Users/shmuel/microglia/plots/"
microglia_obj <- readRDS(path_to_obj)
hippocampus_14_fit <- readRDS(path_to_fitted)
calculate_corralation_by_diagosis(obj = microglia_obj,
                                  fitted = hippocampus_14_fit,
                                  path_to_plots = path_to_plots)










#load the fit model
cortex_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/all_microglia_topic_fit.15.RDS")
hippocampus_14 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_14.rds")


# When looking for markers, we observe that the genes identified as markers in Topic 3 of Roi are also markers in Topic 15 of our model.
roi_fit_15$F["TOP2A",]
roi_fit_15$F["MKI67",]
roi_fit_15$F["MRC1",]
roi_fit_15$F["TMEM163",]
roi_fit_15$F["HSPA1B",]
roi_fit_15$F["HSPA1A",]
roi_fit_15$F["NFKB1",]
roi_fit_15$F["CX3CR1",]
roi_fit_15$F["DAM2",]
roi_fit_15$F["GPNMB",]
fit_15$F[c("TOP2A","MKI67","MRC1","TMEM163","HSPA1B"),]













# adi function
topic_reweight_f<-function(f_matrix){
  topics_amount<-ncol(f_matrix)
  reweighted_f<-f_matrix
  for(topic_idx in 1:topics_amount){
    # create a vector with FALSE in each idx that isn't the topic_idx (put there TRUE)
    logical_vec_curr_idx<-logical(topics_amount)
    logical_vec_curr_idx[topic_idx]<-TRUE
    cur_topic_as_mat<-f_matrix[, logical_vec_curr_idx, drop=FALSE]
    all_other_topics_mat<-f_matrix[, !logical_vec_curr_idx, drop=FALSE]
    # find max per each row (hence MARGIN=1) across all other topics
    max_other_topics<-apply(all_other_topics_mat, 1, function(x){max(x, na.rm=TRUE)})
    reweighted_f[, topic_idx]<-f_matrix[, logical_vec_curr_idx] * (log1p(cur_topic_as_mat) - log1p(max_other_topics))
  }
  return(reweighted_f)
}
reweight_roi_f <- topic_reweight_f(roi_fit_15$F)
"CX3CR1" %in%  rownames(reweight_roi_f[order(-reweight_roi_f[, "k6"]),  ][1:1500, ])  


