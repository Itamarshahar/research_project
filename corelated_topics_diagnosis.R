
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

################################################################################
######         compare topics by DE gene                ##########
################################################################################




#load the fit model
cortex_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/cortex_all_microglia_topic_fit.15.RDS")
hippocampus_14 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_14.rds")
hippocampus_13 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_13.rds")
hippocampus_12 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_12.rds")
hippocampus_11 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_11.rds")
hippocampus_10 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_10.rds")

#hippocampus_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds")
hippocampus_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds")


################################################################################
######                adi function                          ##########
################################################################################

#topic_reweight_f_quantile_v2<-function(f_matrix, upper_quantile, lower_quantile){
  topics_amount<-ncol(f_matrix)
  reweighted_f<-f_matrix
  for(topic_idx in 1:topics_amount){
    # create a vector with FALSE in each idx that isn't the topic_idx (put there TRUE)
    logical_vec_curr_idx<-logical(topics_amount)
    logical_vec_curr_idx[topic_idx]<-TRUE
    cur_topic_as_mat<-f_matrix[, logical_vec_curr_idx, drop=FALSE]
    all_other_topics_mat<-f_matrix[, !logical_vec_curr_idx, drop=FALSE]
    # find quantile per each row (hence MARGIN=1) => across all other topics
    upper_quantile_other_topics<-apply(all_other_topics_mat, 1, function(x){quantile(x, upper_quantile)})
    lower_quantile_other_topics<-apply(all_other_topics_mat, 1, function(x){quantile(x, lower_quantile)})
    log_ratio_ours_vs_upper<-log(cur_topic_as_mat+1e-15) - log(upper_quantile_other_topics+1e-15)
    log_ratio_ours_vs_lower<-log(cur_topic_as_mat+1e-15) - log(lower_quantile_other_topics+1e-15)
    reweighted_f[, topic_idx]<-f_matrix[, logical_vec_curr_idx] * log_ratio_ours_vs_upper *
      ifelse(log_ratio_ours_vs_lower<0, Inf, log_ratio_ours_vs_lower^sign(log_ratio_ours_vs_upper))
  }
  return(reweighted_f)
}
df_postmean_lfsr_scores<-function(reweighted_f, original_f){
  reweighted_gene_score<-reweighted_f %>%
    melt(variable.factor = FALSE) %>%
    rename_with(~sub("Var1", "gene" , .x)) %>%
    rename_with(~sub("Var2", "topic" , .x)) %>%
    rename_with(~sub("value", "reweighted_gene_score" , .x))
  gene_scores <- original_f %>%
    as.data.frame %>%
    rownames_to_column('gene') %>%
    pivot_longer(!gene, names_to='topic', values_to='gene_score') %>%
    dplyr::select(gene, topic, gene_score)
  # create df of gene name, topic, gene id,  reweighted_gene_score & original gene score
  de_gene<-reweighted_gene_score %>%
    dplyr::left_join(gene_scores, by = c("gene", "topic"))
  return(de_gene)
}

# I edited a small bug in the reweighing function, so in case you'll want to use it -
topic_reweight_f_quantile_v2<-function(f_matrix, upper_quantile, lower_quantile){
  topics_amount<-ncol(f_matrix)
  reweighted_f<-f_matrix
  for(topic_idx in 1:topics_amount){
    # create a vector with FALSE in each idx that isn't the topic_idx (put there TRUE)
    logical_vec_curr_idx<-logical(topics_amount)
    logical_vec_curr_idx[topic_idx]<-TRUE
    cur_topic_as_mat<-f_matrix[, logical_vec_curr_idx, drop=FALSE]
    all_other_topics_mat<-f_matrix[, !logical_vec_curr_idx, drop=FALSE]
    # find quantile per each row (hence MARGIN=1) => across all other topics
    upper_quantile_other_topics<-apply(all_other_topics_mat, 1, function(x){quantile(x, upper_quantile)})
    lower_quantile_other_topics<-apply(all_other_topics_mat, 1, function(x){quantile(x, lower_quantile)})
    log_ratio_ours_vs_upper<-log(cur_topic_as_mat+1e-15) - log(upper_quantile_other_topics+1e-15)
    log_ratio_ours_vs_lower<-log(cur_topic_as_mat+1e-15) - log(lower_quantile_other_topics+1e-15)
    reweighted_f[, topic_idx]<- log_ratio_ours_vs_upper *
      ifelse(log_ratio_ours_vs_lower<0, Inf, (f_matrix[, logical_vec_curr_idx]*log_ratio_ours_vs_lower)^sign(log_ratio_ours_vs_upper))
  }
  return(reweighted_f)
}


################################################################################
###### check if there is commune genes in topics that we think they competable 
################################################################################
intersaction_de_between_topic <- function(score_matrix_hippo,
                                          score_matrix_cprtex, 
                                          path_to_plots, 
                                          sum_de=c(50,100,200,300,400), 
                                          scale="NA"){
  file_name <- glue("{path_to_plots}intersaction_de_gene_hippocampus_and_cortex_scale_{scale}.pdf")
  pdf(file_name)
  for (n in sum_de){
    cols <- unique(score_matrix_hippo$topic)
    rows <- unique(score_matrix_cprtex$topic)
    intersaction_matrix <- matrix(0, nrow = length(rows), ncol = length(cols), dimnames = list(rows, cols))
    for (i in 1:length(cols)) {
      for (j in 1:length(rows)){
        result <- sum(
          (score_matrix_hippo[score_matrix_hippo$topic == glue("k{i}"), ] %>%
             arrange(desc(reweighted_gene_score)) %>%
             head(n) %>%
             pull(gene)) %in%
            (score_matrix_cprtex[score_matrix_cprtex$topic == glue("k{j}"), ] %>%
               arrange(desc(reweighted_gene_score)) %>%
               head(n) %>%
               pull(gene))
        )
        intersaction_matrix[glue("k{j}"),glue("k{i}")] <- result/n
      }
    }
    
    
    draw(pheatmap(intersaction_matrix,
                  scale = scale,
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 #col = col_fun,
                 column_title = glue("sum intersaction from {n} de genes hippocampus and cortex"),
                 column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 name = glue("intersection size from{n}"),
                 #rect_gp = gpar(col = "white", lwd = 2),
                 #column_names_rot = 45,
                 cell_fun = HeatmapHelper_add_values_to_display(correlation = intersaction_matrix,
                                                                OnlyPositive = TRUE)
    ))
  }
  dev.off()
  
}


main_flow <- function(obj, hippocampus, cortex, path_to_plots, scale ="NA"){
  
  obj_count <- extract_counts_matrix(obj)
  intersection_gene <- (intersect(colnames(obj_count), rownames(cortex$F)))
  cortex$F <- cortex$F[rownames(cortex$F) %in% intersection_gene,]
  hippocampus$F <- hippocampus$F[rownames(hippocampus$F) %in% intersection_gene,]
  
  score_matrix_hippo <- df_postmean_lfsr_scores(reweighted_f = topic_reweight_f_quantile_v2(f_matrix = hippocampus$F,upper_quantile = 0.9,lower_quantile = 0.45)
                                                ,original_f = hippocampus$F)
  score_matrix_cprtex <- df_postmean_lfsr_scores(reweighted_f = topic_reweight_f_quantile_v2(f_matrix = cortex$F,upper_quantile = 0.9,lower_quantile = 0.45)
                                                 ,original_f = cortex$F)
  
  
  intersaction_de_between_topic(score_matrix_hippo = score_matrix_hippo,
                                score_matrix_cprtex = score_matrix_cprtex,
                                path_to_plots = path_to_plots,
                                sum_de=c(50,100,200,300,400,1000,5000,10000),
                                scale = scale)
}

for (name in c("column", "row", "NA") ){
  main_flow(obj= obj,
            hippocampus =  hippocampus_15,
            cortex = cortex_15,
            #path_to_plots = "/Users/shmuel/microglia/plots/gene_correlation/correlation_de_adi_15/",
            path_to_plots = "/Volumes/habib-lab/shmuel.cohen/microglia/plots/correlation/de/",
            scale = name)  
  
}




