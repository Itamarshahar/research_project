# run example:
# de_hippocampus <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/DE_hippocampus_X15.rds")
# de_cortex <- read.csv("/Volumes/habib-lab/shmuel.cohen/microglia/objects/DE_cortex_X15.csv")
# fit_hippo <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds")
# fit_cortex <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/cortex_all_microglia_topic_fit.15.RDS")
# matrix <- generate_jaccrd_matrix(de_hippocampus, de_cortex)
# matrix <- generate_hypergeomtric_matrix (de_hippocampus, de_cortex, fit_hippo, fit_cortex)


generate_jaccrd_matrix <- function(de_hippocampus, de_cortex){
  extract_des <- function(de_table){
    des <- lapply(1:15, function(i){
      de_table[de_table$topic == glue("k{i}"), ] %>%
        pull(gene)
    }) %>% setNames(., paste("k", 1:15, sep = ""))
    return(des)
  }
  
  des_hippo <- extract_des(de_hippocampus) 
  des_cortex <- extract_des(de_cortex) 
  
  # Compute the Jaccard similarity matrix
  similarity_matrix <- matrix(0, nrow = length(des_hippo), ncol = length(des_cortex),dimnames = list(names(des_hippo), names(des_cortex)))
  for (i in seq_along(des_hippo)) {
    for (j in seq_along(des_cortex)) {
      sets = list(
        des_hippo[[i]],
        des_cortex[[j]]
      )
      similarity_matrix[i, j] <- mlr3measures::jaccard(sets)
    }
  }
  return(similarity_matrix)
}



#calculate hypergeomtric test
generate_hypergeomtric_matrix <- function(de_hippocampus, de_cortex, fit_hippo, fit_cortex){
  extract_des <- function(de_table){
    des <- lapply(1:15, function(i){
      de_table[de_table$topic == glue("k{i}"), ] %>%
        pull(gene)
    }) %>% setNames(., paste("k", 1:15, sep = ""))
    return(des)
  }
  
  des_hippo <- extract_des(de_hippocampus) 
  des_cortex <- extract_des(de_cortex) 
  
  #calculate the bachground genes
  #intersection_gene <- (intersect(rownames(fit_hippo$F), rownames(fit_cortex$F)))
  des_gene <- list()
  for (i in seq_along(des_hippo)){
    des_gene <- c(des_gene,des_hippo[[i]])
  }
  for (i in seq_along(des_cortex)){
    des_gene <- c(des_gene,des_cortex[[i]])
  }
  des_gene <- unique(des_gene)
  #calculate the matrix
  sum_gene <- length(unlist(des_gene))
  hypergeomtric_matrix <- matrix(0, nrow = length(des_hippo), ncol = length(des_cortex), dimnames = list(names(des_hippo), names(des_cortex)))
  
  for (i in seq_along(des_hippo)) {
    for (j in seq_along(des_cortex)) {
      x <- length(intersect(des_hippo[[i]], des_cortex[[j]]))
      m <- length(des_cortex[[j]])
      n <- sum_gene - m
      k <- length(des_hippo[[i]])
      hypergeomtric_matrix[i, j] <- phyper(q = x, m = m, n = n, k = k, lower.tail = FALSE, log.p = FALSE)
    }
  }
  return(hypergeomtric_matrix)
}




save_to_pdf <- function(hypergeomtric_matrix, matrix){
  file_name <- glue("/Users/shmuel/microglia/plots/gene_correlation/correlation_de_adi/p_value.pdf")
  pdf(file_name)
  pheatmap(t(hypergeomtric_matrix),cluster_rows = F, cluster_cols = F,
           cell_fun = HeatmapHelper_add_values_to_display(correlation = t(hypergeomtric_matrix),OnlyPositive = TRUE))
  dev.off()
  
  file_name <- glue("/Users/shmuel/microglia/plots/gene_correlation/correlation_de_adi/jaccard_marix.pdf")
  pdf(file_name)
  pheatmap(t(matrix) ,cluster_rows = F, cluster_cols = F,cell_fun = HeatmapHelper_add_values_to_display(correlation = t(matrix),OnlyPositive = TRUE))
  dev.off()
}




