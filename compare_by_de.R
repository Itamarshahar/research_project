################################################################################
######         compare topics by DE gene                ##########
################################################################################

# install.packages('ggcorrplot')
# install.packages('ltm')

#load the fit model
# cortex_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/cortex_all_microglia_topic_fit.15.RDS")
# hippocampus_14 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_14.rds")
# hippocampus_13 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_13.rds")
# hippocampus_12 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_12.rds")
# hippocampus_11 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_11.rds")
# hippocampus_10 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_10.rds")

#hippocampus_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds")
# hippocampus_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds")


################################################################################
#
################################################################################
find_de_genes <- function(obj, fit) {
  reweighted_f <- topic_reweight_f(fit$F)
  df <- df_pos_des(obj, reweighted_f, fit)
  des_df <- df %>% filter(z_score_log > 1 & percent_cells > 0.01)
  return(des_df)
}


extract_counts_matrix <- function(obj){
  obj_count <- obj@assays$RNA@counts
  obj_count <- t(obj_count)
  col_sums <- colSums(obj_count)
  nonzero_cols <- col_sums != 0
  obj_count <- obj_count[, nonzero_cols]
  return(obj_count)
}


################################################################################
###### check if there is commune genes in topics that we think they competable 
################################################################################
intersaction_de_between_topic <- function(score_matrix_hippo,
                                          score_matrix_cortex,
                                          path_to_plots,
                                          sum_de = c(50, 100, 200, 300, 400),
                                          scale = "NA") {

  file_name <- glue("{path_to_plots}intersaction_de_gene_hippocampus_and_cortex_scale_{scale}.pdf")
  pdf(file_name)
  for (n in sum_de) {
    cols <- unique(score_matrix_hippo$topic)
    rows <- unique(score_matrix_cortex$topic)
    intersaction_matrix <- matrix(0, nrow = length(rows), ncol = length(cols), dimnames = list(rows, cols))
    for (i in 1:length(cols)) {
      for (j in 1:length(rows)) {
        result <- sum(
          (score_matrix_hippo[score_matrix_hippo$topic == glue("k{i}"),] %>%
            arrange(desc(reweighted_gene_score)) %>%
            head(n) %>%
            pull(gene)) %in%
            (score_matrix_cortex[score_matrix_cortex$topic == glue("k{j}"),] %>%
              arrange(desc(reweighted_gene_score)) %>%
              head(n) %>%
              pull(gene))
        )
        intersaction_matrix[glue("k{j}"), glue("k{i}")] <- result / n
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


run_compare_by_de <- function(obj, hippocampus, cortex, path_to_plots, scale = "NA") {
  source("utilsDE.R")
  for (scale in c("column", "row", "NA")) {
    obj_count <- extract_counts_matrix(obj)
    intersection_gene <- (intersect(colnames(obj_count), rownames(cortex$F)))
    cortex$F <- cortex$F[rownames(cortex$F) %in% intersection_gene,]
    hippocampus$F <- hippocampus$F[rownames(hippocampus$F) %in% intersection_gene,]

    score_matrix_hippo <- find_de_genes(obj = obj,
                                 fit = hippocampus)

    score_matrix_cortex <- find_de_genes(obj = obj,
                                  fit = score_matrix_cortex)

    intersaction_de_between_topic(score_matrix_hippo = score_matrix_hippo,
                                  score_matrix_cortex = score_matrix_cortex,
                                  path_to_plots = path_to_plots,
                                  sum_de = c(50, 100, 200, 300, 400, 1000, 5000, 10000),
                                  scale = scale)
  }
}



run_compare_by_de(obj = obj,
                  hippocampus = hippocampus_15,
                  cortex = cortex_15,
                  path_to_plots = "/Volumes/habib-lab/shmuel.cohen/microglia/plots/correlation/de/",
)


# call to main flow
# obj <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_topic_model.rds")
#obj <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/filtered_microglia.rds")
hippocampus_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds")
cortex_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/cortex_all_microglia_topic_fit.15.RDS")
path_to_plots <- "/Volumes/habib-lab/shmuel.cohen/microglia/plots/correlation/de/"

run_compare_by_de(obj = obj,
                  hippocampus = hippocampus_15,
                  cortex = cortex_15,
                  path_to_plots = path_to_plots)