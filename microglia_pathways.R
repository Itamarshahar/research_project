########



generate_k_to_genes_table <- function(de_table){
  k_to_genes_table <- lapply(1:15, function(i){
  de_table[de_table$topic == glue::glue("k{i}"), ] %>%
    pull(gene)
}) %>% setNames(., paste("k", 1:15, sep = ""))
}


generate_de_objects <- function(des, all_genes){
  for(topic in names(des)){
  print(topic)
  result <- clusterProfiler::enrichGO(
                                        gene = des[[topic]],
                                        #gene = des[["k1"]],
                                        OrgDb = "org.Hs.eg.db",
                                        keyType="SYMBOL",
                                        ont="ALL",
                                        pvalueCutoff=0.05,
                                        pAdjustMethod="BH",
                                        universe = unique(all_genes))

   result <- clusterProfiler::simplify(result)
   saveRDS(result, glue("/Volumes/habib-lab/shmuel.cohen/microglia/objects/pathway_topic_{topic}.RDS"))
}
}



generate_de_plots <- function(k_to_genes_table, path_to_objects, path_to_plot="/Volumes/habib-lab/shmuel.cohen/microglia/plots/Plots_for_k=15/pathways/pathways_heatmaps.pdf") {
path_to_plot <- "/Volumes/habib-lab/shmuel.cohen/microglia/plots/Plots_for_k=15/pathways/pathways_heatmaps.pdf"
pdf(path_to_plot, width = 20, height = 8)
for (topic in names(k_to_genes_table)){
  # if (topic %in% c("k13", "k9")) {
  #   next
  # } <-
  path <- file.path(path_to_objects, glue("pathway_topic_", topic, ".RDS"))
  pathway_result <- readRDS(path)
  heatmap <-strsplit(pathway_result@result$geneID, "/") %>% setNames(., pathway_result@result$Description) %>%
    stack() %>% table() %>%
    pheatmap::pheatmap(main = glue("pathway_topic_{topic}"),
                                   fontsize_col = 5,
                                  fontsize_row = 5,
                                   angle_col = 45,
                                    cluster_rows = T,
                                      legend = F,
                                   )
  print(heatmap)

}
dev.off()

}



load_libraries <- function(){
  library(dplyr)
  library(tidyverse)
  library(glue)
}

run_pathways <- function(path_to_fit="/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds",
                         path_to_de="/Volumes/habib-lab/shmuel.cohen/microglia/objects/DE_hippocampus_X15.rds",
                         path_to_plot="/Volumes/habib-lab/shmuel.cohen/microglia/plots/Plots_for_k=15",
                          path_to_objects="/Volumes/habib-lab/shmuel.cohen/microglia/objects/",
                         generate_objects = TRUE
                          ) {
      load_libraries()
      fit <- readRDS(path_to_fit)
      de_table <- readRDS(path_to_de)
      all_genes <- rownames(fit$F)
      k_to_genes_table <- generate_k_to_genes_table(de_table)
      if (generate_objects){
        generate_de_objects(k_to_genes_table, all_genes)
      }
      generate_de_plots(k_to_genes_table, path_to_objects, path_to_plot)
}