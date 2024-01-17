
hippocampus_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds")

score_matrix_hippo <- df_postmean_lfsr_scores(reweighted_f = topic_reweight_f_quantile_v2(f_matrix = hippocampus_15$F,upper_quantile = 0.9,lower_quantile = 0.45)
                                              ,original_f = hippocampus_15$F)



n <- 100
des <- lapply(1:15, function(i){
  score_matrix_hippo[score_matrix_hippo$topic == glue("k{i}"), ] %>%
    arrange(desc(reweighted_gene_score)) %>%
    head(n) %>%
    pull(gene)
}) %>% setNames(., paste("k", 1:15, sep = ""))

all_pathways <- list()
for(topic in names(des)){
  print(topic)
  result <- clusterProfiler::enrichGO( gene = des[[topic]],
                                        OrgDb = "org.Hs.eg.db",
                                        keyType="SYMBOL",
                                        ont="ALL", 
                                        pvalueCutoff=0.05, 
                                        pAdjustMethod="BH", 
                                        universe = unique(score_matrix_hippo[,"gene"]))
   
   result1 <- simplify(result)
   saveRDS(result1, glue("/Volumes/habib-lab/shmuel.cohen/microglia/objects/pathway_topic_{topic}.RDS"))
   
   #heatmap <-strsplit(result@result$geneID, "/") %>% setNames(., result@result$Description) %>% stack() %>% table() %>% pheatmap::pheatmap()

   #print(head(result@result, n = 5))
   
}
clusterProfiler::goplot(result)
clusterProfiler::goplot(result, ont = "ALL", sampleGeneSets = des)

saveRDS(result, "/Volumes/habib-lab/shmuel.cohen/microglia/objects/pathway_topic_15.RDS")


heatmap <-strsplit(result@result$geneID, "/") %>% setNames(., result@result$Description) %>% stack() %>% table() %>% 
  pheatmap::pheatmap(filename = glue("/Volumes/habib-lab/shmuel.cohen/microglia/plots/Plots_for_k=15/pathways/pathways_topic_{topic}"))


#heatmap <-strsplit(result@result$geneID, "/") %>% setNames(., result@result$Description) %>% stack() %>% table() %>% pheatmap::pheatmap()

library("grid")
library("ggplotify")

# Initialize an empty list to store grobs
heatmaplist <- list()

for (topic in names(des)) {
  path <- glue("/Volumes/habib-lab/shmuel.cohen/microglia/objects/pathway_topic_{topic}.RDS")
  
  # Skip certain topics
  if (topic %in% c("k13", "k9")) {
    next
  }
  
  print(glue("/Volumes/habib-lab/shmuel.cohen/microglia/objects/pathway_topic_{topic}.RDS"))
  
  pathway_result <- readRDS(path)
  
  # Create a heatmap using pheatmap
  heatmap_data <- strsplit(pathway_result@result$geneID, "/") %>% 
    setNames(., pathway_result@result$Description) %>% 
    stack() %>% table()
  
  heatmap_plot <- pheatmap::pheatmap(heatmap_data)
  
  # Convert the heatmap to grob
  heatmap_grob <- ggplotify::as.grob(heatmap_plot)
  
  # Store the grob in the list
  heatmaplist[[topic]] <- heatmap_grob
}


# Save the arranged plot
ggsave("/Volumes/habib-lab/shmuel.cohen/microglia/plots/Plots_for_k=15/pathways/pathways_heatmaps.pdf", width = 14, height = 10)

ggsave("/Volumes/habib-lab/shmuel.cohen/microglia/plots/Plots_for_k=15/pathways/pathways_heatmaps.pdf", width = 14, height = 10)
  
  
