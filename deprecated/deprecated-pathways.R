################################################################################
######      calculate the pathways of every topics              ##########
################################################################################


hippocampus_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds")

#choose the 100 top DEs gene by Adi function
score_matrix_hippo <- df_postmean_lfsr_scores(reweighted_f = topic_reweight_f_quantile_v2(f_matrix = hippocampus_15$F,upper_quantile = 0.9,lower_quantile = 0.45)
                                              ,original_f = hippocampus_15$F)



n <- 100
des <- lapply(1:15, function(i){
  score_matrix_hippo[score_matrix_hippo$topic == glue("k{i}"), ] %>%
    arrange(desc(reweighted_gene_score)) %>%
    head(n) %>%
    pull(gene)
}) %>% setNames(., paste("k", 1:15, sep = ""))

#run the pathways and save as rds objects
all_pathways <- list()
for(topic in names(des)){
  result <- clusterProfiler::enrichGO( gene = des[[topic]],
                                        OrgDb = "org.Hs.eg.db",
                                        keyType="SYMBOL",
                                        ont="ALL", 
                                        pvalueCutoff=0.05, 
                                        pAdjustMethod="BH", 
                                        universe = unique(score_matrix_hippo[,"gene"]))
   
   result <- simplify(result)
   saveRDS(result, glue("/Volumes/habib-lab/shmuel.cohen/microglia/objects/pathway_topic_{topic}.RDS"))

}

#do the plot
#clusterProfiler::goplot(result)
#clusterProfiler::goplot(result, ont = "ALL", sampleGeneSets = des)

# read the object and print heatmap of pathays with genes
file_name <- glue("/Volumes/habib-lab/shmuel.cohen/microglia/plots/Plots_for_k=15/pathways/pathways_heatmaps.pdf")
pdf(file_name, width = 20, height = 8)
for (topic in names(des)){
  if (topic %in% c("k13", "k9")) {
    next
  }
  path <- glue("/Volumes/habib-lab/shmuel.cohen/microglia/objects/pathway_topic_{topic}.RDS")
  pathway_result <- readRDS(path)
  heatmap <-strsplit(pathway_result@result$geneID, "/") %>% setNames(., pathway_result@result$Description) %>%
    stack() %>% table() %>%
    pheatmap::pheatmap(main = glue("pathway_topic_{topic}"),
                                   fontsize_col = 4,
                                  fontsize_row = 4,
                                   angle_col = 45,
                                    cluster_rows = F,
                                      legend = F,
                                   )
  print(heatmap)

}
dev.off()






