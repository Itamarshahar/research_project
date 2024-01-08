
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


for(topic in names(des)){
   result <- clusterProfiler::enrichGO( gene = des[[topic]],
                                        OrgDb = "org.Hs.eg.db",
                                        keyType="SYMBOL",
                                        ont="ALL", 
                                        pvalueCutoff=0.05, 
                                        pAdjustMethod="BH", 
                                        universe = unique(score_matrix_hippo[,"gene"]))
   print(head(result@result, n = 5))
   
}
goplot(result)