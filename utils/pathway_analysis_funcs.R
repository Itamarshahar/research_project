##########################################
######## source files & libraries ########
##########################################
setwd("~/HFD/")
library(dplyr)
library(ggplot2)
library(reshape)
library(clusterProfiler)
source("utils/pathways.analysis.gilad.R")

# Based on Gali's function.
# Gets a data frame of de_genes (output of findMarkers) and outputs its pathway.
# Universe - all genes expressed by the cells we are looking at (can take from counts \ data unique(Dimnames[[1]])).
get_pathways_go_bp<-function(de_genes, avg_log2FC_threshold=0.5, p_val_adj_threshold=0.05, universe){
  chosen_de_genes_df<-create_filtered_de_df(de_genes, avg_log2FC_threshold, p_val_adj_threshold)
  universe_ids<-GeneIdMapping("mmu")$ids[universe]
  res<-enrichGO(gene=chosen_de_genes_df$id, universe=universe_ids, OrgDb=org.Mm.eg.db, ont="BP",
  pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE)
  return(res@result)
}

# Based on Roi's function.
# Gets a data frame of de_genes (output of findMarkers) and outputs its pathway.
# Universe - all genes expressed by the cells we are looking at (can take from counts \ data unique(Dimnames[[1]])).
enrichment_analysis_wrapper<-function(de_genes, avg_log2FC_threshold=0.5, p_val_adj_threshold=0.05, qvalue_threshold=0.05, universe){
  chosen_de_genes_df<-create_filtered_de_df(de_genes, avg_log2FC_threshold, p_val_adj_threshold)
  universe_ids<-GeneIdMapping("mmu")$ids[universe]
  results_list<-EnrichmentAnalysis(chosen_de_genes_df, formula="id~cluster", universe=universe_ids, species="mmu")
  # adds a row with db name to the results_list
  for(db_name in names(results_list)){
    results_list[[db_name]]@compareClusterResult['db']<-db_name 
  }
  #concatenate all objects of results_list into 1 table
  results<-do.call("rbind", lapply(results_list, attr, 'compareClusterResult')) %>%
    filter(qvalue < qvalue_threshold) %>% select(-c(desc_id, Cluster))
  return(results)
}

  
# Helper function for get_pathways function.
# Filters the data frame of de_genes according to the 2 given thresholds 
# and adds 2 columns - cluster name and id mapping.
create_filtered_de_df<-function(de_genes, avg_log2FC_threshold, p_val_adj_threshold){
  de_genes<-de_genes %>% mutate(gene=as.character(row.names(de_genes))) %>%
    filter(avg_log2FC > avg_log2FC_threshold & p_val_adj < p_val_adj_threshold) %>%
    mutate(cluster="test", id=GeneIdMapping("mmu")$ids[gene])
  return(de_genes)
}

# Visualize significant results using a barplot.
# X axis - (-log(qvalue)), we want it to be the highest.
# Y axis - the different pathways.
pathways_barplots<-function(enrichment_results, qvalue_threshold=0.05){
  enrichment_results %>% filter(qvalue <= qvalue_threshold) %>%
  ggplot(aes(x=reorder(Description, -1 * log10(qvalue)), y=-1 * log10(qvalue))) + geom_col() + coord_flip() +
    labs(y="-log10(q value)", x="Gene Set") + theme_bw() + theme(panel.border = element_blank(), text=element_text(size=15), axis.text=element_text(size=16),
                                                                 axis.title.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }