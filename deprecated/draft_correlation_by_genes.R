# Here, I want to explore several approaches to address the challenge of comparing genes.
# The issue arises when comparing topics from different cells or when comparing the same cells, but with a focus on gene expression.
# We can't simply check the correlation between the gene scores for each topic because there is no meaningful interpretation 
# of gene scores without considering the context. Therefore, we cannot straightforwardly assess the correlation between the topic vectors.

# First, let's print the correlation between genes and see what we get. Of course, it may not be good enough.
hippocampus_fit_14 <-  readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_14.rds")
cortex_fit_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/all_microglia_topic_fit.15.RDS")
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
col_fun(seq(-20, 20))

correlation_method <- "pearson"
correlation <- claculate_correlation(hippocampus_fit_14$F, cortex_fit_15$F, method = correlation_method)
  file_name <- glue('{correlation_method}_corrlation_between_hippocampus_to_cortex_by_gene')
  pdf(glue("/Users/shmuel/microglia/plots/gene_correlation/correlation_with_cortex/{file_name}.pdf"))
  draw(Heatmap(correlation,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               col = col_fun,
               column_title = glue("The Corralation Between hippocampus(k=14) with cortex(k=15) Topics"),
               column_title_gp = gpar(fontsize = 10, fontface = "bold"),
               name = "Correlation",
               rect_gp = gpar(col = "white", lwd = 2),
               column_names_rot = 45,
               cell_fun = HeatmapHelper_add_values_to_display(correlation = correlation,
                                                              OnlyPositive = TRUE)
  )
  )
  dev.off()
  load_libraries()


# 1
######################################################################################################
########## Focus on a specific topic, the split of this topic, and compare the markers. ##############
######################################################################################################

########################################################
# here we want check the split of topic 6_8 to 6_9 & 4_9
########################################################
fit <-  readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_9.rds")
obj <- LoadH5Seurat("/Users/shmuel/microglia/objects/filtered_microglia.h5Seurat")

#scater plot to show the gene that make the different between the topic
pdf("tmp.pdf", height = 6, width = 6)
ggplot(cbind(y = sort(fit$F[,c(6)]), x = 1:nrow(fit$F)) %>% as.data.frame.matrix() %>% filter(y > 0.00005), aes(x, y)) + geom_point() + geom_hline(yintercept = 0.0025/3)
dev.off()

binary.thres <- 1/10**15

F.binary <- fit$F > binary.thres
combn(colnames(F.binary), 2)

#number of topics that express gene  
pdf("tmp.pdf", height = 3, width = 6)
hist(rowSums(F.binary))
dev.off()

#
pdf("tmp1.pdf", height = 6, width = 6)
ggplot(fit$F %>% as.data.frame.matrix()) + geom_point(aes(k1, k1)) + scale_y_log10()+ scale_x_log10() 
ggplot(fit$F %>% as.data.frame.matrix()) + geom_point(aes(k6, k6)) + scale_y_log10()+ scale_x_log10() 
dev.off()

#heatmap to compare the expression (binary) between  all the topics
F.bin.n <- apply(F.binary[rowSums(F.binary) > 1 & rowSums(F.binary) < 9,], 2, as.numeric)
pdf("tmp2.pdf", height = 30, width = 6)
pheatmap::pheatmap(F.bin.n[rowSums(F.bin.n) > 1 & rowSums(F.bin.n) < 7,])
pheatmap::pheatmap(F.bin.n[rowSums(F.bin.n) < 9,])
dev.off()

#correlation between topic 4 and 6, and correlation between both to the rest
cor(fit$F[rownames(F.binary)[(rowSums(F.bin.n[,c(4, 6)]) == 1)], c(4, 6)], method = "spearman")
cor(fit$F[rownames(F.binary)[(rowSums(F.bin.n[,c(4, 6)]) == 1)], ], method = "spearman")
cor(fit$F[rownames(F.binary)[(rowSums(F.bin.n[,c(4, 8)]) == 1)], ], method = "spearman")

c4 <- rownames(fit$L)[(fit$L[,"k4"] > 0.6)]
c6 <- rownames(fit$L)[(fit$L[,"k6"] > 0.6)]
#add column wich topic cell is belong 
obj@meta.data$max.topic <- (setNames(colnames(fit$L)[apply(fit$L, 1, which.max)], rownames(fit$L)))[rownames(obj@meta.data)]

#color the cell on the umap by max topic 
pdf("tmp2.pdf", height = 6, width = 6)
DimPlot(obj,cells.highlight = c4)
DimPlot(obj,cells.highlight = c6)
DimPlot(obj, group.by = "max.topic")
dev.off()

#find the markers
obj@meta.data$max.topic
mks <-FindMarkers(obj, group.by = "max.topic", ident.1 = "k4", ident.2 = "k6")
mks.1 <-FindMarkers(obj, group.by = "max.topic", ident.1 = c("k4", "k6"))

##############################################################
obj@meta.data$max.topic_value <- (setNames(colnames(fit$F)[apply(fit$L, 1, which.max)], rownames(fit$L)))[rownames(obj@meta.data)]
#find the max of each gene expression
F.max <- fit$F %>% as.data.frame.matrix()
F.max$max_value <- apply(F.max, 1, max)
F.max$avg_value <- apply(F.max, 1, mean)
sum(F.max[,"max_value"])

#find the gene that make the 
fit <-  readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_8.rds")
obj@meta.data$max.topic <- (setNames(colnames(fit$L)[apply(fit$L, 1, which.max)], rownames(fit$L)))[rownames(obj@meta.data)]
mks_6_8 <-FindMarkers(obj, group.by = "max.topic", ident.1 = "k6")

fit <-  readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_9.rds")
obj@meta.data$max.topic <- (setNames(colnames(fit$L)[apply(fit$L, 1, which.max)], rownames(fit$L)))[rownames(obj@meta.data)]
mks_6_9 <-FindMarkers(obj, group.by = "max.topic", ident.1 = "k6")
mks_4_9 <-FindMarkers(obj, group.by = "max.topic", ident.1 = "k4")
mks_6And4_9 <-FindMarkers(obj, group.by = "max.topic", ident.1 = c("k4", "k6"))

common_rows <- Reduce(intersect, list( rownames(mks_6_8), rownames(mks_6And4_9)))

# 2
######################################################################################################
########## choose just the gene with the high score. ##############
######################################################################################################
fit <-  readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_9.rds")
F.max <- fit$F %>% as.data.frame.matrix()
F.max$max_topic <- apply(F.max, 1, which.max)
F.max$max_value <- apply(F.max[, -which(names(F.max) == "max_topic")], 1, max)
i <- 1
top_gene <- list()
for (i in 1:10) {
  top_gene <- c(top_gene, rownames(head(F.max %>% arrange(desc(max_value)) %>% filter(max_topic == i), n = 50)))
}

# 3
######################################################################################################
########## Select only the genes that appear to be differential expressed. ##############
###################################################################################################### 
top_gene <- list()

for (i in 1:12) {
  top_gene <- c(top_gene, rownames(F.max[rownames(de_12$lfsr[de_12$lfsr[,i]==0 & F.max$max_value > 0.002,]),]))
}
length(top_gene)
top_gene <- list()
top_gene <- c(top_gene, rownames(F.max[F.max$max_value > 1e-5,]))

F.max[rownames(de_12$lfsr[de_12$lfsr[,4]==0,]),]

# 4
######################################################################################################
### Select the marker for each topic by using a subset of cells with high expression of that topic.###
###################################################################################################### 
L.max <- fit$L %>% as.data.frame.matrix()
L.max$max_topic <- apply(L.max, 1, which.max)
L.max$min_topic <- apply(L.max, 1, which.min)
L.max$max_value <- apply(L.max[, !(names(L.max) %in% c("max_topic", "min_topic"))], 1, max)
L.max$min_value <- apply(L.max, 1, min)
L.max <- L.max[L.max$max_value>0.7,]
table(L.max$max_topic)

# find the signatures by choose markers between cells thea express at least 70% of topic agianst 
topics_markers <- vector("list", length = 9)
for (i in 1:9){
  topic_mark_cell <- rownames(L.max[L.max$max_topic==i,])
  topic_unmark_cell <- rownames(fit$L[fit$L[,paste0("k",i)] < 0.05, ])
  topics_markers[i] <- FindMarkers(obj, ident.1 = topic_mark_cell,ident.2 = topic_unmark_cell)
  topics_markers[i] <- (topics_markers[i][topics_markers[i][,"p_val"]<0.05,])
}
topic9_mark_cell <- rownames(L.max[L.max$max_topic==9,])
topic9_unmark_cell <- rownames(fit$L[fit$L[,"k9"] < 0.05, ])
markers9 <- FindMarkers(obj, ident.1 = topic9_mark_cell,ident.2 = topic9_unmark_cell)
markers9 <- (markers9[markers9[,"p_val"]<0.05,])
df <- F.max[rownames(markers9),]
df$diff_k9_max_value <- df$k9 - df$max_value
table(tmp$max_topic)
#now we explore the de of this gene 
i <- 9
de_9 <-  readRDS(paste0("/Volumes/habib-lab/shmuel.cohen/microglia/objects/de_psd_k_",i ,"_vsnull.rds"))
postmean <- de_9$postmean[rownames(de_9$postmean) %in% rownames(markers9),]
lfsr <- de_9$lfsr[rownames(de_9$lfsr) %in% rownames(markers9),]
z_score <- de_9$z[rownames(de_9$z) %in% rownames(markers9),]
boxplot(z_score, main= "z score topic 9 signatures vs the other topics")
boxplot(postmean, main= "postmean topic 9 signatures vs the other topics")
boxplot(lfsr, main= "lfsr topic 9 signatures vs the other topics")

length(rownames(de_9$postmean))
unique(rownames(de_9$postmean) %in% rownames(markers9))
ggplot(markers9%>% as.data.frame.matrix()) + geom_point(aes(p_val, avg_log2FC)) #+ scale_y_log10()+ scale_x_log10() 
colnames(markers9)

# 5
######################################################################################################
### "Create a subset of genes that appear to be more significant." ###
######################################################################################################
make_subset_of_gene <- function(fit_l, fit_r){
  l_counts <- fit_l$F
  r_counts <- fit_r$F
  row_sum <- rowSums(l_counts)
  nonzero_rows <- row_sum != 0
  r_counts <- r_counts[nonzero_rows,]
  row_sum <- rowSums(r_counts)
  nonzero_rows <- row_sum != 0
  r_counts <- r_counts[nonzero_rows,]
  #load table of DE anf choose the most significant
  des <- read.csv("/Users/shmuel/Downloads/de_vs_null_genes_enriched.csv") 
  #des <- des[des$lfsr > 0.001,]
  des <- des[order(des$lfsr, decreasing = TRUE), ] #the most p value significant 
  #des <- des[order(des$p.exp, decreasing = TRUE), ] #the most p value significant 
  # Get the top 1000 values from the sorted data frame
  des <- des[1:1000,]
  top_genes_lfc <- des %>% #the most lfc significant 
    group_by(topic) %>%
    arrange(desc(lfc)) %>%
    slice_head(n = 20) %>%
    ungroup() %>%
    pull(gene)
  top_genes_z <- des %>% #the most z score significant 
    group_by(topic) %>%
    arrange(desc(z)) %>%
    slice_head(n = 20) %>%
    ungroup() %>%
    pull(gene)
  #top_genes <- top_genes_lfc[top_genes_lfc %in% top_genes_z]  #choose the intersection 
  
  #another option:
  #fit <-  readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_9.rds")
  F.max <- l_counts %>% as.data.frame.matrix()
  F.max$max_topic <- apply(F.max, 1, which.max)
  F.max$max_value <- apply(F.max[, -which(names(F.max) == "max_topic")], 1, max)
  i <- dim(l_counts)[2]
  de <- readRDS(paste0("/Volumes/habib-lab/shmuel.cohen/microglia/objects/de_psd_k_",i ,"_vsnull.rds"))
  top_genes <- c( rownames(F.max[rownames(de$lfsr[de$lfsr[,i]==0 ,]),]))
  
  print(length(top_genes))
  ###end
  l_counts <- l_counts[rownames(l_counts) %in% top_genes,] #subset of counts
  r_counts <- r_counts[rownames(r_counts) %in% top_genes,] #subset of counts
  return(list(L = l_counts,R = r_counts))
}


