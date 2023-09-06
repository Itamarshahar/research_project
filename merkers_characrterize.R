################################################################################
## This script suppose to validate significant genes of topics
################################################################################

################################################################################
## # Load the object
################################################################################
path_to_de <- "/Users/shmuel/Downloads/de_psd_0.1"
path_to_markers <- "/Users/shmuel/Downloads/cell_marker_24_8_23.csv"
de <-readRDS(path_to_de)
# open the markers file
markers <- read.csv(path_to_markers, header = TRUE, row.names = 1)


################################################################################
## # Filter the gene
################################################################################
K <- 10
significant_genes <- list()
genes <- rownames(de$lfsr)
for (i in 1:K) {
 markers <- intersect(genes[de$lfsr[, i] < 0.001], genes[de$est[, i] > 0])
 significant_genes <- append(significant_genes, list(markers))
}


markers$representative_topic <- c(8, 4,3,10,1,6,7,2)

#marge Nneurons type
markers$Markers[5] <- paste(markers$Markers[5], markers$Markers[6], markers$Markers[7], sep = ",")
markers$representative_topic[5] <- paste(markers$representative_topic[5], markers$representative_topic[6], markers$representative_topic[7], sep = ",")
# Remove the original second and third rows from the data frame
markers <- markers[-c( 6, 7), ]

#rename the row
row_index <- which(rownames(markers) == "Mature_neurons")
# Change the row name to "Nneurons"
rownames(markers)[row_index] <- "Nneurons"

markers$sam_gene <- c(2, 5,5,9,12,5)


markers$gene_in_significant_genes <- c(8, 4,3,10,1,6)
for (i in seq_len(nrow(markers))) {
  if (i==5){
    markers$representative_topic[i]
    list1 <- c(significant_genes[[1]], significant_genes[[6]], significant_genes[[7]])
  }
  else{
    list1 <- significant_genes[[as.numeric(markers$representative_topic[i])]]
  }
  list2 <- markers$Markers[i]
  list2 <- unlist(strsplit(list2, ","))
  intersection <- intersect(list1, list2)
  markers$gene_in_significant_genes[i] <- length(intersection)
}


markers$gene_in_significant_genes <- c(8, 4,3,10,1,6)
for (i in seq_len(nrow(markers))) {
  for (j in 1:K){
    list1 <- significant_genes[[as.numeric(markers$representative_topic[i])]]
    list2 <- markers$Markers[i]
    ist2 <- unlist(strsplit(list2, ","))
    markers$gene_in_significant_genes[i] <- length(intersection)
  }
  if (i==5){
    markers$representative_topic[i]
    list1 <- c(significant_genes[[1]], significant_genes[[6]], significant_genes[[7]])
  }
  else{
    list1 <- significant_genes[[as.numeric(markers$representative_topic[i])]]
  }
  list2 <- markers$Markers[i]
  list2 <- unlist(strsplit(list2, ","))
  intersection <- intersect(list1, list2)
  markers$gene_in_significant_genes[i] <- length(intersection)
}

