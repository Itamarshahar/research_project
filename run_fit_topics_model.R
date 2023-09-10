args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  cat(" please provide 5 args: <path_to_the_input_data> <path_to_output_data> <K> <NC> <TRANSPOS(1/0)>")
} else {
  input_path <- args[1]
  output_path <- args[2]
  K <- as.integer(args[3])
  NC <- as.integer(args[4])
  TRANSPOSE <-as.integer(args[5])
}

if (!requireNamespace("fastTopics", quietly = TRUE)) {
  # Install the package
  install.packages("fastTopics", repos = cran_mirror)
}
library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
set.seed(1)
obj <- LoadH5Seurat(input_path)

################################################################################
## fetch the data: 
################################################################################

counts <- obj@assays$RNA@counts
if (TRANSPOSE == 1) {
  counts <- t(counts)
} 
################################################################################
## fit the topic model 
################################################################################
print(paste0("Starting the fit with K=", K))

# Remove the gene that have only zero (in all cells)
col_sums <- colSums(counts)
nonzero_cols <- col_sums != 0
counts <- counts[, nonzero_cols]

fit <- fit_topic_model(counts, 
                       k = K,
                       control.init = list(nc = NC),
                       control.main = list(nc = NC),
                       control.refine = list(nc = NC)
)
print("Finished the fit")

saveRDS(fit,  paste0(output_path, "fitted_topic_model_k_", K, ".rds"))
