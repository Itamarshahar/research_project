require_packages <- c("Matrix", "dplyr", "Seurat", "ggpubr")
lapply(require_packages, require, character.only = TRUE)  # Load the required packages

# Source the external file with microglia markers
source("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/microglia_gene_makers.R", local = TRUE)

# Load the objects
obj <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/filtered_microglia_with_topic.rds")
fit <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds")

# Function to get relevant markers
get_relevant_microglia_markers <- function(with_category = FALSE, common_genes = NULL) {
  markers <- get_microglia_markers()[["microglia"]]
  relevant_markers <- c()

  for (category in names(markers)) {
    if (!is.null(markers[[category]]$genes)) {
      filtered_genes <- if (!is.null(common_genes)) {
        markers[[category]]$genes[markers[[category]]$genes %in% common_genes]
      } else {
        markers[[category]]$genes
      }

      if (with_category) {
        relevant_markers <- c(relevant_markers, paste(category, filtered_genes, sep = "_"))
      } else {
        relevant_markers <- c(relevant_markers, filtered_genes)
      }
    }
  }
  return(relevant_markers)
}

check_and_create_directory <- function(filepath) {
  dir <- dirname(filepath)  # Get the directory part of the file path
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)  # Create directory if it doesn't exist
    message("Directory created: ", dir)
  } else {
    message("Directory already exists: ", dir)
  }
}
# Function to check if file exists
check_file_exists <- function(filepath) {
  if (file.exists(filepath)) {
    message("File exists: ", filepath)
    return(TRUE)
  } else {
    message("File does not exist: ", filepath)
    return(FALSE)
  }
}

# Wrapper function for heatmap generation
wrapper_heatmap_topic_gene_expression <- function(obj, fit, filename="/Volumes/habib-lab/shmuel.cohen/microglia/plots/heatmap_topic_gene_expression.pdf", override_file=FALSE) {
  count_mat <- t(obj[['RNA']]@counts)
  l_mat <- t(fit$L)

  genes_of_interest <- as.vector(get_relevant_microglia_markers())
  genes_of_interest <- genes_of_interest[genes_of_interest %in% colnames(count_mat)]

  if (check_file_exists(filename) && !override_file) {
    message("File already exists. wrapper_heatmap_topic_gene_expression didn't run")
    return(NULL)
  }
  check_and_create_directory(filename)
  # Preprocessing count matrix
  preprocessed_count_mat <- count_mat[, genes_of_interest]

  # Calculate the genes_topics matrix
  genes_topics_mat <- l_mat %*% preprocessed_count_mat
  genes_topics_mat <- sweep(genes_topics_mat, 1, colSums(t(l_mat)), "/")
  colnames(genes_topics_mat) <- genes_of_interest

  # Plot heatmap
  title <- "Scaled Average Gene Expression Per Topic"
  pheatmap::pheatmap(t(genes_topics_mat[,]),
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     display_numbers = FALSE,
                     number_format = "%.1f",
                     scale = "row",
                     show_rownames = TRUE,
                     main = title,
                     filename = filename,
                     fontsize_row = 3
  )
}

# Example usage

#wrapper_heatmap_topic_gene_expression(obj, fit, filename="/Volumes/habib-lab/shmuel.cohen/microglia/plots/heatmap_topic_gene_expression2.pdf", override_file = TRUE)
