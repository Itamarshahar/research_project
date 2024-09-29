# require_packages <- c("Matrix", "dplyr", "Seurat", "ggpubr")
# lapply(require_packages, require, character.only = TRUE)  # Load the required packages

# Source the external file with microglia markers
#
# # Load the objects
# obj <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/filtered_microglia_with_topic.rds")
# fit <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds")

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

flatten_atlas_simple <- function(atlas, genes_of_interest) {
  result <- data.frame(Gene = character(), Category = character(), stringsAsFactors = FALSE)
  for (category in names(atlas)) {
    for (subcategory in names(atlas[[category]])) {
      if (subcategory == "genes") {
        for (gene in atlas[[category]][[subcategory]]) {
          if (gene %in% genes_of_interest) {
            result <- rbind(result, data.frame(Gene = gene, Category = category, stringsAsFactors = FALSE))
          }
        }
      }
    }
  }
  return(result)
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

check_file_exists <- function(filepath) {
  if (file.exists(filepath)) {
    message("File exists: ", filepath)
    return(TRUE)
  } else {
    message("File does not exist: ", filepath)
    return(FALSE)
  }
}

generate_distinct_colors <- function(annotation, seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)

  # Get unique values from the annotation
  unique_values <- unique(annotation)
  n <- length(unique_values)

  # Generate distinct colors
  hues <- seq(15, 375, length = n + 1)
  colors <- hcl(h = hues, l = 65, c = 100)[1:n]

  # Create a named vector of colors
  color_vector <- setNames(colors, unique_values)

  return(color_vector)
}
# Wrapper function for heatmap generation
wrapper_heatmap_topic_gene_expression <- function(obj, fit, filename="/Volumes/habib-lab/shmuel.cohen/microglia/plots/heatmap_topic_gene_expression.pdf", override_file=FALSE) {
  require_packages <- c("Matrix", "dplyr", "Seurat", "ggpubr")
  lapply(require_packages, require, character.only = TRUE)
  source("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/microglia_gene_makers.R", local = TRUE)

  count_mat <- t(obj[['RNA']]@counts)
  l_mat <- t(fit$L)

  genes_of_interest <- as.vector(get_relevant_microglia_markers())
  genes_of_interest <- genes_of_interest[genes_of_interest %in% colnames(count_mat)]

  check_and_create_file(filename)
  # Preprocessing count matrix
  preprocessed_count_mat <- count_mat[, genes_of_interest]

  # Calculate the genes_topics matrix
  genes_topics_mat <- l_mat %*% preprocessed_count_mat
  genes_topics_mat <- sweep(genes_topics_mat, 1, colSums(t(l_mat)), "/")

  title <- "Scaled Average Gene Expression Per Topic"
  # atlas <- get_microglia_markers()
  # heatmap_annotation_simple <- flatten_atlas_simple(atlas$microglia, genes_of_interest)
  # message(dim(heatmap_annotation_simple))
  # message(dim(t(genes_topics_mat[,])))

  # Create a data frame with the correct structure for pheatmap
  # annotation_df <- data.frame(
  # Category = heatmap_annotation_simple$Category,
  # row.names = heatmap_annotation_simple$Gene
# )

# Ensure annotation_df only includes rows that are in the heatmap matrix
# annotation_df <- annotation_df[rownames(genes_topics_mat), , drop = FALSE]
annotation_colors <- generate_distinct_colors(annotation_df$Category)
# Now use this in your pheatmap function
# pheatmap::pheatmap(t(genes_topics_mat[,]),
#                    # cluster_rows = TRUE,
#                    cluster_rows = FALSE,
#                    # cluster_cols = FALSE,
#                    cluster_cols = TRUE,
#                    display_numbers = FALSE,
#                    number_format = "%.1f",
#                    scale = "row",
#                    show_rownames = TRUE,
#                    main = title,
#                    filename = filename,
#                    fontsize_row = 10,
#                    annotation_row = annotation_df,
#                    annotation_colors = list(Category = annotation_colors),
#                    # annotation_names_row = TRUE,
#                    angle_col = 45,
# )
  mat <- t(genes_topics_mat[,])
  mat_scaled <- t(scale(t(mat)))  # Scale by row

# Generate colors for all categories in the annotation
all_categories <- lapply(annotation_df, function(col) generate_distinct_colors(col))

# Create the annotation
ha_row <- HeatmapAnnotation(df = annotation_df,
                            col = all_categories,
                            which = "row",
                            show_annotation_name = TRUE)

# Create the main heatmap
ht <- Heatmap(mat_scaled,
              name = "Expression",
              cluster_rows = FALSE,
              cluster_columns = TRUE,
              show_column_dend = FALSE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 10),
              column_names_rot = 45,
              right_annotation = ha_row,
              column_title = title,
              )
  pdf(filename)
  draw(ht)
  dev.off()

}

# Example usage

wrapper_heatmap_topic_gene_expression(obj, fit, filename="/Volumes/habib-lab/shmuel.cohen/microglia/plots/plot_for_3rd_year_project/VlnPlot_Diagnosis_SampleID.pdf", override_file = TRUE)
