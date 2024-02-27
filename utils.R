#' Title csv_to_markers
#' the function takes a path to csv and returns a markers as a DataFrame 
#' @param path 
#'
#' @return the markers as a DataFrame, list("Oligo" : list("marker1", "marker2"...))

csv_to_markers <- function(path) {
  cell_markers_read <- read.csv(path)
  # Convert the Markers column to a list of vectors
  cell_markers <- lapply(strsplit(cell_markers_read$Markers, ","), trimws)
  # Create a named list
  cell_markers <- setNames(cell_markers, cell_markers_read$Cell_Type)
  return(cell_markers)
}


markers_to_csv <- function(markers, path = "cell_markers.csv") { # markers is list(key:list(char))

  # Convert the cell_markers list to a data frame with one row per cell type
  cell_markers_df <- data.frame(Cell_Type = names(markers),
                                Markers = sapply(markers, paste, collapse = ","))

  write.csv(cell_markers_df, path, row.names = FALSE)
  return(path)
}


################################################################################
## Visualizations   
################################################################################


dot_plot <- function(data,
                     columns,
                     group.vector,
                     order = TRUE,
                     group.levels = NULL,
                     do.return.order = F,
                     alpha.threshold = 0.1,
                     mean_threshold = 0
) {
  alpha <- function(x) { mean(x > alpha.threshold) }
  group.mean <- function(x) { mean(x[x > mean_threshold]) }

  columns <- columns[columns %in% colnames(data)]
  t1 <- reshape2::melt(data.frame(data.frame(data[, columns], groups = as.factor(group.vector)) %>%
                                    group_by(groups) %>%
                                    summarise_all(funs(alpha))))
  t2 <- reshape2::melt(data.frame(data.frame(data[, columns], groups = as.factor(group.vector)) %>%
                                    group_by(groups) %>%
                                    summarise_all(funs(group.mean))))

  data.summarized <- merge(t1,
                           t2,
                           by = c('groups', 'variable'),
                           suffixes = c('.alpha', '.group.mean'))
  data.summarized$variable <- make.names(data.summarized$variable)
  # genes.of.interest <- make.names(genes.of.interest)
  if (order) {
    hclustering <- hclust(as.dist(1 - cor(data[, columns])), method = 'ward.D2')
    data.summarized$variable <- factor(data.summarized$variable, levels = make.names(columns[hclustering$order]))
  } else {
    data.summarized$variable <- factor(data.summarized$variable, levels = make.names(columns))
  }
  data.summarized <- data.summarized[order(data.summarized$variable),]
  if (!is.null(group.levels)) {
    levels(data.summarized$groups) <- group.levels
  }
  P <- ggplot(data.summarized, aes(groups, variable)) +
    geom_point(aes(size = value.alpha, color = value.group.mean), stroke = 0) +
    scale_color_viridis(direction = -1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(size = paste0('% cells have more than ', alpha.threshold),
         color = paste0('mean x > ', mean_threshold)) +
    coord_flip()
  if (do.return.order) { return(columns[hclustering$order]) }
  return(P)
}

dotplot_topics <- function(obj,
                           topic_columns,
                           group.by,
                           ...) {
  data <- FetchData(obj, c(topic_columns, group.by))
  dot_plot(data, topic_columns, group.vector = data[[group.by]], ...)
}


#' Combine Topics and Metadata in a Seurat Object
#'
#' This function combines topics generated from a topic modeling fit and adds them
#' as columns to the metadata of a Seurat object.
#'
#' @param obj A Seurat object.
#' @param fit A topic modeling fit object.
#' @param path_to_new_obj Path where the new Seurat object will be saved (optional).
#'
#' @return A modified Seurat object with topics added to metadata.
#'
#' @examples
#' # Load your Seurat object and topic modeling fit
#' seurat_obj <- LoadSeuratObject("path/to/seurat_obj.rds")
#' lda_fit <- LDA(topic_matrix, ...)
#'
#' # Combine topics and metadata
#' seurat_obj <- combine_topics_and_meta_data(seurat_obj, lda_fit)
#'
#' # Optionally, save the modified Seurat object
#' # combine_topics_and_meta_data(seurat_obj, lda_fit, "path/to/save/")
#'
#' @export
combine_topics_and_meta_data <- function(obj, fit, path_to_new_obj = NULL) {
  obj@meta.data <- cbind(obj@meta.data, fit$L)
  if (!is.null(path_to_new_obj)) {
    SaveH5Seurat(obj,
                 filename = paste0(path_to_new_obj, "obj_and_topics_", K),
                 overwrite = TRUE)
  }
  return(obj)
}


claculate_correlation <- function(mat1, mat2, reorder = TRUE, method = 'pearson', ignore_missing_rows = TRUE) {
  # We we want to get a distance matrix in order to run hclust on it.
  # So we compare the correlation between all the topics even from the same fit

  colnames(mat1) <- paste("1", colnames(mat1), sep = "_")
  colnames(mat2) <- paste("2", colnames(mat2), sep = "_")

  if (ignore_missing_rows && nrow(mat1) != nrow(mat2)) {
    common_rownames <- intersect(rownames(mat1), rownames(mat2))
    message("Correlation between ", length(common_rownames), " rows")
    mat1 <- mat1[common_rownames,]
    mat2 <- mat2[common_rownames,]
  }

  mat <- cbind(mat1, mat2)
  correlation <- cor(mat, mat, method = method)

  if (isTRUE(reorder) && !matrixcalc::is.square.matrix(correlation)) {
    warning("Reorder isn't supported on non square correlation matrix")
    reorder <- FALSE
  }

  if (isTRUE(reorder)) {
    hclustering <- hclust(as.dist((1 - correlation) / 2), method = 'ward.D2')
    correlation <- correlation[hclustering$order, hclustering$order]
  }

  # remove unwanted correlations between topics from the same fit
  correlation <- correlation[!colnames(correlation) %in% colnames(mat1),
                             !colnames(correlation) %in% colnames(mat2)]

  colnames(correlation) <- gsub("^1_", "", colnames(correlation))
  rownames(correlation) <- gsub("^2_", "", rownames(correlation))
  return(correlation)
}


HeatmapHelper_add_values_to_display <- function(correlation, OnlyPositive = FALSE) {
  if (OnlyPositive == TRUE) {

    cell_fun = function(j, i, x, y, width, height, fill) {
    {
      if (correlation[i, j] > 0)
        grid.text(sprintf("%.2f", correlation[i, j]), x, y, gp = gpar(fontsize = 10, fontfamily = 'Helvetica')) }
    }

  }
  else {

    cell_fun = function(j, i, x, y, width, height, fill) {
    {
      grid.text(sprintf("%.1f", correlation[i, j]), x, y, gp = gpar(fontsize = 10, fontfamily = 'Helvetica')) }
    }

  }

  return(cell_fun) }


generate_ScatterPlot_DE <- function(de_enriched,
                                    path_to_plots,
                                    K,
                                    p.exp = F,
                                    file_name = "DeStatsScatter.pdf") {
  if (p.exp) {
    measures <- c('gene_score', 'lfc', 'lfsr', 'p.exp') # gene_score the F value, lfc is the postmean, lsfr (is the same), p.exp is the number of cells that expression the gene
  }
  else { measures <- c('gene_score', 'lfc', 'lfsr')
  }
  # gene_score the F value, lfc is the postmean, lsfr (is the same), p.exp is the number of cells that expression the gene
  lapply(X = combn(measures, m = 2, simplify = FALSE),
         FUN = function(column, K = 10) {
           p <- ggplot(de_enriched %>% mutate(lfsr = lfsr + 1e-30), aes_string(x = column[[1]], y = column[[2]])) +
             geom_pointdensity() +
             scale_color_viridis() +
             facet_wrap(~topic)
           if (column[[1]] %in% c("lfsr", 'gene_score', 'p.exp'))
             p <- p + scale_x_log10()
           if (column[[2]] %in% c("lfsr", 'gene_score', 'p.exp'))
             p <- p + scale_y_log10()
           p
         }
  ) %>%
    plot_grid(plotlist = ., ncol = 1) %>%
    ggsave(filename = file.path(path_to_plots, file_name),
           height = (as.integer(K) / 2) * 6,
           width = 10,
           plot = ., limitsize = FALSE)
}


# Define the function
save_to_google_drive <- function(file_path, folder_link) {
  library(googledrive)

  # Authenticate and set up your Google Drive
  drive_auth()

  # Extract the folder ID from the provided link
  extracted_id <- sub(".*/folders/(.*)", "\\1", folder_link)

  # Set the extracted folder ID as the target folder
  target_folder <- drive_get(as_id(extracted_id))

  # Upload the PDF file to the target folder
  drive_upload(file_path, name = basename(file_path), path = target_folder$id)

}
