#' Title csv_to_markers
#' the function takes a path to csv and returns a markers as a DataFrame 
#' @param path 
#'
#' @return the markers as a DataFrame, list("Oligo" : list("marker1", "marker2"...))

csv_to_markers <- function(path){
  cell_markers_read <- read.csv(path)
  # Convert the Markers column to a list of vectors
  cell_markers <- lapply(strsplit(cell_markers_read$Markers, ","), trimws)
  # Create a named list
  cell_markers <- setNames(cell_markers, cell_markers_read$Cell_Type)
  return (cell_markers)
}


markers_to_csv <- function(markers, path="cell_markers_tmp.csv"){ # markers is list(key:list(char))
  
  # Convert the cell_markers list to a data frame with one row per cell type
  cell_markers_df <- data.frame(Cell_Type = names(markers),
                                Markers = sapply(markers, paste, collapse = ","))
  
  write.csv(cell_markers_df, path, row.names = FALSE)
  return (path)
}





################################################################################
## Visualizations   
################################################################################


dot_plot <- function(data, 
                     columns,
                     group.vector,
                     order=TRUE,
                     group.levels=NULL,
                     do.return.order=F,
                     alpha.threshold = 0.1,
                     mean_threshold = 0
){
  alpha <- function(x) {mean(x > alpha.threshold)}
  group.mean <- function(x) {mean(x[x > mean_threshold])}
  
  columns <- columns[columns %in% colnames(data)]
  t1 <- reshape2::melt(data.frame(data.frame(data[,columns],groups = as.factor(group.vector)) %>%
                                    group_by(groups) %>%
                                    summarise_all(funs(alpha))))
  t2 <- reshape2::melt(data.frame(data.frame(data[,columns],groups = as.factor(group.vector)) %>%
                                    group_by(groups) %>%
                                    summarise_all(funs(group.mean))))
  
  data.summarized <- merge(t1,
                           t2,
                           by=c('groups','variable'),
                           suffixes=c('.alpha','.group.mean'))
  data.summarized$variable <- make.names(data.summarized$variable)
  # genes.of.interest <- make.names(genes.of.interest)
  if (order){
    hclustering <- hclust(as.dist(1-cor(data[,columns])),method='ward.D2')
    data.summarized$variable <- factor(data.summarized$variable,levels=make.names(columns[hclustering$order]))
  } else {
    data.summarized$variable <- factor(data.summarized$variable,levels=make.names(columns))
  }
  data.summarized <- data.summarized[order(data.summarized$variable),]
  if (!is.null(group.levels)){
    levels(data.summarized$groups) <- group.levels
  }
  P <- ggplot(data.summarized,aes(groups,variable)) +
    geom_point(aes(size=value.alpha,color=value.group.mean),stroke=0) +
    scale_color_viridis(direction=-1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(size=paste0('% cells have more than ', alpha.threshold),
         color=paste0('mean x > ', mean_threshold)) +
    coord_flip()
  if (do.return.order){return(columns[hclustering$order])}
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
