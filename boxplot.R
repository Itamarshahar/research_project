#box plot- edit of code from Roi, 
comparisons_specific <- list(
  'X15.k1' = list(c('HA', 'MCI'), c('HA', 'AD'), c('HA', 'SuperAgers'), c('HA', 'Young CTRL'), c('MCI', 'AD'), c('MCI', 'SuperAgers'), c('MCI', 'Young CTRL'), c('AD', 'SuperAgers'), c('AD', 'Young CTRL'), c('SuperAgers', 'Young CTRL')),
  'X15.k2' = list(c('HA', 'MCI'), c('HA', 'AD'), c('HA', 'SuperAgers'), c('HA', 'Young CTRL'), c('MCI', 'AD'), c('MCI', 'SuperAgers'), c('MCI', 'Young CTRL'), c('AD', 'SuperAgers'), c('AD', 'Young CTRL'), c('SuperAgers', 'Young CTRL')),
  'X15.k3' = list(c('HA', 'MCI'), c('HA', 'AD'), c('HA', 'SuperAgers'), c('HA', 'Young CTRL'), c('MCI', 'AD'), c('MCI', 'SuperAgers'), c('MCI', 'Young CTRL'), c('AD', 'SuperAgers'), c('AD', 'Young CTRL'), c('SuperAgers', 'Young CTRL')),
                   
  'X15.k15' = list(c('HA', 'MCI'), c('HA', 'AD')),
  'X15.k9' = list(c('HA', 'SuperAgers'))
  )


generate_comparisons_specific_values <- function(k, comparisons=NULL){
  if (is.null(comparisons)) {
    comparisons <- list(
      c('HA', 'MCI'), c('HA', 'AD'), c('HA', 'SuperAgers'), c('HA', 'Young CTRL'), 
      c('MCI', 'AD'), c('MCI', 'SuperAgers'), c('MCI', 'Young CTRL'), 
      c('AD', 'SuperAgers'), c('AD', 'Young CTRL'), c('SuperAgers', 'Young CTRL'))
    }
  
  comparisons_specific = list()
  for (k in paste0("X15.k", 1:15)) {
    comparisons_specific[[k]] <- comparisons
  }
  return (comparisons_specific)
}
# Extend pairwise comparisons for relevant 'X15.k' variables



add_topics_to_metadata <- function(obj, fits) {
  columns_list <- lapply(seq_along(fits), FUN = function(fit_index) {
    1:ncol(fits[[fit_index]]$L) %>% paste0(names(fits)[[fit_index]], ".k", .) %>% make.names()
  })
  names(columns_list) <- names(fits)
  for (fit_index in seq_along(fits)) {
    data <- data.frame(fits[[fit_index]]$L)
    cols <- columns_list[[fit_index]]
    obj <- AddMetaData(obj, metadata = data, col.name = cols)
  }
  
  list(obj = obj, columns_list = columns_list)
}


print_box_plot <- function(obj, fit, topic_dir, comparisons_specific, columns_list=NULL, size = 60) {
  k <- 1
  obj_with_list <- add_topics_to_metadata(obj, list(X15 = fit))
  obj <- obj_with_list[["obj"]]
  if (is.null(columns_list)) {
  columns_list <- obj_with_list[["columns_list"]]
  }
  print(columns_list)
  topic_k <- names(obj@meta.data)[names(obj@meta.data) %>% str_detect("X\\d+.k\\d+")] %>% str_extract(., "\\d+") %>% unique
  #print(columns_list)
  df <- obj@meta.data %>%
    dplyr::select(columns_list[[k]], Diagnosis_SampleID, Diagnosis) %>%
    dplyr::group_by(Diagnosis_SampleID, Diagnosis) %>%
    summarise(across(columns_list[[k]], mean))
  
  df$Diagnosis <- factor(df$Diagnosis, levels = c("AD", "MCI", "HA", "SuperAgers", "Young CTRL"))
  df$healthy <- df$Diagnosis %in% c("HA", "SuperAgers", "Young CTRL")
  
  withr::with_pdf(
    file.path(topic_dir, 1, paste0("MeanBoxPlot_", ncol(fit$L), ".pdf")),
    width = 10,
    height = as.integer(k) * size,
    {
      lapply(columns_list[[k]], function(col) {
        p <- ggboxplot(
          df,
          x = "Diagnosis",
          y = col,
          color = "Diagnosis",
          palette = "jco",
          add = "dotplot"
        )
        # Add p-value
        p + stat_compare_means(comparisons = comparisons_specific[[col]],
                               method = "wilcox.test")
      }) %>% plot_grid(ncol = 1, plotlist = .) %>% print()
    }
  )
}

generate_columns_list <- function(k, selected_topics = list(1)) {
  columns_list <- list()
  for (selected_topic in selected_topics) {
    if (length(columns_list) == 0) {
      columns_list[[1]] <- paste("X", k, ".k", selected_topic)
    } else {
      columns_list[[1]] <- c(columns_list[[1]], paste("X", k, ".k", selected_topic))
    }
  }
  columns_list[[1]] <- gsub(" ", "", columns_list[[1]])  # Remove spaces
  return(columns_list)
}


columns_list <- generate_columns_list(k=15, list(3, 5, 6, 7, 8, 9, 11, 13 ,14 ,15))
print_box_plot(obj = obj,
               fit = hippocampus_15,
               topic_dir = topic_dir,
               comparisons_specific = comparisons_specific, 
               columns_list=columns_list, 
               size = 30)

