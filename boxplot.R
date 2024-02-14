#box plot- edit of code from Roi, 

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


print_box_plot <- function(obj, topic_dir, comparisons_specific){
  k <- 1
  obj_with_list <- add_topics_to_metadata(obj, list(X15=fit))
  obj <- obj_with_list[["obj"]]
  columns_list <- obj_with_list[["columns_list"]]
  topic_k <- names(obj@meta.data)[names(obj@meta.data) %>% str_detect("X\\d+.k\\d+")] %>% str_extract(., "\\d+") %>% unique
  
  df <- obj@meta.data %>%
    dplyr::select(columns_list[[k]], Diagnosis_SampleID, Diagnosis) %>%
    dplyr::group_by(Diagnosis_SampleID, Diagnosis) %>%
    summarise(across(columns_list[[k]], mean))
  
  df$Diagnosis <- factor(df$Diagnosis, levels = c("AD", "MCI", "HA", "SuperAgers", "Young CTRL"))
  df$healthy <- df$Diagnosis %in% c("HA", "SuperAgers", "Young CTRL")
  
  withr::with_pdf(file.path(topic_dir, 1, "MeanBoxPlot.pdf"), width = 10, height = as.integer(k) * 4*15, {
    lapply(columns_list[[k]], function (col) {
      p <- ggboxplot(df, x = "Diagnosis", y = col,
                     color = "Diagnosis", palette = "jco",
                     add = "dotplot")
      #  Add p-value
      p + stat_compare_means(comparisons =comparisons_specific[[col]] ,
                             #ref.group = ".all."
      ) + stat_summary(fun.data = function(x) data.frame(y=max(x) * 1.17, label = paste("Mean=",round(mean(x), digits = 2))), geom="text")
    }) %>% plot_grid(ncol = 1, plotlist = .) %>% print()
  })
  
}


print_box_plot(obj, topic_dir,comparisons_specific )
