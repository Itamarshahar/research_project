#box plot- edit of code from Roi, 
topic_k <- names(obj@meta.data)[names(obj@meta.data) %>% str_detect("k\\d+")] %>% str_extract(., "\\d+") %>% unique

columns_list <- sapply(topic_k, FUN = function(k) {
  1:as.integer(k) %>% paste0(k, ".k", .) %>% make.names()
}, simplify = FALSE)

df <- obj@meta.data %>%
  dplyr::select(columns_list[[k]], Diagnosis_SampleID, Diagnosis) %>%
  dplyr::group_by(Diagnosis_SampleID, Diagnosis) %>%
  summarise(across(columns_list[[k]], mean))

df$Diagnosis <- factor(df$Diagnosis, levels = c("AD", "MCI", "HA", "SuperAgers", "Young CTRL"))
df$healthy <- df$Diagnosis %in% c("HA", "SuperAgers", "Young CTRL")

topic_dir <- "/Volumes/habib-lab/shmuel.cohen/microglia/plots/Plots_for_k=15/boxplot"

withr::with_pdf(file.path(topic_dir, 1, "MeanBoxPlot.pdf"), width = 10, height = as.integer(k) * 4*6, {
    lapply(columns_list[[k]], function (col) {
    p <- ggboxplot(df, x = "Diagnosis", y = col,
                   color = "Diagnosis", palette = "jco",
                   add = "dotplot")
    #  Add p-value
    p + stat_compare_means(comparisons =comparisons_specific[[col]] ,
                           #ref.group = ".all."
    ) + stat_summary(fun.data = function(x) data.frame(y=max(x) * 1.17, label = paste("Mean=",round(mean(x), digits = 2))), geom="text")
  }) %>% plot_grid(ncol = 2, plotlist = .) %>% print()
})

comparisons_specific <- list(
  'X18.k8' = list(c('HA', 'SuperAgers'), c('HA', 'Young CTRL')),
  'X18.k2' = list(c('HA', 'SuperAgers'), c('SuperAgers', 'Young CTRL'))
)

add_topics_to_metadata <- function(obj, fits) {
  columns_list <- lapply(seq_along(fits), FUN = function(fit_index) {
    1:ncol(fits[[fit_index]]) %>% paste0(names(fits)[[fit_index]], ".k", .) %>% make.names()
  })
  names(columns_list) <- names(fits)
  for (fit_index in seq_along(fits)) {
    data <- data.frame(fits[[fit_index]]$L)
    cols <- columns_list[[fit_index]]
    obj <- AddMetaData(obj, metadata = data, col.name = cols)
  }
  
  list(obj = obj, columns_list = columns_list)
}

add_topics_to_metadata(obj, topic_models)



variable_names <- c("nCount_RNA", "nFeature_RNA", "orig.ident", "sample_number", 
                    "SampleID", "Diagnosis", "Age", "pmi", "Gender", "BRAAK", 
                    "Overall.CAA.Score", "C.score", "Thal", "Lewy.Body.Disease", 
                    "APOE.Genotype", "apoe_2", "apoe_3", "apoe_4", "percent.mt", 
                    "nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.2", "SCT_snn_res.0.4", 
                    "SCT_snn_res.0.8", "SCT_snn_res.1.2", "SCT_snn_res.1.6", 
                    "seurat_clusters", "RNA_snn_res.0.8", "RNA_snn_res.0.2", 
                    "SampleID_Diagnosis", "Diagnosis_SampleID", 
                    paste0("X15.k", 1:15))
variable_names_list <- as.list(variable_names)


colnames(obj@meta.data) <- variable_names
