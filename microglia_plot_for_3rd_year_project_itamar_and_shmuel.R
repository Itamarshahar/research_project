path_to_obj <- "/Volumes/habib-lab/shmuel.cohen/microglia/objects/filtered_microglia_with_topic.rds"
path_to_fit <- "/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds"

# path of the working directory
working_directory_path <- "/Volumes/habib-lab/shmuel.cohen/microglia/plots/plot_for_3rd_year_project"

#path for figures to be saved
# figure_heatmap_path <- file.path(working_directory_path, "heatmap_topic_gene_expression.pdf")
figure_heatmap_path <- file.path(working_directory_path, "heatmap_topic_gene_expression_annotated.pdf")
figure_violin_path <- file.path(working_directory_path, "VlnPlot_Diagnosis_SampleID.pdf")


# check_and_create_file <- function(filepath) {
#   if (!file.exists(filepath)) {
#     file.create(filepath)
#     cat("File created:", filepath, "\n")
#   } else {
#     cat("File already exists, the file will be override:", filepath, "\n")
#   }
# }



wrapper_violin_plot <- function(obj, filename, K=15) {
  check_and_create_file(filename)
  de_VlnPlot <- list()
  for (i in 1:K) {
    de_VlnPlot[[i]] <- VlnPlot(obj,
                               features = paste0("k", i),
                               pt.size = 0,
                               group.by = "SampleID_Diagnosis",
                               split.by = "Diagnosis",
    ) +
      NoLegend() +
      labs(x = NULL) +
      theme(axis.text.x = element_text(size = 8))
  }

  combined_plot <- arrangeGrob(grobs = de_VlnPlot,
                               ncol = 1,
                               nrow = K,
  )
  ggsave(combined_plot,
         file = filename,
         width = K,
         height = K * 2.4)
}

microglia_plot_for_3rd_year_project_itamar_and_shmuel <- function() {
   # Load the objects
  obj <- readRDS(path_to_obj)
  fit <- readRDS(path_to_fit)
  wrapper_heatmap_topic_gene_expression(obj, fit, filename = figure_heatmap_path, override_file = TRUE)
  wrapper_violin_plot(obj, filename = figure_violin_path)



}

# Heatmap
