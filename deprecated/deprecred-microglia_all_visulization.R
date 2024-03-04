#visualization flow

run_all_visualization <- function(path_to_obj, path_to_plots, path_to_cortex="/Volumes/habib-lab-1/shmuel.cohen/microglia/objects/cortex_all_microglia_topic_fit.15.RDS" ){
  obj <- readRDS(path_to_obj)
  fit_files_paths <- c(
    "/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_8.rds",
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_9.rds",
    "/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_10.rds",
    "/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_11.rds",
    "/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_12.rds",
    "/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_13.rds",
    "/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_14.rds",
    "/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds"
    
  )
  
  #run visualization topic model
  source("~/Desktop/project/research_project/microglia_visualization_topic_model.R")
  for (fit_file in fit_files_paths) {
    fit <- readRDS(fit_file)
    k <- as.integer(dim(fit$L)[2])
    path_to_plots_for_k <- glue(path_to_plots, "Plots_for_k={k}/")
    run_main_flow(obj, fit, path_to_plots_for_k) #from ~/Desktop/project/research_project/microglia_visualization_topic_model.R"
  }
  
  #run correlation between the topics
  source("~/Desktop/project/research_project/microglia_topic_correlation.R")
  run_topic_evaluation(fit_files_paths, glue(path_to_plots, "correlation/")) #from "~/Desktop/project/research_project/microglia_topic_correlation.R"
  
  #check correlation vs cortex topics
  source("~/Desktop/project/research_project/correlation_by_genes.R")
  cortex_fit_15 <- readRDS(path_to_cortex)
  correlation_with_cortex(obj = obj,
                          cortex_fit_15 = cortex_fit_15,
                          fits_list = fit_files_paths,
                          path_to_plots = glue(path_to_plots, "correlation/cortex/"),
                          path_to_predicted = "/Volumes/habib-lab/shmuel.cohen/microglia/objects/predict_1000t.RDS")
    
}