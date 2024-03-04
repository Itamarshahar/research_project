

load_libraries <- function() {
  library(Seurat)
  library(withr)
  library(SeuratDisk)
  #library(SeuratData)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(plotly)
  library(hrbrthemes)
  source("utils.R")
  library(Matrix)
  library(fastTopics)
  library(cowplot)
  library(MASS)
}

run_preprocess <- function(obj, path_to_plots) {
  # add the SampleID_Diagnosis column to the meta data
  obj@meta.data["SampleID_Diagnosis"] <- paste(obj$SampleID, obj$Diagnosis, sep = "_")
  obj [["SampleID_Diagnosis"]] <- paste(obj$SampleID, obj$Diagnosis, sep = "_")
  # adding the MT precantage to the meta data
  obj [["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-", assay = "RNA")

  return(obj)
}

run_transformation <- function(obj) {
  obj <- SCTransform(obj,
                     #conserve.memory = TRUE,
                     vars.to.regress = "percent.mt",
  )
  return(obj)
}
run_generate_box_plot <- function (obj, fit){
  source("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/boxplot.R")
  generate_box_plot(obj, fit)
}

run_de_flow <- function(hippocampus_exist=FALSE){
  generate_de <- function (obj, fit, assay="RNA", organism="hsa", z_score_log_val =1, percent_cells_val = 0.01){
  reweighted_f<-topic_reweight_f(fit$F)
  df<-df_pos_des(obj, reweighted_f, fit, assay=assay, organism=organism)
  des_df<-df %>% dplyr::filter(z_score_log > z_score_log_val & percent_cells > percent_cells_val)
  return (des_df)
}
  if (!hippocampus_exist) {
   de_hippocampus_X15 <- generate_de(obj = readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/filtered_microglia.rds"),
                                    fit = readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds"))
    saveRDS(de_hippocampus_X15, "/Volumes/habib-lab/shmuel.cohen/microglia/objects/DE_hippocampus_X15_.rds")
  }
  else {
          print("de object exists")
  }

}

run_correlation_with_cortex_flow <- function() {
  path <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/microglia_topic_correlation_with_cortex.R"
  source(path)
  run_corralation_with_cortex(fit_hippo="/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds",
                                        fit_cortex="/Volumes/habib-lab/shmuel.cohen/microglia/objects/cortex_all_microglia_topic_fit.15.RDS",
                                        de_hippocampus="/Volumes/habib-lab/shmuel.cohen/microglia/objects/DE_hippocampus_X15.rds",
                                        de_cortex="/Volumes/habib-lab/shmuel.cohen/microglia/objects/DE_cortex_X15.csv")
}

main <- function(path_to_obj, run_filler= "NA", run_subset ="NA") {
  path_to_obj <- "/Users/shmuel/SuperAgerRemoveSample7264-2Microglia.h5seurat"
  path_to_plots <- "/Volumes/habib-lab/shmuel.cohen/microglia/plots/"
  path_to_objs <- "/Volumes/habib-lab/shmuel.cohen/microglia/objects/"
  
  load_libraries()
  obj <- LoadH5Seurat(path_to_obj)
  obj <-set_sample_and_diagnosis_order(obj)
  
  #filder microglia
  #if run_subset{}
  source("~/Desktop/project/research_project/microglia_filter_redundence.R")
  obj_subset <- get_filtered_obj(obj, path_to_plots, path_to_objs)
  obj_subset <- run_preprocess(obj_subset) # todo maybe run before make subset 
  # obj_subset <- run_transformation(obj_subset)
  
  #### here we need to run the fit_topic_model but it takes time so we skip it ####

  #list of the links to the fit files
  fit_files_paths <- c(
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_8.rds",
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_9.rds",
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_10.rds",
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_11.rds",
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_12.rds",
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_13.rds",
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_14.rds",
    "/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds"
    
  )
  
  #run visualization topic model
  source("~/Desktop/project/research_project/microglia_visualization_topic_model.R")
  for (fit_file in fit_files_paths) {
    fit <- readRDS(fit_file)
    k <- as.integer(dim(fit$L)[2])
    path_to_plots_for_k <- glue(path_to_plots, "Plots_for_k={k}/")
    run_main_flow(obj_subset, fit, path_to_plots_for_k) #from ~/Desktop/project/research_project/microglia_visualization_topic_model.R"
  }
  
  #run correlation between the topics
  source("~/Desktop/project/research_project/microglia_topic_correlation.R")
  run_topic_evaluation(fit_files_paths, path_to_plots) #from "~/Desktop/project/research_project/microglia_topic_correlation.R"
  
  #check correlation vs cortex topics
  #todo - remove or enter to condition

  run_de_flow()
  run_correlation_with_cortex_flow()
  # generate boxplot
  run_generate_box_plot()
  run_generate_box_plot()

  #run pathways
  #...
  
}


main()

###########3

###################draft
obj <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/filtered_microglia.rds")

