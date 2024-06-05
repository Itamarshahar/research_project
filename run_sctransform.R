# check_dependencies <- function() {
#   required_packeges <- c("BiocManager", "Seurat")
#   for (package in required_packeges) {
#     if (!requireNamespace(package, quietly = TRUE)) {
#       install.packages(package)
#     }
#   }
#   # if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   #   install.packages("BiocManager")
#   # }
#   if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
#     BiocManager::install("glmGamPoi")
#   }
# }

load_libraries <- function() {
  library(Seurat)
  library(withr)
  # library(SeuratDisk)
  #library(SeuratData)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(plotly)
  # library(hrbrthemes)
  # source("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/utils/utils.R")
  library(Matrix)
  # library(fastTopics)
  library(cowplot)
  library(MASS)
  library(rhdf5)
}

# check_dependencies()
print("Loading libraries")
load_libraries()
print("Libraries loaded")
path_to_obj <- "../objects/super_agers_astrocytes.rds"
print("Reading object")
obj <- readRDS(path_to_obj)
print("Object read")
print("Running SCTransform")
obj <- SCTransform(obj = obj,
                       variable.features.n = 5000)
print("SCTransform done")
print("Saving object")
saveRDS(obj, "../objects/super_agers_astrocytes_sctransform.rds")
print("Object saved to ../objects/super_agers_astrocytes_sctransform.rds")