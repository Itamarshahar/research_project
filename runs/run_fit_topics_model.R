args = commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("please provide 5 args: <path_to_the_input_data> <path_to_output_data> <K> <NC> <TRANSPOS(1/0)>")
} else {
  input_path <- args[1]
  output_path <- args[2]
  K <- as.integer(args[3])
  NC <- as.integer(args[4])
  TRANSPOSE <- as.integer(args[5])
}
left_path <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/SuperAgers_n=20/5_percent"
objects <- "objects"
obj_name <- "five_prec"
seurat_postfix <- ".h5seurat"

input_path <- glue("{left_path}/{objects}/{obj_name}{seurat_postfix}")

create_folders <- function(base_path) {
  # Create the main folder with the name "k{K}"
  main_folder_name <- glue("k{K}")
  main_folder_path <- file.path(base_path, main_folder_name)

  # Create the "k{K}" folder
  dir.create(main_folder_path)

  # Create the "objects" and "plots" subfolders inside "k{K}"
  #subfolder_names <- c(, "plots")
  objects_path <- file.path(main_folder_path, "objects")
  plots_path <- file.path(main_folder_path, "plots")

  # Create the subfolders
  dir.create(objects_path)
  dir.create(plots_path)

  # Return the path to the "objects" folder
  objects_folder_path <- objects_path
  return(objects_folder_path)
}
output_path <- create_folders(left_path)

#output_path <- glue("{left_path}/{objects}/")

K <- 20
NC <- 4
TRANSPOSE <- 1

if (!requireNamespace("fastTopics", quietly = TRUE)) {
  # Install the package
  install.packages("fastTopics", repos = cran_mirror)
}
cat(length(args))

library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
if (!requireNamespace("Seurat", quietly = TRUE)
) {
  # Install the package
  install.packages("Seurat", dependencies = TRUE, repos = 'http://cran.rstudio.com/')
}

library(Seurat)
library(SeuratDisk)

set.seed(1)
obj <- LoadH5Seurat(input_path)

################################################################################
## fetch the data: 
################################################################################

counts <- obj@assays$RNA@counts
if (TRANSPOSE == 1) {
  counts <- t(counts)
}
################################################################################
## fit the topic model 
################################################################################
print(paste0("Starting the fit with K=", K))

# Remove the gene that have only zero (in all cells)
col_sums <- colSums(counts)
nonzero_cols <- col_sums != 0
counts <- counts[, nonzero_cols]

fit <- fit_topic_model(counts,
                       k = K,
                       control.init = list(nc = NC),
                       control.main = list(nc = NC),
                       control.refine = list(nc = NC)
)

saveRDS(fit, file.path(output_path, glue("fitted_k_", K, ".rds")))
print("Finished the fit")
