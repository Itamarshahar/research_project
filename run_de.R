################################################################################
## run DE
################################################################################
# Receive variables from the user (Script Only) -----------------------------------------
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  stop("please provide 5 args: <path_to_obj> <path_to_de> <path_to_plots> <K> <type_of_de>")

}
    else {
    path_to_obj <- args[1] # Obj
    path_to_fit <- args[2] # RDS obj with the fitted topic model
    path_to_output <- args[3] # path to the output (will be saved as RDS)
    LFC_type <- args[4] # the type of the calculation of the LogFoldChange (vsnull, k, le) more info in the "FastTopic packege")
    NC <- args[5] # the number of cores to use

}

set.seed(1)
counts <- obj@assays$RNA@counts
row_sum <- rowSums(counts)
nonzero_rows <- row_sum != 0
counts <- counts[nonzero_rows,]
counts <- t(counts)

de_psd <- de_analysis(fit = fit,
                  X = counts,
                  pseudocount = 0.1,
                  lfc.stat  = LFC_type,
                  control = list(ns = 1e4,nc = 6),
                  verbose = TRUE,
                  )
saveRDS(de_psd, paste0(path_to_output, "de_psd_", LFC_type, ".rds"))

