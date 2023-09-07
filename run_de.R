################################################################################
## run DE
################################################################################

set.seed(1)
counts <- obj@assays$RNA@counts
row_sum <- rowSums(counts)
nonzero_rows <- row_sum != 0
counts <- counts[nonzero_rows,]
counts <- t(counts)

de_psd <- de_analysis(fit = fit,
                  X = counts,
                  pseudocount = 0.1,
                  lfc.stat  = "le",
                  control = list(ns = 1e4,nc = 6),
                  verbose = TRUE,
                  )
setwd("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/databases/")
saveRDS(de_psd, "de_night")

shmuel