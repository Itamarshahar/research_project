library(pheatmap)
# install.packages('wordspace')
library(wordspace)
matrix <- matrix(c(
  0.06, 0.04, 0.05, 0.06, 0.03,
  0.03, 0.03, 0.04, 0.05, 0.03,
  0.07, 0.16, 0.01, 0.01, 0.01,
  0.02, 0.02, 0.01, 0.02, 0.01,
  0.01, 0.02, 0.02, 0.03, 0.03,
  0.17, 0.12, 0.05, 0.05, 0.01,
  0.05, 0.06, 0.17, 0.18, 0.13,
  0.04, 0.03, 0.02, 0.04, 0.01,
  0.08, 0.12, 0.11, 0.07, 0.41,
  0.03, 0.04, 0.05, 0.03, 0.02,
  0.01, 0.02, 0.09, 0.05, 0.01,
  0.09, 0.06, 0.14, 0.05, 0.06,
  0.24, 0.18, 0.08, 0.05, 0.03,
  0.04, 0.02, 0.04, 0.2, 0.17,
  0.08, 0.09, 0.12, 0.13, 0.02
  ), nrow = 15, byrow = TRUE)

column_sums <- rowSums(matrix)
matrix <- matrix / column_sums
rownames(matrix) <- paste0("k", 1:15)
colnames(matrix) <- c("AD", "MCI", "HA", "SuperAgers", "Young CTRL")
pheatmap(matrix,
         cluster_rows = T,
         cluster_cols = F,
         display_numbers = TRUE,
         color = colorRampPalette(c("white", "red"))(50),
         main = "Heatmap of Mean Values",
         fontsize = 20,
         angle_col= 45)



# "k1" <-
# Mean= 0.06
# Mean= 0.04
# Mean= 0.05
# Mean= 0.06
# Mean= 0.03
# k2 <-
#   Mean= 0.03
# Mean= 0.03
# Mean= 0.04
# Mean= 0.05
# Mean= 0.03
#
#
# k3
# Mean= 0.07
# Mean= 0.16
# Mean= 0.01
# Mean= 0.01
# Mean= 0.01
#
#
# k4
# Mean= 0.02
# Mean= 0.02
# Mean= 0.01
# Mean= 0.02
# Mean= 0.01
#
#
# k5
# Mean= 0.01
# Mean= 0.02
# Mean= 0.02
# Mean= 0.03
# Mean= 0.03
#
# k6
# Mean = 0.17
# Mean= 0.12
# Mean= 0.05
# Mean= 0.05
# Mean= 0.01
#
# k7
# Mean= 0.05
# Mean= 0.06
# Mean= 0.17
# Mean= 0.18
# Mean= 0.13
#
# k8
# Mean= 0.04
# Mean= 0.03
# Mean= 0.02
# Mean= 0.04
# Mean= 0.01
#
#
# k9
#
#
# Mean= 0.08
# Mean= 0.12
# Mean= 0.11
# Mean= 0.07
# Mean= 0.41
#
# k10
# Mean= 0.03
# Mean= 0.04
# Mean= 0.05
# Mean= 0.03
# Mean= 0.02
#
# k11
# Mean= 0.01
# Mean= 0.02
# Mean= 0.09
# Mean= 0.05
# Mean= 0.01
#
# k12
# Mean= 0.09
# Mean= 0.06
# Mean= 0.14
# Mean= 0.05
# Mean= 0.06
#
#
# k13
# Mean= 0.24
# Mean= 0.18
# Mean= 0.08
# Mean= 0.05
# Mean= 0.03
#
# k14
# Mean= 0.04
# Mean= 0.02
# Mean= 0.04
# Mean= 0.2
# Mean= 0.17
#
# k15
# Mean= 0.08
# Mean= 0.09
# Mean= 0.12
# Mean= 0.13
# Mean= 0.02