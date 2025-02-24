library(tidyr)
library(tibble)
library(mombf)

tf_10kb <- read.csv("~/projects/cell-lines/data/tert_tf_jaspar_10kup.csv", header = TRUE, check.names=FALSE)
tf_1kb <- read.csv("~/projects/cell-lines/data/tert_tf_jaspar_1kup.csv", header = TRUE, check.names=FALSE)
tf_500b <- read.csv("~/projects/cell-lines/data/tert_tf_jaspar_500up.csv", header = TRUE, check.names=FALSE)
d <- read.csv("~/projects/cell-lines/data/rnaseq_merged_rsem_tpm_20250117.csv")

# Formatting
expr <- d[-(2:3), -(2:3)] %>%
  t() %>%
  as.data.frame() %>%
  setNames(make.unique(as.character(.[1,]))) %>%
  .[-1,] %>% remove_rownames() %>%
  column_to_rownames(var = "model_name")

# TFs present within 10kb or 1kb upstream of TERT genomic region  (retrieve from UCSC - Jaspar)
tf_ls_10kb <- unique(tf_10kb$TFName)
tf_ls_1kb <- unique(tf_1kb$TFName)
tf_ls_500b <- unique(tf_500b$TFName)

common <- Reduce(intersect, list(tf_ls_10kb, tf_ls_1kb, tf_ls_500b))
unique_10kb <- setdiff(tf_ls_10kb, union(tf_ls_1kb, tf_ls_500b))

tf_expr_10kb <- expr[ , colnames(expr) %in% tf_ls_10kb]
tf_expr <- expr[ , colnames(expr) %in% tf_ls_1kb]

# convert table to numeric
convert_table_to_numeric <- function(table) {
  numeric_matrix <- as.matrix(table)
  numeric_matrix <- apply(numeric_matrix, 2, as.numeric)
  rownames(numeric_matrix) <- rownames(table)
  colnames(numeric_matrix) <- colnames(table)
  return(numeric_matrix)
}

# x (TFs)
x_10kb <- convert_table_to_numeric(tf_expr_10kb)
x <- convert_table_to_numeric(tf_expr)

# y (TERT)
y <- expr[["TERT"]]

identical(length(y), nrow(x_10kb))
identical(length(y), nrow(x))
sum(is.na(y))
sum(is.na(x_10kb)) 
sum(is.na(x))

# mombf 
fit_10kb <- modelSelection(y, x_10kb, 
                      family = "normal",
                      priorCoef = momprior(), 
                      verbose = FALSE)
sampl_10kb <- rnlp(msfit = fit_10kb)
beta_10kb <- colMeans(sampl_10kb)[c(-1, -ncol(sampl_10kb))]

fit <- modelSelection(y, x, 
                      family = "normal",
                      priorCoef = momprior(), 
                      verbose = FALSE)
sampl <- rnlp(msfit = fit)
beta <- colMeans(sampl)[c(-1, -ncol(sampl))]

# model output 
output_10kb <- as.data.frame(coef(fit_10kb))
output_b_10kb <- transform(output_10kb[-c(1, nrow(output_10kb)), ], beta = beta_10kb)

output <- as.data.frame(coef(fit))
output_b <- transform(output[-c(1, nrow(output)), ], beta = beta)

write.csv(output_b_10kb, "~/projects/cell-lines/mombf/mombf-tf_10kbup_jaspar_beta.csv")
write.csv(output_b, "~/projects/cell-lines/mombf/mombf-tf_1kbup_jaspar_beta.csv")

# Filter significant TFs
filter_10kb <- output_10kb[output_10kb$margpp >= 0.5, ]
beta.fil_10kb <- beta_10kb[names(beta_10kb) %in% rownames(filter_10kb)]

filter <- output[output$margpp >= 0.5, ]
beta.fil <- beta[names(beta) %in% rownames(filter)]

write.csv(filter_10kb, "~/projects/cell-lines/mombf/mombf-tf_10kbup_jaspar_beta_sig.csv")
write.csv(filter, "~/projects/cell-lines/mombf/mombf-tf_1kbup_jaspar_beta_sig.csv")

# Plot filtered beta
model_coef_plot <- function(coef, title, filename) {
  regulators <- names(coef)
  png(filename, width = 10, height = 7, units = "in", res = 300)
  par(mar = c(7, 4, 4, 2) + 0.5, cex.axis = 0.8, mgp = c(3.5, 0.8, 0))
  bar_positions <- barplot(coef, xlab = "TFs", ylab = "Coefficients",
                           names.arg = rep("", length(coef)), las = 1, cex.names = 0.6,
                           main = title, ylim = c(min(coef) - 0.1 * abs(min(coef)), max(coef) + 0.1 * abs(max(coef))))
  text(x = bar_positions, y = par("usr")[3] - 0.1 * abs(par("usr")[3]), labels = regulators, srt = 45, adj = 1, xpd = TRUE, cex = 0.6)
  dev.off()  # Close the PNG device
}

model_coef_plot(beta.fil_10kb, "TERT regulation model_10kb", "~/projects/cell-lines/mombf/mombf-tf_10kbup_jaspar.png")
model_coef_plot(beta.fil, "TERT regulation model_1kb", "~/projects/cell-lines/mombf/mombf-tf_1kbup_jaspar.png")

# Check for GABPA
output_b_10kb["GABPA", ] # 0 beta
output_b["GABPA", ] # 0 beta
