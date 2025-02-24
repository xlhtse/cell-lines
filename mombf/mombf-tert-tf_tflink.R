library(readr)
library(tidyr)
library(tibble)
library(mombf)

tf <- read_delim("~/projects/cell-lines/data/tert_tf_tflink.tsv")
d <- read.csv("~/projects/cell-lines/data/rnaseq_merged_rsem_tpm_20250117.csv")

# Formatting
expr <- d[-(2:3), -(2:3)] %>%
  t() %>%
  as.data.frame() %>%
  setNames(make.unique(as.character(.[1,]))) %>%
  .[-1,] %>% remove_rownames() %>%
  column_to_rownames(var = "model_name")

# TFs of TERT 
tf_ls <- unique(tf$Name.TF)

tf_expr <- expr[ , colnames(expr) %in% tf_ls]

# convert table to numeric
convert_table_to_numeric <- function(table) {
  numeric_matrix <- as.matrix(table)
  numeric_matrix <- apply(numeric_matrix, 2, as.numeric)
  rownames(numeric_matrix) <- rownames(table)
  colnames(numeric_matrix) <- colnames(table)
  return(numeric_matrix)
}

# x
x <- convert_table_to_numeric(tf_expr)


# y
y <- expr[["TERT"]]

identical(length(y), nrow(x))
sum(is.na(y))
sum(is.na(x))


# mombf 
fit <- modelSelection(y, x, 
                      family = "normal",
                      priorCoef = momprior(), 
                      verbose = FALSE)
sampl <- rnlp(msfit = fit)
beta <- colMeans(sampl)[c(-1, -ncol(sampl))]

# model output 
output <- as.data.frame(coef(fit))
output_b <- transform(output[-c(1, nrow(output)), ], beta = beta)

write.csv(output_b, "~/projects/cell-lines/mombf/mombf-tf_tflink_beta.csv")

# Filter significant TFs
filter <- output[output$margpp >= 0.5, ]
beta.fil <- beta[names(beta) %in% rownames(filter)]

write.csv(filter, "~/projects/cell-lines/mombf/mombf-tf_tflink_beta_sig.csv")

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

model_coef_plot(beta.fil, "TERT regulation model", "~/projects/cell-lines/mombf/mombf-tf_tflink.png")

# Check for GABPA
any(names(beta) == "GABPA")
output_b["GABPA", ] # 0 beta
