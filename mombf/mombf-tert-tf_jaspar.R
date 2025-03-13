library(tidyr)
library(tibble)
library(mombf)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 1) {
  stop("Please provide exactly one argument: <file_path>")
}

<<<<<<< HEAD
# Assign argument to variable
file_path <- args[[1]]

# Print the file path
cat("File path:", file_path, "\n")
=======
# TFs present within 10kb or 1kb upstream of TERT genomic region  (retrieve from UCSC - Jaspar)
tf_ls_10kb <- unique(tf_10kb$TFName[tf_10kb$strand == "-"])
tf_ls_1kb <- unique(tf_1kb$TFName[tf_1kb$strand == "-"])
tf_ls_500b <- unique(tf_500b$TFName[tf_500b$strand == "-"])

print(length(tf_ls_10kb))
print(length(tf_ls_1kb))
print(length(tf_ls_500b))

common <- Reduce(intersect, list(tf_ls_10kb, tf_ls_1kb, tf_ls_500b))
unique_10kb <- setdiff(tf_ls_10kb, union(tf_ls_1kb, tf_ls_500b))
print(unique_10kb)
>>>>>>> origin/main

# Read the CSV file
tf <- read.csv(file_path, header = TRUE, check.names = FALSE)

# Load expression data from cell model passport (modified from rnaseq_merged_rsem_tpm_20250117.csv)
expr <- read.csv("../data/expr.csv")

# TFs present within TERT genomic region (retrieve from UCSC - Jaspar)
tf_ls <- unique(tf$TFName)
tf_expr <- expr[ , colnames(expr) %in% tf_ls]

# convert table to numeric
convert_table_to_numeric <- function(table) {
  numeric_matrix <- as.matrix(table)
  numeric_matrix <- apply(numeric_matrix, 2, as.numeric)
  rownames(numeric_matrix) <- rownames(table)
  colnames(numeric_matrix) <- colnames(table)
  return(numeric_matrix)
}

# x (TFs)
x <- convert_table_to_numeric(tf_expr)

# y (TERT)
y <- expr[["TERT"]]

# data check
if (!identical(length(y), nrow(x))) {
  stop("length(y) is not equal to nrow(x)")
} 
if (sum(is.na(y)) != 0) {
  stop("There is value missing in y")
}
if (sum(is.na(x)) != 0) {
  stop("There is value missing in x")
} else {
  cat("data check pass\n")
}

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

# Generate output file name with suffix "_beta"
output_file <- sub("\\.csv$", "_beta.csv", basename(file_path))

write.csv(output_b, output_file)

# Filter significant TFs
filter <- output_b[output_b$margpp >= 0.5, ]
beta.fil <- beta[names(beta) %in% rownames(filter)]

# Generate filtered output file name with suffix "_beta_sig"
filtered_output_file <- sub("\\.csv$", "_beta_sig.csv", basename(file_path))

write.csv(filter, filtered_output_file)

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

# Generate plot file name with suffix "_beta.png"
plot_file <- sub("\\.csv$", ".png", basename(file_path))

model_coef_plot(beta.fil, "TERT regulation model", plot_file)

# Check for GABPA
output_b["GABPA", ] # 0 beta