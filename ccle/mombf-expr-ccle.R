library("tidyr")
library("dplyr")
library("viridis")
library("tibble")
library("gridExtra")
library("stringr")
library("depmap") #or library("ExperimentHub")

# Load depmap database
tpm <- depmap_TPM()

# or 
# eh <- ExperimentHub()
# query(eh, "depmap")
# eh_tpm <- eh[["EH7292"]]


# Formatting
expr <- tpm %>%
  select(cell_line, gene_name, rna_expression) %>%
  pivot_wider(names_from = gene_name, values_from = rna_expression) %>%
  column_to_rownames(var = "cell_line")


# Import gene list references & extract TFs (https://humantfs.ccbr.utoronto.ca/download.php)
gene_ref <- read.csv(file = "~/projects/cell-lines/hTFs.csv", row.names = 2)[, -1]
tf <- gene_ref[gene_ref$Is.TF. == "Yes", 1]

# y
tert <- expr[["TERT"]]


# x
tf.activities <- as.matrix(expr[ , colnames(expr) %in% tf])


# backup parameters
tf.candidates <- colnames(tf.activities)
cell_ls <- rownames(tf.activities)

# check
str(tert)
str(tf.activities)
identical(length(tert), nrow(tf.activities))
sum(is.na(tert))
sum(is.na(tf.activities)) 


# Mombf
fit <- modelSelection(tert, tf.activities, 
                      family = "normal",
                      priorCoef = momprior(), 
                      verbose = FALSE)
sampl <- rnlp(msfit=fit)
beta <- colMeans(sampl)[c(-1, -ncol(sampl))]

# plot beta 
model_coef_plot <-function(coef,title){
  gene_names <- tf.candidates
  par(mar = c(5, 4, 4, 2) + 0.1,cex.axis=0.8)
  barplot(coef, xlab = "Gene Names", ylab = "Coefficients",
          names.arg = gene_names, las = 2, cex.names = 0.6,
          main = title)
}

model_coef_plot(beta,"TERT regulation model")


# save mombf output
output <- as.data.frame(coef(fit))
write.csv(output, file = "~/projects/cell-lines/ccle/mombf-tert-ccle.csv")


# filter with margpp > 0.5
filter <- output[output$margpp >= 0.5, ]
tf.candidates.fil <- rownames(filter)[-c(1, nrow(filter))]

beta.fil <- beta[names(beta) %in% tf.candidates.fil]

model_coef_plot.fil <- function(coef, title, filename) {
  gene_names <- tf.candidates.fil
  # Open a PNG device with 300 DPI resolution
  png(filename, width = 10, height = 7, units = "in", res = 300)
  par(mar = c(5, 4, 4, 2) + 0.5, cex.axis = 0.8, mgp = c(3.5, 0.8, 0))
  # Set ylim to include both positive and negative values
  ylim_range <- range(coef, na.rm = TRUE)
  barplot(coef, xlab = "Gene Names", ylab = "Coefficients",
          names.arg = gene_names, las = 2, cex.names = 0.6,
          main = title, ylim = ylim_range)
  dev.off()  # Close the PNG device
}

model_coef_plot.fil(beta.fil, "TERT regulation model", "~/projects/cell-lines/ccle/mombf-tert-ccle.png")


# check control GAPBA 
any(rownames(output) == "GABPA") 
any(tf.candidates.fil == "GABPA") # margpp of GABPA < 0.5
# GABPA not sig. 
ctr <- output[rownames(output) == "GABPA", ] 
print(ctr)

