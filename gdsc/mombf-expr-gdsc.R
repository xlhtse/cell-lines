library(mombf)

# Import data
d <- read.csv(file = "~/projects/cell-lines/gdsc/expr/rnaseq_tpm_20220624.csv", header = TRUE)

# data formatting
expr <- d[-(1:4), -1] 
colnames(expr) <- d[1, -1]
expr <- t(expr)
colnames(expr) <- expr[1,]
expr <- expr[-1,]

expr_num <- as.data.frame(apply(expr, 2, as.numeric))
rownames(expr_num) <- rownames(expr)

# Import gene list references & extract TFs
# https://humantfs.ccbr.utoronto.ca/download.php
# gene_ref from PMID:29425488, contains TFs and co-factors
gene_ref <- read.csv(file = "~/projects/cell-lines/hTFs.csv", row.names = 2)[, -1]
tf <- gene_ref[gene_ref$Is.TF. == "Yes", 1]


# y
tert <- expr_num[ ,"TERT", drop = FALSE]


# x
tf.activities <- expr_num[ , colnames(expr_num) %in% tf]

# backup parameters
tf.candidates <- colnames(tf.activities)
cell_ls <- rownames(tf.activities)

# check
identical(length(tert), nrow(tf.activities))
sum(is.na(tert))
sum(is.na(tf.activities)) # tf.activities contains NA values

# filter rows with NA
tf_clean <- as.matrix(tf.activities[complete.cases(tf.activities), ])
tert_clean <- tert[rownames(tert) %in% rownames(tf_clean), "TERT"]


# Mombf
# mombf model fit -- function adapted from celine
mombf_model <- function(target, tf.activities) {
  if (length(tf.candidates) == 1) {
    # fallback to univariate model when there is only one candidate
    ulinear_model(target, tf.activities);
  } else {
    require(mombf)
    msfit <- modelSelection(target, tf.activities,
                            family="normal", priorCoef=momprior(), verbose=FALSE);
    sampl <- rnlp(msfit=msfit);
    # remove intercept and phi parameter
    beta <- colMeans(sampl)[c(-1, -ncol(sampl))];
    beta
  }
}

beta <- mombf_model(tert_clean, tf_clean)

# or
fit <- modelSelection(tert_clean, tf_clean, 
                      family = "normal",
                      priorCoef = momprior(), 
                      verbose = FALSE)
sampl <- rnlp(msfit=fit)
beta <- colMeans(sampl)[c(-1, -ncol(sampl))]


# plot beta -- function adapted from celine
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
write.csv(output, file = "~/projects/cell-lines/gdsc/mombf-tert-gdsc.csv")


# filter with margpp > 0.5
filter <- output[output$margpp >= 0.5, ]
tf.candidates.fil <- rownames(filter)[-c(1, nrow(filter))]

beta.fil <- beta[names(beta) %in% tf.candidates.fil]

model_coef_plot.fil <- function(coef, title, filename) {
  gene_names <- tf.candidates.fil
  png(filename)  # Open a PNG device
  par(mar = c(5, 4, 4, 2) + 0.5, cex.axis = 0.8, mgp = c(3.5, 0.8, 0))
  barplot(coef, xlab = "Gene Names", ylab = "Coefficients",
          names.arg = gene_names, las = 2, cex.names = 0.6,
          main = title)
  dev.off()  # Close the PNG device
}

model_coef_plot.fil(beta.fil, "TERT regulation model", "~/projects/cell-lines/gdsc/mombf-tert-gdsc.png")


# check for GABPA control
any(rownames(output) == "GABPA") 
any(tf.candidates.fil == "GABPA") # margpp of GABPA < 0.5
# GABPA not sig. 
ctr <- output[rownames(output) == "GABPA", ] 
print(ctr)
