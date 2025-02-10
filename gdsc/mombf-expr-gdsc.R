library(mombf)

# Import data
d <- read.csv(file = "~/projects/cell-lines/gdsc/expr/rnaseq_tpm_20220624.csv", header = TRUE)

expr <- d[ -(1:4), -1] 
colnames(expr) <- d[1, -1]

# Import gene list references & extract TFs
# https://humantfs.ccbr.utoronto.ca/download.php
# gene_ref from PMID:29425488, contains TFs and co-factors
gene_ref <- read.csv(file = "~/projects/cell-lines/hTFs.csv", row.names = 2)[, -1]
tf <- gene_ref[gene_ref$Is.TF. == "Yes", 1]


# y
tert <- expr[expr[, 1] == "TERT", ] 
rownames(tert) <- tert[, 1]
target <- t(tert[, -1]) # transpose the table


# x
# data formatting/cleaning
tf.activities <- expr[expr[, 1] %in% tf, ]
rownames(tf.activities) <- tf.activities[, 1]
tf.activities <- tf.activities[ ,-1]

tf.candidates <- rownames(tf.activities)
cell_ls <- colnames(tf.activities)

# convert data to number
tf.activities <- data.frame(lapply(tf.activities, 
                                   function(x) as.numeric(as.character(x))))
str(tf.activities)
rownames(tf.activities) <- tf.candidates
colnames(tf.activities) <- cell_ls

# transpose the table
tf.activities <- t(tf.activities)

identical(rownames(target), rownames(tf.activities))
sum(is.na(target))
sum(is.na(tf.activities)) # tf.activities contains NA values

# filter rows with NA
tf_clean <- tf.activities[complete.cases(tf.activities), ]
target_clean <- target[rownames(target) %in% rownames(tf_clean),]

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

beta <- mombf_model(target_clean, tf_clean)


# or
fit <- modelSelection(target_clean, tf_clean, 
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


# filter data
output <- as.data.frame(coef(fit))

# filter with margpp > 0.5
filter <- output[output$margpp >= 0.5, ]
filter <- output[output$margpp >= 0.5, ]
tf.candidates.fil <- rownames(filter)[-c(1, nrow(filter))]

beta.fil <- beta[names(beta) %in% tf.candidates.fil]

model_coef_plot.fil <-function(coef,title){
  gene_names <- tf.candidates.fil
  par(mar = c(5, 4, 4, 2) + 0.5,cex.axis=0.8, mgp = c(3.5, 0.8, 0))
  barplot(coef, xlab = "Gene Names", ylab = "Coefficients",
          names.arg = gene_names, las = 2, cex.names = 0.6,
          main = title)
}

model_coef_plot.fil(beta.fil,"TERT regulation model")


### positive control GAPBA exist
any(tf == "GABPA")
any(expr == "GABPA")
any(tf.candidates == "GABPA") 
any(colnames(tf_clean) == "GABPA") 
any(rownames(output) == "GABPA") 
any(tf.candidates.fil == "GABPA") # margpp of GABPA < 0.5

ctr <- output[rownames(output) == "GABPA", ] 

