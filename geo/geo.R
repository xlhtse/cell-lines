if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)

# Load file
ESC <- read.csv("../filter/Embryonic-stem-cell.csv")
CS <- read.csv("../filter/Cancer-cell-line.csv")

# Function to get platform ID and GSE ID from GSM ID
get <- function(gsm_id) {
  gsm <- tryCatch({
    getGEO(gsm_id, GSEMatrix = FALSE)
  }, error = function(e) {
    return(list(platform_id = NA, series_id = NA))
  })
  
  if (inherits(gsm, "GSM")) {
    platform_id <- Meta(gsm)$platform_id
    series_id <- Meta(gsm)$series_id
  } else {
    platform_id <- NA
    series_id <- NA
  }
  
  return(list(platform_id = platform_id, series_id = series_id))
}

# Grab GSE IDs given by Cellosaurus without GSM
grabGSE <- function(geo_column) {
  sapply(geo_column, function(geo) {
    paste(grep("^GSE", unlist(strsplit(geo, "\\|")), value = TRUE), collapse = "|")
  })
}
data$c.gse <- grabGSE(data$geo)

# Grab GSM IDs and create the 'gsm' column
grabGSM <- function(geo_column) {
  sapply(geo_column, function(geo) {
    paste(grep("^GSM", unlist(strsplit(geo, "\\|")), value = TRUE), collapse = "|")
  })
}
data$gsm <- grabGSM(data$geo)

# Function to query GSE IDs and platform IDs from the 'gsm' column
query <- function(gsm_column) {
  gsm_ids <- unique(unlist(strsplit(paste(gsm_column, collapse = "|"), "\\|")))
  results <- lapply(gsm_ids, get)
  
  gse_list <- lapply(results, function(x) x$series_id)
  gpl_list <- lapply(results, function(x) x$platform_id)
  
  gsm_to_gse <- setNames(gse_list, gsm_ids)
  gsm_to_gpl <- setNames(gpl_list, gsm_ids)
  
  gse <- sapply(gsm_column, function(gsm) {
    if (gsm == "") return(NA)
    ids <- unlist(strsplit(gsm, "\\|"))
    unique_gse <- unique(unlist(lapply(ids, function(id) gsm_to_gse[[id]])))
    if (length(unique_gse) == 0) unique_gse <- NA
    return(paste(unique_gse, collapse = "|"))
  })
  
  gpl <- sapply(gsm_column, function(gsm) {
    if (gsm == "") return(NA)
    ids <- unlist(strsplit(gsm, "\\|"))
    unique_gpl <- unique(unlist(lapply(ids, function(id) gsm_to_gpl[[id]])))
    if (length(unique_gpl) == 0) unique_gpl <- NA
    return(unique_gpl)
  })
  
  return(list(gse = gse, gpl = gpl))
}

# Fetch GSE and GPL 
## very slow
fetchESC <- query(ESC$gsm)
fetchCS <- query(CS$gsm)

# Add gse column to data
ESC$gse <- fetchESC$gse
CS$gse <- fetchCS$gse


# Create a new table with platform IDs as columns
ESC_GPL <- fetchESC$gpl
CS_GPL <- fetchCS$gpl
ESC_platform <- unique(unlist(ESC_GPL))
CS_platform <- unique(unlist(CS_GPL))

ESC_result <- data.frame(id = ESC$id)
for (pid in ESC_platform[!is.na(ESC_platform)]) {
  ESC_result[[pid]] <- sapply(seq_len(nrow(ESC)), function(i) {
    if (!is.null(ESC_GPL[[i]]) && pid %in% ESC_GPL[[i]]) {
      ESC$gse[i]
    } else {
      NA
    }
  })
}

CS_result <- data.frame(id = CS$id)
for (pid in CS_platform[!is.na(CS_platform)]) {
  CS_result[[pid]] <- sapply(seq_len(nrow(CS)), function(i) {
    if (!is.null(CS_GPL[[i]]) && pid %in% CS_GPL[[i]]) {
      CS$gse[i]
    } else {
      NA
    }
  })
}
write.csv("ESC_result.csv")
write.csv("CS_result.csv")