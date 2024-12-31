if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)
library(tidyr)
library(dplyr)
library(parallel)


# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("No input file provided.")

# Read input file
input_file <- args[1]
d <- read.csv(input_file)


# Grab GSM IDs from the "geo" column 
grabGSM <- function(geo_column) {
  sapply(geo_column, function(geo) {
    paste(grep("^GSM", unlist(strsplit(geo, "\\|")), value = TRUE), collapse = "|")
  })
}
d$gsm <- grabGSM(d$geo)

# Select cell names and GSM IDs, unlist
ls <- d %>%
  select(id, gsm) %>%
  separate_rows(gsm, sep = "\\|") 

# Functions to fetch GSE and GPL IDs from GSM IDs
get <- function(gsm_id) {
  gsm <- tryCatch(getGEO(gsm_id, GSEMatrix = FALSE), error = function(e) return(NULL))
  if (is.null(gsm)) return(list(platform_id = NA, series_id = NA))
  list(platform_id = Meta(gsm)$platform_id, series_id = Meta(gsm)$series_id)
}

query <- function(gsm_column) {
  gsm_ids <- unique(gsm_column)
  results <- lapply(gsm_ids, get)
  cl <- makeCluster(6)
  clusterExport(cl, c("get", "getGEO", "Meta"))
  gse_list <- setNames(lapply(results, `[[`, "series_id"), gsm_ids)
  gpl_list <- setNames(lapply(results, `[[`, "platform_id"), gsm_ids)
  gse <- sapply(gsm_column, function(gsm) paste(unique(gse_list[[gsm]]), collapse = "|"))
  gpl <- sapply(gsm_column, function(gsm) unique(gpl_list[[gsm]]))
  list(gse = gse, gpl = gpl)
}

# Fetch the GSE and GPL IDs
fetch <- query(ls$gsm)
ls$gse <- fetch$gse
ls$gpl <- fetch$gpl
ls <- ls %>%
  separate_rows(gpl, sep = "\\|") %>%
  separate_rows(gse, sep = "\\|")


# Functions to fetch platform description from GPL IDs
get_pl <- function(gpl_id) {
  gpl <- tryCatch(getGEO(gpl_id, GSEMatrix = FALSE), error = function(e) return(NULL))
  if (is.null(gpl)) return(list(gpl = gpl_id, platform = NA))
  list(gpl = gpl_id, platform = Meta(gpl)$title)
}

query_pl <- function(gpl_column) {
  gpl_ids <- unique(gpl_column)
  cl <- makeCluster(6)
  clusterExport(cl, c("get_pl", "getGEO", "Meta"))
  results_pl <- lapply(gpl_ids, get_pl)
  pl_list <- setNames(lapply(results_pl, `[[`, "platform"), gpl_ids)
  sapply(gpl_column, function(gpl) paste(unique(pl_list[[gpl]]), collapse = "|"))
}

ls$platform <- query_pl(ls$gpl)
ls <- ls %>%
  separate_rows(platform, sep = "\\|")



# Functions to fetch experiment type from GSE IDs
get_exp <- function(gse_id) {
  gse <- tryCatch(getGEO(gse_id, GSEMatrix = FALSE), error = function(e) return(NA))
  if (inherits(gse, "GSE")) {
    exp_type <- Meta(gse)$type
    if (is.null(exp_type)) exp_type <- NA
  } else {
    exp_type <- NA
  }
  return(as.character(exp_type))
}

query_exp <- function(gse_column) {
  gse_ids <- unique(gse_column)
  cl <- makeCluster(6)
  clusterExport(cl, c("get_exp", "getGEO", "Meta"))
  results <- parLapply(cl, gse_ids, get_exp)
  stopCluster(cl)
  gse_to_type <- setNames(results, gse_ids)
  type <- sapply(gse_column, function(gse) {
    if (is.na(gse) || gse == "") return(NA)
    ids_gse <- unlist(strsplit(gse, "\\|"))
    type_values <- unlist(lapply(ids_gse, function(id) gse_to_type[[id]]))
    type_values <- ifelse(is.na(type_values), NA, as.character(type_values))
    return(paste(type_values, collapse = "|"))
  })
  return(type)
}

ls$exp <- query_exp(ls$gse)
ls <- ls %>%
  separate_rows(exp, sep = "\\|")
platform_info <- ls[ ,c("gpl","platform","exp")] %>% distinct()


head(ls)  
head(platform_info)

# Export 
output_file <- sub("\\.csv$", "_info.csv", input_file)
write.csv(ls, file = output_file)
file_name <- sub("\\.csv$", "_platform_info.csv", input_file)
write.csv(platform_info, file = "platform_info.csv")



