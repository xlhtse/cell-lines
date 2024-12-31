library(dplyr)

# Read the CSV file
d <- read.csv("../filter/Cancer-cell-line.csv")

# Function to extract the disease
disease <- function(disease) {
  parts <- unlist(strsplit(as.character(disease), "[;|]"))
  if (length(parts) >= 3) {
    return(trimws(parts[3]))
  } else {
    return(NA)
  }
}

# Apply the function to create a new column 'third_info'
d$dis <- sapply(d$disease, disease)

# Split the data frame based on the 'dis' column
CS_subset <- split(d, d$dis)

# Create a new folder named 'CS'
dir.create("CS", showWarnings = FALSE)

# Save each group to a separate CSV file
for (name in names(CS_subset)) {
  file_name <- gsub("[ /]", "-", name)
  write.csv(CS_subset[[name]], paste0("CS-subset/", file_name, ".csv"), row.names = FALSE)
}
