library(filenamer)
library(io)

d <- qread("../rds/cellosaurus.rds");
omics <- qread("../rds/cellosaurus_omics.rds");

<<<<<<< HEAD
=======

>>>>>>> 218ed2e335e15b95579560fc8b85770d9ec95b4a
# Filter cell lines that are human species and with GEO data
h.cell <- d[d$species == "NCBI_TaxID=9606; ! Homo sapiens (Human)" & 
                 !is.na(d$geo), c("id", "accession", "synonyms", "category", "disease", "origin", "geo", "encode")]

# Subset omics to keep only the cells present in h.cell
id_omics <- intersect(h.cell$id, row.names(omics))
h.omics <- omics[row.names(omics) %in% id_omics, ]

# Summarise omics 
sum <- apply(h.omics, 1, function(x) paste(colnames(h.omics)[x == 1], collapse = "|"))
sum <- data.frame(id = row.names(h.omics), omics = sum, row.names = NULL)

<<<<<<< HEAD
# Combine all information and export csv file
annotation <- merge(h.cell, sum, by = "id", all.x = TRUE)
write.csv(annotation, file = "../filter/annotation.csv")

# generate list for each cell line type and export csv file
subset <- split(annotation, annotation$category)
=======
# Combine all information 
annotation <- merge(h.cell, sum, by = "id", all.x = TRUE)

# generate list for each cell line type
subset <- split(annotation, c$category)
>>>>>>> 218ed2e335e15b95579560fc8b85770d9ec95b4a

for (category_name in names(subset)) {
  category <- gsub(" ", "-", category_name)  # Replace whitespace with underscores
  assign(category, subset[[category_name]])  # Assign subset as a new table with the modified name
  write.csv(subset[[category_name]], file = paste0(category, ".csv"), row.names = FALSE)  # Export the table as a CSV file
}
<<<<<<< HEAD
=======





>>>>>>> 218ed2e335e15b95579560fc8b85770d9ec95b4a
