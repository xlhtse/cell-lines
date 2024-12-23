library(io)
library(arrow)
library(feather)
library(data.table)

# read feather file 
# d <- qread("../feather/cellosaurus.feather")

# read csv file
d <- fread("../feather/cellosaurus.csv")

table(d$category)

# extract accession

extract_xref <- function(d, name) {
	xrefs <- strsplit(d$xref, "|", fixed=TRUE);
	accessions <- lapply(xrefs,
		function(xref) {
			idx <- grep(name, xref);
			if (length(idx) > 0) {
				# extract accession number
				sub(sprintf("%s=([^|]*)", name), "\\1", xref[idx])
			} else {
				NA
			}
		}
	);
	# collapse into a single string
	unlist(lapply(accessions,
		function(x) {
			if (is.na(x[1])) {
				NA
			} else {
				paste(x, collapse="|")
			}
		}
	))
}


d$geo <- extract_xref(d, "GEO");
d$encode <- extract_xref(d, "ENCODE");

# check parsed results
d$geo[!is.na(d$geo)][1:10]
d$encode[!is.na(d$encode)][1:10]

# extract omics availability

extract_comment <- function(d, name) {
	comments <- strsplit(d$comment, "|", fixed=TRUE);
	accessions <- lapply(comments,
		function(cc) {
			idx <- grep(name, cc);
			if (length(cc) > 0) {
				# extract accession number
				sub(sprintf(" ?%s: ([^|]*)\\.", name), "\\1", cc[idx])
			} else {
				NA
			}
		}
	);
	accessions
}

omics <- extract_comment(d, "Omics");
omics.all <- unique(unlist(omics));
omics.all

# create omics availability matrix
idx.omics.avail <- unlist(lapply(omics, length)) > 0;
samples <- d$id[idx.omics.avail];
omics.avail <- matrix(0,
	nrow = length(samples), ncol = length(omics.all),
	dimnames = list(samples, omics.all)
);
omics.f <- omics[idx.omics.avail];
for (i in 1:nrow(omics.avail)) {
	omics.avail[ i, omics.f[[i]] ] <- 1;
}

table(rowSums(omics.avail))
sort(colSums(omics.avail))

qwrite(d, "cellosaurus.rds")
qwrite(omics.avail, "cellosaurus_omics.rds")

