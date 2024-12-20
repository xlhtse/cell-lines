library(io)
library(arrow)

d <- qread("../feather/cellosaurus.feather")

table(d$category)

# extract GEO accession

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

d$encode[!is.na(d$encode)][1:10]

qwrite(d, "cellosaurus.rds")

