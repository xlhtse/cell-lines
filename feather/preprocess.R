library(io)
library(arrow)

d <- qread("cellosaurus.feather")

# cell lines with GEO access
idx <- grep("GEO", d$xref);
d[idx[1:10], ]

table(d$category)
