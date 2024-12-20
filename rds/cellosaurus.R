
# ontologyIndex::get_ontology uses > 60g RAM
library(ontologyIndex)
x <- get_ontology("../data/cellosaurus.obo", extract_tags="everything");

