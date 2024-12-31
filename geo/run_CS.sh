#!/bin/bash

# Specify the folder containing the input files
INPUT_FOLDER="../filter/CS"

# Loop through each CSV file in the folder
for INPUT_FILE in "$INPUT_FOLDER"/*.csv
do
  # Execute the R script with the input file as an argument
  Rscript get-info.R "$INPUT_FILE"
done