#!/bin/bash 

# Specify the input file 
INPUT_FILE="../filter/Embryonic-stem-cell.csv" 

# Execute the R script with the input file as an argument 
Rscript get-info.R "$INPUT_FILE"