#!/bin/bash
# Download data 

set -euo pipefail
IFS=$'\n\t'

# OBO format file does not contain:
#   the STR profile data
#   the age at sampling
#   the full reference records (authors, title, journal)

# curl -O https://ftp.expasy.org/databases/cellosaurus/cellosaurus.obo
curl -O https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt
curl -O https://ftp.expasy.org/databases/cellosaurus/cellosaurus_refs.txt
curl -O https://ftp.expasy.org/databases/cellosaurus/cellosaurus_xrefs.txt
<<<<<<< HEAD
=======



# Cell model passports 
# Expression - RNA-Seq processed data
# Merged rnaseq_merged_20250117.zip
curl -O https://cog.sanger.ac.uk/cmp/download/rnaseq_merged_20250117.zip
unzip rnaseq_merged_20250117.zip
rm rnaseq_merged_20250117.zip

# GDSC 
# CFE in Cell Lines (https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//CFEs_in_Cell_Lines.html)
# Beta values for all CpG islands across all cell lines
# access date: 20250313 (new)
curl -O https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources///Data/preprocessed/methylation/METH_CELL_DATA.txt.zip
unzip METH_CELL_DATA.txt.zip
rm METH_CELL_DATA.txt.zip
# annotation file for mapping between cell line cosmic identifiers and methylation data sample identifiers
curl -o cell-annot.xlsx https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources////Data/otherAnnotations/methSampleId_2_cosmicIds.xlsx

# TERT genomic information retrieve from UCSU genome browswer (20250221)
# region chr5 1253167-1295068 + 10 kb upstream  ==> i.e. chr5:1253167-1305068
# Database: hg38    Primary Table: jaspar2024 Data last updated: 2024-02-09
curl -o tert_tf_jaspar_10kup.csv https://usegalaxy.org/api/datasets/f9cad7b01a472135f28d0e486221e470/display?to_ext=csv
# region chr5 1253167-1295068 + 1 kb upstream  ==> i.e. chr5:1253167-1296068
# Database: hg38    Primary Table: jaspar2024 Data last updated: 2024-02-09
curl -o tert_tf_jaspar_1kup.csv https://usegalaxy.org/api/datasets/f9cad7b01a472135c4dc4e450abcaf41/display?to_ext=csv
# region chr5 1253167-1295068 + 500 b upstream  ==> i.e. chr5:1253167-1295568
# Database: hg38    Primary Table: jaspar2024 Data last updated: 2024-02-09
curl -o tert_tf_jaspar_500up.csv https://usegalaxy.org/api/datasets/f9cad7b01a47213542b734b2fe3f29a1/display?to_ext=csv
# region chr5 1294968 + 1 kb upstream   ==> i.e. chr5:1294968-1296068
# Database: hg38    Primary Table: jaspar2024 Data last updated: 2024-02-09
curl -o tert_tf_jaspar_100d-1kbup.csv https://usegalaxy.org/api/datasets/f9cad7b01a47213557431e7bcf875ce7/display?to_ext=csv



# TERT TFLink 
# data retrieve date: 20250221
curl -o tert_tf_tflink.tsv https://cdn.netbiol.org/tflink/proteinDownload_tfs/TFLink_tfs_of_O14746.tsv
>>>>>>> 601d57c (update get.sh)
