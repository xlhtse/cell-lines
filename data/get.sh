#!/bin/bash
# Download data 

set -euo pipefail
IFS=$'\n\t'


# Cellosaurus 
# last modified 2024-12-19 11:02

# OBO format file does not contain:
#   the STR profile data
#   the age at sampling
#   the full reference records (authors, title, journal)

# curl -O https://ftp.expasy.org/databases/cellosaurus/cellosaurus.obo
curl -O https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt
curl -O https://ftp.expasy.org/databases/cellosaurus/cellosaurus_refs.txt
curl -O https://ftp.expasy.org/databases/cellosaurus/cellosaurus_xrefs.txt



# Cell model passports 
# Expression - RNA-Seq processed data
# Merged rnaseq_merged_20250117.zip
curl -O https://cog.sanger.ac.uk/cmp/download/rnaseq_merged_20250117.zip
unzip rnaseq_merged_20250117.zip
rm rnaseq_merged_20250117.zip



# TERT genomic information retrieve from UCSU genome browswer (20250221)
# region chr5 1253167-1295068 + 10 kb upstream  ==> i.e. chr5: 1253167-1305068
# Database: hg38    Primary Table: jaspar2024 Data last updated: 2024-02-09
curl -o tert_tf_jaspar_10kup.csv https://usegalaxy.org/api/datasets/f9cad7b01a472135f28d0e486221e470/display?to_ext=csv
# region chr5 1253167-1295068 + 1 kb upstream  ==> i.e. chr5: 1253167-1296068
# Database: hg38    Primary Table: jaspar2024 Data last updated: 2024-02-09
curl -o tert_tf_jaspar_1kup.csv https://usegalaxy.org/api/datasets/f9cad7b01a472135c4dc4e450abcaf41/display?to_ext=csv
# region chr5 1253167-1295068 + 100 b upstream  ==> i.e. chr5: 1253167-1295568
# Database: hg38    Primary Table: jaspar2024 Data last updated: 2024-02-09
curl -o tert_tf_jaspar_500up.csv https://usegalaxy.org/api/datasets/f9cad7b01a47213542b734b2fe3f29a1/display?to_ext=csv


# TERT TFLink 
# data retrieve date: 20250221
curl -o tert_tf_tflink.tsv https://cdn.netbiol.org/tflink/proteinDownload_tfs/TFLink_tfs_of_O14746.tsv
