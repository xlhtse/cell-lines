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
