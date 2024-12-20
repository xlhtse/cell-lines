#!/bin/bash
# Download data 

set -euo pipefail
IFS=$'\n\t'

#curl -O https://ftp.expasy.org/databases/cellosaurus/cellosaurus.obo
curl -O https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt
