# Cell model passports (GDSC)
# https://www.cancerrxgene.org/downloads/bulk_download

#!/bin/bash

# Create directories
mkdir -p gdsc

# Array of URLs
urls=(
    "https://cog.sanger.ac.uk/cmp/download/mutations_wes_vcf_20221010.zip"
        "https://cog.sanger.ac.uk/cmp/download/mutations_wgs_vcf_20221123.zip"
	    "https://cog.sanger.ac.uk/cmp/download/mutations_tgs_vcf_20211124.zip"
	        "https://cog.sanger.ac.uk/cmp/download/mutations_all_20230202.zip"
		    "https://cog.sanger.ac.uk/cmp/download/rnaseq_all_20220624.zip"
		        "https://cog.sanger.ac.uk/cmp/download/WES_pureCN_CNV_genes_20221213.zip"
			    "https://cog.sanger.ac.uk/cmp/download/WGS_purple_CNV_genes_20230303.zip"
			        "https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/preprocessed/methylation/METH_CELL_DATA.txt.zip"
				    "https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/suppData/TableS2J.xlsx"
				        "https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/BEMs/CellLines/CellLines_METH_BEMs.zip"
				)

				# Array of filenames
				filenames=(
				    "gdsc/mut-annotated-wes-vfc.zip"
				        "gdsc/mut-annotated-wgs-vfc.zip"
					    "gdsc/mut-annotated-tgs-vfc.zip"
					        "gdsc/mut-all.zip"
						    "gdsc/rnaseq.zip"
						        "gdsc/copy-no-wes-pureCN.zip"
							    "gdsc/copy-no-wgs-PURPLE.zip"
							        "gdsc/met-beta.zip"
								    "gdsc/met-CpG-ls.zip"
								        "gdsc/met-iCpG-BEM.zip"
								)

								# Loop through URLs and filenames
								for i in "${!urls[@]}"; do
									    curl -o "${filenames[$i]}" "${urls[$i]}"
								    done



								    # CCLE dataset download
								    # https://depmap.org/portal/data_page/?tab=allData

								    mkdir -p ccle

								    # Array of URLs
								    urls=(
								        "https://ndownloader.figshare.com/files/51065489"
									    "https://ndownloader.figshare.com/files/51065492"
									        "https://ndownloader.figshare.com/files/51065495"
										    "https://ndownloader.figshare.com/files/51065732"
										        "https://ndownloader.figshare.com/files/51065324"
											    "https://ndownloader.figshare.com/files/51065303"
											        "https://ndownloader.figshare.com/files/51065297"
											)

											# Array of filenames
											filenames=(
											    "ccle/expr-prot-coding-genes-tpm-logp1.csv"
											        "ccle/expr-prot-coding-genes-tpm-logp1batchcorr.csv"
												    "ccle/expr-prot-coding-genes-tpm-logp1stranded.csv"
												        "ccle/somatic-mut.csv"
													    "ccle/copy-no.csv"
													        "ccle/copy-no-absolute.csv"
														    "ccle/cell-annot.csv"
													    )

													    # Loop through URLs and filenames
													    for i in "${!urls[@]}"; do
														        curl -o "${filenames[$i]}" "${urls[$i]}"
														done 
