###### R ENVIRONMENT INIT #####
# First time, initiate the conda environment
# conda env create --file environment.yaml
# all other times after this, conda activate metagenomics-env

# missing packages in yaml
BiocManager::install(pkgs = c("phyloseq","Biostrings","microbiome"))
devtools::install_github("microsud/microbiomeutilities")
