###### R ENVIRONMENT INIT #####
# First time, initiate the conda environment
# conda create -n microbiomes-env -c conda-forge r-base=4.0.2 r-devtools r-biocmanager r-renv 
# all other times after this, conda activate metagenomics-env

# library(renv)
# renv::init()
# renv::snapshot()

renv::restore()