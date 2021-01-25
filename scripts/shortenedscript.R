#script to produce heatmap for top 20 fungal families in reindeer samples WITH CONTAMINANTS REMOVED

#load packages
library(usethis)
use_git_config(user.name = "Imogen-D", user.email = "imogen.dumville@gmail.com")
library(dplyr)
library(phyloseq)
library(tidyr)
library(ggplot2)
library(microbiomeutilities)
library(decontam)
source("scripts/ancom_v2.1.R")
library(vegan)

library(nlme)
library(tidyverse)
library(compositions)

full_otu <- read.delim("./data/reindeer_kraken2_otu_table_merged_201129-otu.fungi.txt",na.strings = c("","NA"), stringsAsFactors=FALSE) %>% 
  select(which(colSums(.) > 0))  # remove empty taxa
# filter(rowSums(.) > 0) # remove empty samples

colnames(full_otu) <- str_replace(colnames(full_otu), "X", "")
OTU <- otu_table(full_otu, taxa_are_rows = FALSE)

#reading metadata table
metadata <- read.delim("./data/Sample_processing_masterlist.txt", stringsAsFactors=FALSE) %>% 
  select(-starts_with("X")) %>% 
  distinct(Seq.label,.keep_all=T) %>%  # remove duplicates, the tidyr way
  filter(!is.na(Seq.label)) %>% 
  mutate(Reindeer.ecotype=ifelse(Sample.R_cat == "Reindeer_pre","pre-historic", # fill in some of the blank ecotypes
                                 ifelse(Sample.R_cat %in% c("ExtBlank","LibBlank","Swab"),"blank",
                                        ifelse(Sample.R_cat == "" & is.na(Sample.R_cat),"unknown",Reindeer.ecotype))))
rownames(metadata)<-metadata$Seq.label # add rownames

#reading and formatting taxonomy table, made with classification
OTUtaxonomyformatted <- read.csv("./data/OTUtaxonomyformattedwcont.csv", row.names=1, stringsAsFactors=FALSE) %>% # read in taxa table saved from taxize 
  rename_all(str_to_title)  # make the column names into title case
taxotable <- tax_table(as.matrix(OTUtaxonomyformatted))

sampledata <- sample_data(metadata[sample_names(OTU),]) # only take the samples that are present in the OTU table

# making full phyloseq data format
phydata <- phyloseq(OTU, sampledata,taxotable)

##### CONTAMINANT FILTERING #####
# create single variable for control samples
controls<-phydata@sam_data$Sample.R_cat %in% c("ExtBlank","LibBlank","Swab")

# both prevalance and frequency based filtering
both.contaminants <- isContaminant(phydata, method="either",neg = controls,conc = "Seq.copies.in.pool",threshold = c(0.1,0.5))
table(both.contaminants$contaminant)

# use phyloseq to prune taxa instead
phywocont <- prune_taxa(both.contaminants$contaminant==FALSE,phydata)

# only pulled the reindeer samples for the graphic without contaminants, 
# need to do the same for the contaminant only graphic

keep_samples<-rownames(sample_data(phywocont)[grepl(sample_data(phywocont)$Seq.label, pattern = "Rt")])
ecophycont <- prune_samples(keep_samples,phywocont)

meta_data = as.data.frame(sample_data(phywocont))
feature_table = as.data.frame(otu_table(phywocont))
main_var = "Reindeer.ecotype"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL
feature_table = otu_table(phywocont)
feature_table <- t(as.data.frame(feature_table))
meta_data = sample_data(phywocont)
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, main_var, struc_zero = NULL, p_adj_method, alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start # around 30s


