##Data Quality, Read Counts##
##So need help on working out how to determine cross contaminants
## related to (rownames(top_samples) %in% rownames(top_blanks)
##Visulisations/tidy up of code

###DON'T RUN THIS CODE - PRODUCES THE RDS OBJECTS###

library(usethis)
use_git_config(user.name = "Imogen-D", user.email = "imogen.dumville@gmail.com")

library(tidyverse)
library(phyloseq)
library(reshape2)
library(microbiomeutilities)
library(decontam)
source("scripts/ancom_v2.1.R")
library(vegan)
library(nlme)
library(compositions)
library(pairwiseAdonis)
library(reshape2)

##### ABUNDANCES -  BRACKEN output #####
## phyloseq object with abundance data
full_otu <- read.delim("./data/kraken2_otu_table_merged_210216-otu.fungi.txt",na.strings = c("","NA"), stringsAsFactors=FALSE)

## OTU table
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

# only take the samples that are present in the OTU table
sampledata <- sample_data(metadata[sample_names(OTU),])

# reading and formatting taxonomy table, made with classification
##ASK ADRIAN FOR THIS TAXO INFO INFORMATION
OTUtaxonomyformatted <- read.csv("../taxa.csv",header=T,na.strings = c("NA","")) %>% 
  filter(kingdom == "Fungi") %>% 
  rename_all(str_to_title) %>%    # make the column names into title case
  select(Tax_id,Phylum:Species) # important to subset here!!!
  
taxotable <- OTUtaxonomyformatted %>% 
  filter(Tax_id %in% taxa_names(OTU))
taxotable <- tax_table(as.matrix(taxotable)) # only include the taxa that are found in the OTU table

# only use this to make row names if taxa match the OTU table
rownames(taxotable)<-OTUtaxonomyformatted$Tax_id[OTUtaxonomyformatted$Tax_id %in% taxa_names(OTU)]

# making full phyloseq data format
phydata <- phyloseq(OTU, sampledata,taxotable)

#Only extraction # not with reindeer = BE103 and bear swabs BS003, BS005 - remove
#sample_names(phydata) != 
wanted <- !(sample_names(phydata) %in% c("BE103", "BS003", "BS005"))
phydata <- prune_samples(wanted, phydata)

# saving a "base" object, for later use
saveRDS(phydata,file = "data/phyloseq-otu-base.rds")

##Making phyloseq object for read information
reads <- read.table("./data/kraken2_otu_table_merged_210216-reads.fungi.txt",na.strings = c("","NA"), stringsAsFactors=FALSE) %>% 
  select(which(colSums(.) > 0))

readcounts <- otu_table(reads, taxa_are_rows = FALSE)

# remake taxa table with read taxa
read_taxotable <- OTUtaxonomyformatted %>% 
  filter(Tax_id %in% taxa_names(OTU))
read_taxotable <- tax_table(as.matrix(read_taxotable)) # only include the taxa that are found in the OTU table

# only use this to make row names if taxa match the OTU table
rownames(read_taxotable)<-OTUtaxonomyformatted$Tax_id[OTUtaxonomyformatted$Tax_id %in% taxa_names(OTU)]

# fix column names too
colnames(read_taxotable) <- colnames(OTUtaxonomyformatted)

# slight mismatch in samples between otu and reads
read_sampledata <- sample_data(metadata[sample_names(readcounts),])

# making full phyloseq data format
readphydata <- phyloseq(readcounts, read_sampledata,read_taxotable)

#Only extraction # not with reindeer = BE103 and bear swabs BS003, BS005 - remove
readwanted <- !(sample_names(readphydata) %in% c("BE103", "BS003", "BS005"))
readphydata <- prune_samples(readwanted, readphydata)

saveRDS(readphydata,file = "data/phyloseq-read-base.rds")