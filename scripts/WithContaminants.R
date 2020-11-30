#script to produce heatmap for top 20 fungal families in reindeer samples

library(usethis)
use_git_config(user.name = "Imogen-D", user.email = "imogen.dumville@gmail.com")
library(dplyr)
library(phyloseq)
library(rentrez)
library(taxize)
library(stringr)
library(tidyr)
library(data.table)
library(ggplot2)
library(microbiomeutilities)

setwd("~/MEME/Uppsala_Katja_Project/Metagenomics") #for local script

#making OTU table
full_otu <- read.delim("./data/reindeer_kraken2_otu_table_merged_201112-otu.fungi.txt",na.strings = c("","NA"), row.names=1, stringsAsFactors=FALSE) %>% 
  replace(., is.na(.), 0) %>% 
  select(which(colSums(.) > 0)) # remove taxa summing to zero
colnames(full_otu) <- str_replace(colnames(full_otu), "X", "")
OTU <- otu_table(full_otu, taxa_are_rows = FALSE)

#reading metadata table
metadata <- read.delim("./data/reindeer_sample_metadate_merged.txt", stringsAsFactors=FALSE) %>% 
  distinct(Seq.label,.keep_all=T) %>%  # remove duplicates, the tidyr way
  filter(!is.na(Seq.label))
rownames(metadata)<-metadata$Seq.label # add rownames

#Script to produce taxonomy table (in combination with allOTUcode.R)

set_entrez_key("ee1b29805250345f705302e643b1bfc4e007")

#making TaxonomyTable
OTUtaxa <- classification(colnames(full_otu), db = "ncbi")

bound1<-bind_rows(as_tibble(cbind(OTUtaxa))) %>%
  select(kingdom,phylum,class,order,family,genus,species)
rownames(bound1)<-names(OTUtaxa)

write.csv(bound1, "./data/OTUtaxonomyformatted.csv")

#reading and formatting taxonomy table
OTUtaxonomyformatted <- read.csv("./data/OTUtaxonomyformatted.csv", row.names=1, stringsAsFactors=FALSE) %>% # read in taxa table saved from taxize 
  rename_all(str_to_title)  # make the column names into title case
taxotable <- tax_table(as.matrix(OTUtaxonomyformatted)) # matrix required for tax table
sampledata <- sample_data(metadata[sample_names(OTU),]) # only take the samples that are present in the OTU table

#maing full phyloseq data format
phydata <- phyloseq(OTU, sampledata,taxotable)

#subsetting for quicker analysis removing those without ecotype information
noeco <- (which(is.na(metadata$Reindeer.ecotype)))
ecotypemeta <- sample_data(metadata[-c(noeco),])
ecophy <- phyloseq(OTU, ecotypemeta,taxotable)

#aggregating at family level
OTUfam <- microbiome::aggregate_taxa(ecophy, "Family")

#creating heatmap with clr transformation, top 20 families
pdf(file = "./images/heatmap20familyclr.pdf", height = 5, width = 10)
plot_taxa_heatmap(OTUfam, subset.top = 20, transformation = "clr",
                  taxonomic.level = "Family", border_color = "grey60",
                  VariableA = "Reindeer.ecotype")
dev.off()


OTUfam <- microbiome::aggregate_taxa(ecophy, "Family")

#creating heatmap with clr transformation, top 20 families
pdf(file = "./images/heatmap20familyclr.pdf", height = 5, width = 10)
plot_taxa_heatmap(OTUfam, subset.top = 20, transformation = "clr",
                  taxonomic.level = "Family", border_color = "grey60",
                  VariableA = "Reindeer.ecotype")
dev.off()