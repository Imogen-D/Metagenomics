#script to produce heatmap for top 20 fungal families in reindeer samples WITH CONTAMINANTS REMOVED

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
library(decontam)
library(tibble)

setwd("~/MEME/Uppsala_Katja_Project/Metagenomics") #for local script

full_otu <- read.delim("./data/reindeer_kraken2_otu_table_merged_201129-otu.fungi.txt",na.strings = c("","NA"), stringsAsFactors=FALSE)

#making OTU table
#couldn't do this filtering (and didn't seem to need to) as removed heaps of taxa
#full_otu <- read.delim("./data/reindeer_kraken2_otu_table_merged_201129-otu.fungi.txt",na.strings = c("","NA"), row.names=1, stringsAsFactors=FALSE) %>% 
 # replace(., is.na(.), 0) %>% 
  #select(which(colSums(.) > 0)) # remove taxa summing to zero
colnames(full_otu) <- str_replace(colnames(full_otu), "X", "")
OTU <- otu_table(full_otu, taxa_are_rows = FALSE)


#reading metadata table
metadata <- read.delim("./data/Sample_processing_masterlist.txt", stringsAsFactors=FALSE) %>% 
  distinct(Seq.label,.keep_all=T) %>%  # remove duplicates, the tidyr way
  filter(!is.na(Seq.label))
rownames(metadata)<-metadata$Seq.label # add rownames



set_entrez_key("ee1b29805250345f705302e643b1bfc4e007")
#making TaxonomyTable
OTUtaxa <- classification(colnames(full_otu), db = "ncbi")
bound1<-bind_rows(as_tibble(cbind(OTUtaxa))) %>%
  select(kingdom,phylum,class,order,family,genus,species)
rownames(bound1)<-names(OTUtaxa) 
write.csv(bound1, "./data/OTUtaxonomyformattedwcont.csv")

#reading and formatting taxonomy table
OTUtaxonomyformatted <- read.csv("./data/OTUtaxonomyformattedwcont.csv", row.names=1, stringsAsFactors=FALSE) %>% # read in taxa table saved from taxize 
  rename_all(str_to_title)  # make the column names into title case
taxotable <- tax_table(as.matrix(taxa_names %>% column_to_rownames("tax_id"))) # matrix required for tax table
sampledata <- sample_data(metadata[sample_names(OTU),]) # only take the samples that are present in the OTU table

# making full phyloseq data format
phydata <- phyloseq(OTU, sampledata,taxotable)


controls<-phydata@sam_data$Sample.R_cat %in% c("ExtBlank","LibBlank","Swab")
contaminants <- isContaminant(phydata, method="prevalence",neg=controls,batch="Ext.batch")
#contaminants <- isContaminant(phydata, method="frequency",conc = "Seq.copies.in.pool")
sum(contaminants$contaminant==TRUE)

nocontan<-rownames(contaminants[which(contaminants$contaminant != TRUE),])
phydata.nocontan<-phyloseq(OTU[,nocontan], sampledata,taxotable[nocontan,])
 
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