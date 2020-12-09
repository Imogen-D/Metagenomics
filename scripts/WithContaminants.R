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

full_otu <- read.delim("./data/reindeer_kraken2_otu_table_merged_201129-otu.fungi.txt",na.strings = c("","NA"), stringsAsFactors=FALSE) %>% 
  select(which(colSums(.) > 0)) %>% # remove empty taxa
  filter(rowSums(.) > 0) # remove empty samples

colnames(full_otu) <- str_replace(colnames(full_otu), "X", "")
OTU <- otu_table(full_otu, taxa_are_rows = FALSE)

#reading metadata table
metadata <- read.delim("./data/Sample_processing_masterlist.txt", stringsAsFactors=FALSE) %>% 
  select(-starts_with("X")) %>% 
  distinct(Seq.label,.keep_all=T) %>%  # remove duplicates, the tidyr way
  filter(!is.na(Seq.label))
rownames(metadata)<-metadata$Seq.label # add rownames


#reading and formatting taxonomy table, made with classification

OTUtaxonomyformatted <- read.csv("./data/OTUtaxonomyformattedwcont.csv", row.names=1, stringsAsFactors=FALSE) %>% # read in taxa table saved from taxize 
  rename_all(str_to_title)  # make the column names into title case
taxotable <- tax_table(as.matrix(OTUtaxonomyformatted))

sampledata <- sample_data(metadata[sample_names(OTU),]) # only take the samples that are present in the OTU table

# making full phyloseq data format
phydata <- phyloseq(OTU, sampledata,taxotable)

controls<-phydata@sam_data$Sample.R_cat %in% c("ExtBlank","LibBlank")
phydata@sam_data$Ext.batch<-as.factor(phydata@sam_data$Ext.batch)

controls<-phydata@sam_data$Sample.R_cat %in% c("ExtBlank","LibBlank","Swab")
contaminants <- isContaminant(phydata, method="prevalence",neg=controls, threshold = 0.5)
contaminants <- isContaminant(phydata, method="either",neg = controls,conc = "Seq.copies.in.pool")
conts <- c(which(contaminants$contaminant==TRUE))


otuswocont <- full_otu[,-c(conts)]
OTUwocont <- otu_table(otuswocont, taxa_are_rows = FALSE)
sampledatawocont <- sample_data(metadata[sample_names(OTUwocont),]) # only take the samples that are present in the OTU table

phywocont <- phyloseq(OTUwocont, sampledatawocont, taxotable)



#only pulled the reindeer samples for the graphic without contaminants, 
#need to do the same for the contaminant only graphic
reindeersamples <- sample_data(sampledata[c(str_which(sampledata$Seq.label, pattern = "Rt")),])

#could maybe just subset out ecotype samples again, but then only 8??
#noeco <- (which(metadata$Reindeer.ecotype == "")) #39, some with weird names
#ecotypemeta <- sample_data(metadata[-c(noeco),])

ecophycont <- phyloseq(OTUwocont, reindeersamples, taxotable)

#aggregating at family level
OTUfamcont <- microbiome::aggregate_taxa(ecophycont, "Family")

#creating heatmap with clr transformation, top 20 families
pdf(file = "./images/heatmap20familycont.pdf", height = 5, width = 10)
plot_taxa_heatmap(OTUfamcont, subset.top = 20, transformation = "clr",
                  taxonomic.level = "Family", border_color = "grey60",
                  VariableA = "Reindeer.ecotype")
dev.off()
 

#looking only at contamination taxa
otusWcont <- full_otu[,c(conts)]
OTUWcont <- otu_table(otusWcont, taxa_are_rows = FALSE)
sampledataWcont <- sample_data(metadata[sample_names(OTUWcont),]) # only take the samples that are present in the OTU table

phyWcont <- phyloseq(OTUWcont, sampledataWcont, taxotable)

OTUfamcont <- microbiome::aggregate_taxa(phyWcont, "Family")

#creating heatmap with clr transformation, top 20 families OF CONTAMINATION
pdf(file = "./images/TOPCONTAMINANTS.pdf", height = 5, width = 10)
plot_taxa_heatmap(OTUfamcont, subset.top = 20, transformation = "clr",
                  taxonomic.level = "Family", border_color = "grey60",
                  VariableA = "Reindeer.ecotype")
dev.off()
