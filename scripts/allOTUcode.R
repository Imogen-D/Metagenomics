#Analysing OTU and metagenomic data for all reindeer samples; extracting fungal taxa
#phyloseq
#12112020 OTU without fungi so new OTU table


library(usethis)
use_git_config(user.name = "Imogen-D", user.email = "imogen.dumville@gmail.com")
library(dplyr)
library(phyloseq)
library(rentrez)
library(taxize)
library(stringr)
library(tidyr)
library(data.table)

setwd("~/MEME/Uppsala_Katja_Project/Metagenomics") #for local script

set_entrez_key("ee1b29805250345f705302e643b1bfc4e007")
full_otu <- read.delim("./data/reindeer_kraken2_otu_table_merged_201112-otu.fungi.txt", row.names=1, stringsAsFactors=FALSE)
metadata <- read.delim("./data/reindeer_sample_metadate_merged.txt", stringsAsFactors=FALSE)
#nodupmeta <- metadata[-c(which(duplicated(metadata$Seq.label))),]
metadata <- metadata[-c(which(is.na(metadata$Seq.label))),] #removing unsequenced data
#is.na(nodupmeta$Seq.label) #14

#nodupmeta <- nodupmeta[-14,]


rownames(metadata) <- metadata$Seq.label

colnames(full_otu) <- str_replace(colnames(full_otu), "X", "")


#maybe tidy up metadata?? 
OTU <- otu_table(full_otu, taxa_are_rows = FALSE)

#making TaxonomyTable
#OTUtaxa <- classification(colnames(full_otu), db = "ncbi")
#bound1 <- bind_rows(OTUtaxa, .id = "column_label")
#write.csv(bound1, "./data/OTUtaxonomy.csv")

#OTUtaxonomy <- read.csv("~/MEME/Uppsala_Katja_Project/Metagenomics/data/OTUtaxonomy.csv", row.names=1, stringsAsFactors=FALSE)
#OTUtaxonomy <- OTUtaxonomy[,c(1:3)] #remove id column, otherwise messes with wide function
#wide <- pivot_wider(OTUtaxonomy, names_from = rank, values_from = name, id_cols = column_label)
#wide <- wide[,-c(2,4)]
#wide <- as.data.frame(wide, row.names = column_label)
#wide <- mutate_all(wide, as.character)
#rownames(wide) <- wide$column_label
#wide <- wide[,-1]
#write.csv(wide, "./data/OTUtaxonomyformatted.csv")

OTUtaxonomyformatted <- read.csv("~/MEME/Uppsala_Katja_Project/Metagenomics/data/OTUtaxonomyformatted.csv", row.names=1, stringsAsFactors=FALSE)
OTUtaxonomyformatted <- as.matrix(OTUtaxonomyformatted) #matrix required for tax table

taxotable <- tax_table(OTUtaxonomyformatted)
sampledata <- sample_data(metadata)

phydata <- phyloseq(OTU, sampledata, taxotable)

#subsetting for quicker analysis
noeco <- (which(is.na(metadata$Reindeer.ecotype)))
ecotypemeta <- sample_data(metadata[-c(noeco),]) #only those with ecotype data
ecophy <- phyloseq(OTU, ecotypemeta, taxotable)

topN <- 100
most_abundant_taxa <- sort(taxa_sums(ecophy), TRUE)[1:topN]
GP100 <- prune_taxa(names(most_abundant_taxa), ecophy)

plot_bar(GP100, "Reindeer.ecotype", fill="family") # facet_grid=~Reindeer.ecotype

#heatmap function not working yet
plot_heatmap(GP100, sample.label = "Reindeer.ecotype", taxa.label = "family") #SampleType could be metadata - bind 1 + 2 and other data


