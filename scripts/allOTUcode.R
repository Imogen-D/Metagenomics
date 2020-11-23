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
library(ggplot2)
library(microbiomeutilities)

setwd("~/MEME/Uppsala_Katja_Project/Metagenomics") #for local script

#set_entrez_key("ee1b29805250345f705302e643b1bfc4e007")
full_otu <- read.delim("./data/reindeer_kraken2_otu_table_merged_201112-otu.fungi.txt",na.strings = c("","NA"), row.names=1, stringsAsFactors=FALSE) %>% 
  replace(., is.na(.), 0) %>% 
  select(which(colSums(.) > 0)) # remove taxa summing to zero

colnames(full_otu) <- str_replace(colnames(full_otu), "X", "")

metadata <- read.delim("./data/reindeer_sample_metadate_merged.txt", stringsAsFactors=FALSE) %>% 
  distinct(Seq.label,.keep_all=T) %>%  # remove duplicates, the tidyr way
  filter(!is.na(Seq.label))
rownames(metadata)<-metadata$Seq.label # add rownames

OTU <- otu_table(full_otu, taxa_are_rows = FALSE)

#making TaxonomyTable
OTUtaxa <- classification(colnames(full_otu), db = "ncbi")

bound1<-bind_rows(as_tibble(cbind(OTUtaxa))) %>%
  select(kingdom,phylum,class,order,family,genus,species)
rownames(bound1)<-names(OTUtaxa)

write.csv(bound1, "./data/OTUtaxonomyformatted.csv")

OTUtaxonomyformatted <- read.csv("./data/OTUtaxonomyformatted.csv", row.names=1, stringsAsFactors=FALSE) %>% # read in taxa table saved from taxize 
  rename_all(str_to_title)  # make the column names into title case

taxotable <- tax_table(as.matrix(OTUtaxonomyformatted)) # matrix required for tax table
sampledata <- sample_data(metadata[sample_names(OTU),]) # only take the samples that are present in the OTU table

phydata <- phyloseq(OTU, sampledata,taxotable)

#subsetting for quicker analysis
noeco <- (which(is.na(metadata$Reindeer.ecotype)))
ecotypemeta <- sample_data(metadata[-c(noeco),]) #only those with ecotype data
ecophy <- phyloseq(OTU, ecotypemeta,taxotable)

topN <- 100
most_abundant_taxa <- sort(taxa_sums(ecophy), TRUE)[1:topN]
GP100 <- prune_taxa(names(most_abundant_taxa), ecophy)

plot_heatmap(GP100,sample.label = "Seq.label",taxa.label = "Family", max.label = ntaxa(GP100)) + facet_grid(~ Reindeer.ecotype,scales="free",space="free_x")


#3000 OTU
top3000 <- 3000
most_abundant_taxa_3000 <- sort(taxa_sums(ecophy), TRUE)[1:top3000]
GP3000 <- prune_taxa(names(most_abundant_taxa_3000), ecophy)
plot_heatmap(GP3000, sample.label = "Seq.label", taxa.label = "Family", max.label = ntaxa(GP3000)) + facet_grid(~ Reindeer.ecotype, scales="free",space="free_x")

#pdf(file = "./images/100eco.pdf", height = 50, width = 50)
#plot_heatmap(GP100, sample.label = "Reindeer.ecotype", taxa.label = "species") #SampleType could be metadata - bind 1 + 2 and other data
#dev.off()

pdf(file = "./images/heatmap100wide.pdf", height = 25, width = 50)

plot_taxa_heatmap(ecophy,
                  subset.top = 100, taxonomic.level = "Family", 
                  transformation = "log10", VariableA = "Reindeer.ecotype")
dev.off()

OTUfam <- microbiome::aggregate_taxa(ecophy, "Family")

pdf(file = "./images/heatmap1000family.pdf", height = 25, width = 50)
plot_taxa_heatmap(OTUfam, subset.top = 1000, transformation = "log10",
                  taxonomic.level = "Family",
                  VariableA = "Reindeer.ecotype")
dev.off()

pdf(file = "./images/heatmap20familyclr.pdf", height = 5, width = 10)
plot_taxa_heatmap(OTUfam, subset.top = 20, transformation = "clr",
                  taxonomic.level = "Family", border_color = "grey60",
                  VariableA = "Reindeer.ecotype")
dev.off()