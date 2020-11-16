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

setwd("~/MEME/Uppsala_Katja_Project/Metagenomics") #for local script

set_entrez_key("ee1b29805250345f705302e643b1bfc4e007")

full_otu <- read.delim("./data/reindeer_kraken2_otu_table_merged_201112-otu.fungi.txt", row.names=1, stringsAsFactors=FALSE) %>% 
  replace(., is.na(.), 0) %>% 
  select(which(colSums(.) > 0)) # remove taxa summing to zero

colnames(full_otu) <- str_replace(colnames(full_otu), "X", "")

metadata <- read.delim("./data/reindeer_sample_metadate_merged.txt", stringsAsFactors=FALSE) %>% 
  distinct(Seq.label,.keep_all=T) %>%  # remove duplicates, the tidyr way
  filter(!is.na(Seq.label))
rownames(metadata)<-metadata$Seq.label # add rownames

#maybe tidy up metadata?? 
OTU <- otu_table(full_otu, taxa_are_rows = FALSE)

#making TaxonomyTable
OTUtaxa <- classification(colnames(full_otu), db = "ncbi")
bound1<-do.call(dplyr::bind_rows,OTUtaxa)
write.csv(bound1, "./data/OTUtaxonomy.csv")

wide <- bound1 %>% distinct(name,rank,id,.keep_all = T) %>% 
  pivot_wider(id_cols = id, names_from = rank, values_from = name) %>% 
  mutate_all(as.character) %>% 
  mutate(genus.species=paste0(genus," ",species))
rownames(wide)<-wide$id # add rownames
  
write.csv(wide, "./data/OTUtaxonomyformatted.csv")

OTUtaxonomyformatted <- read.csv("./data/OTUtaxonomyformatted.csv", row.names=1, stringsAsFactors=FALSE)
OTUtaxonomyformatted <- as.matrix(OTUtaxonomyformatted) #matrix required for tax table

taxotable <- tax_table(OTUtaxonomyformatted)
sampledata <- sample_data(metadata[sample_names(OTU),]) # only take the samples that are present in the OTU table

phydata <- phyloseq(OTU, sampledata,taxotable)

#subsetting for quicker analysis
noeco <- (which(is.na(metadata$Reindeer.ecotype)))
ecotypemeta <- sample_data(metadata[-c(noeco),]) #only those with ecotype data
ecophy <- phyloseq(OTU, ecotypemeta,taxotable)

topN <- 100
most_abundant_taxa <- sort(taxa_sums(ecophy), TRUE)[1:topN]
GP100 <- prune_taxa(names(most_abundant_taxa), ecophy)

plot_bar(GP100,x = "Seq.label", fill="species")+theme(legend.position="none") #+facet_grid=~Reindeer.ecotype

#heatmap function not working yet
plot_heatmap(GP100, sample.label = "Seq.label", taxa.label = "species") #SampleType could be metadata - bind 1 + 2 and other data


