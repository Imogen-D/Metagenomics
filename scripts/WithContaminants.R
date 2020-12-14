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

# prevalence based filtering
prev.contaminants <- isContaminant(phydata, method="prevalence",neg=controls, threshold = 0.5)
table(prev.contaminants$contaminant)

# frequency based filtering
freq.contaminants <- isContaminant(phydata, method="frequency",conc= "Seq.copies.in.pool", threshold = 0.1)
table(freq.contaminants$contaminant)

# both prevalance and frequency based filtering
both.contaminants <- isContaminant(phydata, method="either",neg = controls,conc = "Seq.copies.in.pool",threshold = c(0.1,0.5))
table(both.contaminants$contaminant)

# overlap between frequency and both
length(intersect(rownames(freq.contaminants[freq.contaminants$contaminant==TRUE,]),
          rownames(both.contaminants[both.contaminants$contaminant==TRUE,])))

# overlap between prevalance and both
length(intersect(rownames(prev.contaminants[prev.contaminants$contaminant==TRUE,]),
rownames(both.contaminants[both.contaminants$contaminant==TRUE,])))

# each of these methods exclude a non-overlapping set of contaminants (?)
length(intersect(intersect(rownames(prev.contaminants[prev.contaminants$contaminant==TRUE,]),
          rownames(freq.contaminants[freq.contaminants$contaminant==TRUE,])),
          rownames(both.contaminants[both.contaminants$contaminant==TRUE,])))

# use phyloseq to prune taxa instead
phywocont <- prune_taxa(both.contaminants$contaminant==FALSE,phydata)

# only pulled the reindeer samples for the graphic without contaminants, 
# need to do the same for the contaminant only graphic

# could maybe just subset out ecotype samples again, but then only 8??
# noeco <- (which(metadata$Reindeer.ecotype == "")) #39, some with weird names
# ecotypemeta <- sample_data(metadata[-c(noeco),])

keep_samples<-rownames(sample_data(phywocont)[grepl(sample_data(phywocont)$Seq.label, pattern = "Rt")])
ecophycont <- prune_samples(keep_samples,phywocont)

#aggregating at family level
OTUfamcont <- microbiome::aggregate_taxa(ecophycont, "Family")

#creating heatmap with clr transformation, top 20 families
pdf(file = "./images/heatmap20familycont.pdf", height = 5, width = 12)
plot_taxa_heatmap(OTUfamcont, subset.top = 20, transformation = "clr",
                  taxonomic.level = "Family", border_color = "grey60",
                  VariableA = "Reindeer.ecotype")
dev.off()

#looking only at contamination taxa
phyWcont <- prune_taxa(both.contaminants$contaminant==TRUE,phydata)

OTUfamcont <- microbiome::aggregate_taxa(phyWcont, "Family")

#creating heatmap with clr transformation, top 20 families OF CONTAMINATION
pdf(file = "./images/TOPCONTAMINANTS.pdf", height = 5, width = 12)
plot_taxa_heatmap(OTUfamcont, subset.top = 20, transformation = "clr",
                  taxonomic.level = "Family", border_color = "grey60",
                  VariableA = "Sample.R_cat")
dev.off()
