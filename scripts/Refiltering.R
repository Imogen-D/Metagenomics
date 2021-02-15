##Data Quality, Read Counts##

library(dplyr)
library(tidyverse)
library(phyloseq)
library(usethis)
use_git_config(user.name = "Imogen-D", user.email = "imogen.dumville@gmail.com")
library(reshape2)
library(dplyr)
library(phyloseq)
library(tidyr)
library(ggplot2)
library(microbiomeutilities)
library(decontam)
source("scripts/ancom_v2.1.R")
library(vegan)
library(devtools)
library(nlme)
library(tidyverse)
library(compositions)
library(pairwiseAdonis)


##Making object for ABUNDANCES -  BRACKEN output
full_otu <- read.delim("./data/kraken2_otu_table_merged_210203-otu.fungi.txt",na.strings = c("","NA"), stringsAsFactors=FALSE) %>% 
  select(which(colSums(.) > 0)) %>%   # remove empty taxa
  filter(rowSums(.) > 0) # remove empty samples

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
phydata <- phyloseq(OTU, sampledata, taxotable)

##Making phyloseq object for read information
reads <- read.delim("./data/kraken2_otu_table_merged_210203-reads.fungi.txt",na.strings = c("","NA"), stringsAsFactors=FALSE) %>% 
  select(which(colSums(.) > 0)) %>%   # remove empty taxa
  filter(rowSums(.) > 0) # remove empty samples

colnames(reads) <- str_replace(colnames(reads), "X", "")
readcounts <- otu_table(reads, taxa_are_rows = FALSE)
readphydata <- phyloseq(readcounts, sampledata, taxotable) #need to redo taxonomy table

#Only extraction # not with reindeer = BE103 and bear swabs BS003, BS005 - remove
#sample_names(phydata) != 
wanted <- !(sample_names(phydata) %in% c("BE103", "BS003", "BS005"))
phydata <- prune_samples(wanted, phydata)
wanted <- !(sample_names(readphydata) %in% c("BE103", "BS003", "BS005"))
readphydata <- prune_samples(wanted, readphydata)

#NOT PRUNING
#ABUNDANCE Fitlering

# abundance filtering function
abundance_filter_f <- function(count.mat, cutoff) {
  #count.mat must be in format of samples = rows and taxa = columns
  #row.names must be sample IDs
  #cutoff must be proportion (e.g. 0.001 = 0.1% relative abundance)
  count.filt <- count.mat
  prop.mat <- as.data.frame(prop.table(as.matrix(count.filt), margin = 1))
  prop.mat$Sample <- row.names(prop.mat)
  prop.mat.m <- melt(prop.mat, by="Sample", value.name = "Proportion", variable.name = "Taxon")
  samples <- as.character(prop.mat$Sample)
  for (s in samples) {
    exclude.taxa <- as.vector(subset(prop.mat.m, Sample==s & Proportion < cutoff)$Taxon)
    count.filt[s,exclude.taxa] = 0
  }
  return(count.filt)
}

phy.filt <- abundance_filter_f(as.data.frame(otu_table(phydata)), 0.0005) %>%
  select(which(colSums(.) > 0)) %>%   # remove empty taxa
  filter(rowSums(.) > 0)
phydata.filt <- phyloseq(otu_table(phy.filt, taxa_are_rows = FALSE), sampledata, taxotable)

phyread.filt <- abundance_filter_f(as.data.frame(otu_table(readphydata)), 0.0005) %>%
  select(which(colSums(.) > 0)) %>%   # remove empty taxa
  filter(rowSums(.) > 0)
phydata.filt.read <- phyloseq(otu_table(read.filt, taxa_are_rows = FALSE), sampledata, taxotable)

filteredout <- which(!sample_names(phydata) %in% sample_names(phydata.filt))
s <- sample_names(phydata)
s[filteredout] #"Rt022" "Rt066" "BE207" "BE210" "BkL2"  "BL105" "BL106" "BL107" "BL110" "BL112"

read.filteredout <- which(!sample_names(readphydata) %in% sample_names(phydata.filt.read))
sr <- sample_names(readphydata)
sr[read.filteredout] #"BE210" "BkL2"  "BL105" "BL106" "BL107" "BL110" "BL112" "Rt022" "Rt066"

##### CONTAMINANT FILTERING #####
# create single variable for control samples
controls<-phydata.filt@sam_data$Sample.R_cat %in% c("ExtBlank","LibBlank","Swab")
both.contaminants <- isContaminant(phydata.filt, method="either",neg = controls,conc = "Seq.copies.in.pool",threshold = c(0.1,0.5))
table(both.contaminants$contaminant)
phywocont <- prune_taxa(both.contaminants$contaminant==FALSE,phydata.filt)
phyWcont <- prune_taxa(both.contaminants$contaminant==TRUE,phydata.filt)


#decontam on reads data
readcontrols<-phydata.filt.read@sam_data$Sample.R_cat %in% c("ExtBlank","LibBlank","Swab")
both.contaminants.read <- isContaminant(phydata.filt.read, method="either",neg = readcontrols,conc = "Seq.copies.in.pool",threshold = c(0.1,0.5))
table(both.contaminants.read$contaminant)

decontaminantsreadphy <- prune_taxa(both.contaminants.read$contaminant==FALSE, phydata.filt.read)
readphyWcont <- prune_taxa(both.contaminants.read$contaminant==TRUE,phydata.filt.read)





meta_data <- data.frame(sample_data(phywocont)) %>% rownames_to_column("SampleID")
rt.samples <- meta_data %>% filter(grepl("^Rt",SampleID))
phywocont.rt <- prune_samples(rt.samples$SampleID,phywocont)
##remove taxa without any reads - not required so none are blank only
##which(colSums(otu.rt) == 0) ##0

##just extraction blanks - whats not there?
blanksamples <- meta_data %>% filter(grepl("Blank",Sample.R_cat))
phywocont.blanks <- prune_samples(blanksamples$SampleID,phywocont)

otu.blanks <- data.frame(otu_table(phywocont.blanks)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
colnames(otu.blanks) <- str_replace(colnames(otu.blanks), "X", "")
blank <- colnames(otu.blanks)
phyblanks <- prune_taxa(blank, phywocont.blanks)

##just swabs
swabsamples <- meta_data %>% filter(grepl("Swab",Sample.R_cat))
phywocont.swabs <- prune_samples(swabsamples$SampleID, phywocont)
otu.swabs <- data.frame(otu_table(phywocont.swabs)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
colnames(otu.swabs) <- str_replace(colnames(otu.swabs), "X", "")
swab <- colnames(otu.swabs)
physwabs <- prune_taxa(swab, phywocont.swabs)

##how many overlap between swabs and blanks?
overlapphy <- prune_taxa(taxa_names(physwabs), phyblanks) #383 taxa


##so now need to do all of above but with read data
##only reindeer
meta_data_read <- data.frame(sample_data(decontaminantsreadphy)) %>% rownames_to_column("SampleID")
rt.samples.read <- meta_data_read %>% filter(grepl("^Rt",SampleID))
phywocont.rt.read <- prune_samples(rt.samples.read$SampleID,decontaminantsreadphy)

otu.rt.read <- data.frame(otu_table(phywocont.rt.read)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
colnames(otu.rt.read) <- str_replace(colnames(otu.rt.read), "X", "")
rt.read <- colnames(otu.rt.read)
physamples.read <- prune_taxa(rt.read, phywocont.rt.read)

###10 taxa removed - rt samples only


##just extraction blanks - whats not there?
blanksamples.read <- meta_data_read %>% filter(grepl("Blank", Sample.R_cat))
phywocont.blanks.read <- prune_samples(blanksamples.read$SampleID,decontaminantsreadphy)

otu.blanks.read <- data.frame(otu_table(phywocont.blanks.read)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
colnames(otu.blanks.read) <- str_replace(colnames(otu.blanks.read), "X", "")
blank.read <- colnames(otu.blanks.read)
phyblanks.read <- prune_taxa(blank.read, phywocont.blanks.read)
#266 taxa not in blanks


##just swabs
swabsamples.read <- meta_data_read %>% filter(grepl("Swab",Sample.R_cat))
phywocont.swabs.read <- prune_samples(swabsamples.read$SampleID, decontaminantsreadphy)
otu.swabs.read <- data.frame(otu_table(phywocont.swabs.read)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
colnames(otu.swabs.read) <- str_replace(colnames(otu.swabs.read), "X", "")
swab.read <- colnames(otu.swabs.read)
physwabs.read <- prune_taxa(swab.read, phywocont.swabs.read)
##341 taxa not in swabs

##how many overlap between swabs and blanks?
overlapphy.read <- prune_taxa(taxa_names(physwabs.read), phyblanks.read) #288 taxa

##okay so number of taxa aren't the same but will trim dataframes
#facets for samples only, also in blanks, also in swabs, also in blanks + swabs = 4
samplesonly <- data.frame(otu_table(phywocont.rt))
samplesonly.read <- data.frame(otu_table(phywocont.rt.read))
inblanks <- data.frame(otu_table(phyblanks))
inblanks.read <- data.frame(otu_table(phyblanks.read))
inswabs <- data.frame(otu_table(physwabs))
inswabs.read <- data.frame(otu_table(physwabs.read))
inswabsandblanks <- data.frame(otu_table(overlapphy))
inswabsandblanks.read <- data.frame(otu_table(overlapphy.read))
colnames(samplesonly) <- str_replace(colnames(samplesonly), "X", "")
colnames(samplesonly.read) <- str_replace(colnames(samplesonly.read), "X", "")
colnames(inblanks) <- str_replace(colnames(inblanks), "X", "")
colnames(inblanks.read) <- str_replace(colnames(inblanks.read), "X", "")
colnames(inswabs) <- str_replace(colnames(inswabs), "X", "")
colnames(inswabs.read) <- str_replace(colnames(inswabs.read), "X", "")
colnames(inswabsandblanks) <- str_replace(colnames(inswabsandblanks), "X", "")
colnames(inswabsandblanks.read) <- str_replace(colnames(inswabsandblanks.read), "X", "")


#trim in blansk to only those not in colanmes(Swabs and blanks)
inblanks.read <- inblanks.read[,!(colnames(inblanks.read)) %in% colnames(inswabsandblanks.read)]

#same for swabs
inswabs.read <- inswabs.read[,!(colnames(inswabs.read)) %in% colnames(inswabsandblanks.read)]

#PLOTS
boxplot(colSums(samplesonly.read), colSums(inblanks.read), colSums(inswabsandblanks.read), colSums(inswabs.read), names = c("Samples", "Blanks", "S and B", "Swabs"))
boxplot(colSums(samplesonly.read), main = "Samples")
boxplot(colSums(inblanks.read), main = "Blanks")
boxplot(colSums(inswabsandblanks.read), main = "SwabsAndBlanks")
boxplot(colSums(inswabs.read), main = "Swabs")
barchart(colSums(samplesonly.read), colSums(inblanks.read), colSums(inswabsandblanks.read), colSums(inswabs.read), names = c("Samples", "Blanks", "S and B", "Swabs"))

#plotting contaminants
pdf(file = "./images/TOPCONTAMINANTS_1502.pdf", height = 5, width = 12)
plot_taxa_heatmap(phyWcont, subset.top = 14, transformation = "clr",
                  taxonomic.level = "Genus", border_color = "grey60",
                  VariableA = "Sample.R_cat")
dev.off()

pdf(file = "./images/TOPREADCONTAMINANTS_1502.pdf", height = 5, width = 12)
plot_taxa_heatmap(phyWcont, subset.top = 13, transformation = "clr",
                  taxonomic.level = "Genus", border_color = "grey60",
                  VariableA = "Sample.R_cat")
dev.off()

pdf(file = "./images/TOPTAXA_1502.pdf", height = 5, width = 12)
plot_taxa_heatmap(phywocont, subset.top = 20, transformation = "clr",
                  taxonomic.level = "Genus", border_color = "grey60",
                  VariableA = "Sample.R_cat")
dev.off()

pdf(file = "./images/TOPREADTAXA_1502.pdf", height = 5, width = 12)
plot_taxa_heatmap(decontaminantsreadphy, subset.top = 20, transformation = "clr",
                  taxonomic.level = "Genus", border_color = "grey60",
                  VariableA = "Sample.R_cat")
dev.off()

pdf(file = "./images/TOPTAXAinsamples_1502.pdf", height = 5, width = 12)
plot_taxa_heatmap(phywocont.rt, subset.top = 20, transformation = "clr",
                  taxonomic.level = "Genus", border_color = "grey60",
                  VariableA = "Sample.R_cat")
dev.off()


df <- data.frame(colSums(otu_table(phywocont.rt)))
top_samples <- df %>% slice_max(df, n = 20) #list of top 20 abundant taxa

##now can see where from, if cross contam in swabs and/or blanks?
dfblanks <- data.frame(colSums(otu_table(phywocont.blanks)))
top_blanks <- dfblanks %>% slice_max(dfblanks, n = 20)

dfswabs <- data.frame(colSums(otu_table(phywocont.swabs)))
top_swabs <- dfswabs %>% slice_max(dfswabs, n = 20)

rownames(top_samples) %in% rownames(top_blanks)
rownames(top_samples) %in% rownames(top_swabs)
