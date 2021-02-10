##Data Quality, Read Counts##

library(dplyr)
library(tidyverse)
library(phyloseq)
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
library(devtools)
library(nlme)
library(tidyverse)
library(compositions)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)


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

ntaxa(phydata)
#6492
ntaxa(readphydata)
#6795
nrow(both.contaminants)
#5491

sum(!taxa_names(phydata) %in% taxa_names(readphydata))
#1001

#pruning to common taxa
prunedreadphy <- prune_taxa(taxa_names(phydata), readphydata)
ntaxa(prunedreadphy)
#[1] 5491
prunedphy <- prune_taxa(taxa_names(prunedreadphy), phydata)
ntaxa(prunedphy)
#[1] 5491

##fitlered to above 5% threshold
# transform to relative abundance
# only OTUs greater than 5% relative abundance are kept.
phydata.ra  <- transform_sample_counts(prunedphy, function(x) x / sum(x) )
### from here prune taxa - > remake phylo with absolute #'s not relative
phydata.filt <- filter_taxa(prune_taxa(taxa_names(phydata.ra), prunedphy), function(x) mean(x) >= 0.05, TRUE)

#filter read data
phydata.filt.read <- prune_taxa(taxa_names(phydata.filt), prunedreadphy)


##### CONTAMINANT FILTERING #####
# create single variable for control samples
controls<-phydata.filt@sam_data$Sample.R_cat %in% c("ExtBlank","LibBlank","Swab")
both.contaminants <- isContaminant(phydata.filt, method="either",neg = controls,conc = "Seq.copies.in.pool",threshold = c(0.1,0.5))
table(both.contaminants$contaminant)
phywocont <- prune_taxa(both.contaminants$contaminant==FALSE,phydata.filt)

#decontam using abundance data
decontaminantsreadphy <- prune_taxa(both.contaminants$contaminant==FALSE, phydata.filt.read)

#sample <- data.frame(sample_data(phywocont))
#readsample <- data.frame(sample_data(decontaminantsreadphy))

##only reindeer
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
#contains 3 less samples Rt017, 11 and 13 ah

##remove taxa without any reads - not required so none are blank only
##which(colSums(otu.rt) == 0) ## all with reads

##just extraction blanks - whats not there?
blanksamples.read <- meta_data_read %>% filter(grepl("Blank", Sample.R_cat))
phywocont.blanks.read <- prune_samples(blanksamples.read$SampleID,decontaminantsreadphy)

otu.blanks.read <- data.frame(otu_table(phywocont.blanks.read)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
colnames(otu.blanks.read) <- str_replace(colnames(otu.blanks.read), "X", "")
blank.read <- colnames(otu.blanks.read)
phyblanks.read <- prune_taxa(blank.read, phywocont.blanks.read)
#one extra blank

##just swabs
swabsamples.read <- meta_data_read %>% filter(grepl("Swab",Sample.R_cat))
phywocont.swabs.read <- prune_samples(swabsamples.read$SampleID, decontaminantsreadphy)
otu.swabs.read <- data.frame(otu_table(phywocont.swabs.read)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
colnames(otu.swabs.read) <- str_replace(colnames(otu.swabs.read), "X", "")
swab.read <- colnames(otu.swabs.read)
physwabs.read <- prune_taxa(swab.read, phywocont.swabs.read)

##how many overlap between swabs and blanks?
overlapphy.read <- prune_taxa(taxa_names(physwabs.read), phyblanks.read) #380 taxa

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
inblanks <- inblanks[,!(colnames(inblanks)) %in% colnames(inswabsandblanks)]
inblanks.read <- inblanks.read[,!(colnames(inblanks.read)) %in% colnames(inswabsandblanks.read)]

#same for swabs
inswabs <- inswabs[,!(colnames(inswabs)) %in% colnames(inswabsandblanks)]
inswabs.read <- inswabs.read[,!(colnames(inswabs.read)) %in% colnames(inswabsandblanks.read)]

###PLOTS
plot(colSums(samplesonly), colSums(samplesonly.read), xlab = "Abundance", ylab = "Read counts", main = "samples only")

inswabs <- inswabs[,colnames(inswabs) %in% colnames(inswabs.read)]
inswabs.read <- inswabs.read[,colnames(inswabs.read) %in% colnames(inswabs)]
plot(colSums(inswabs), colSums(inswabs.read), xlab = "Abundance", ylab = "Read counts", main = "swabs and samples") #x and y differ

inblanks <- inblanks[,colnames(inblanks) %in% colnames(inblanks.read)]
inblanks.read <- inblanks.read[,colnames(inblanks.read) %in% colnames(inblanks)]
plot(colSums(inblanks), colSums(inblanks.read), xlab = "Abundance", ylab = "Read counts", main = "blanks and samples") #x and y differ

inswabsandblanks <- inswabsandblanks[,colnames(inswabsandblanks) %in% colnames(inswabsandblanks.read)]
inswabsandblanks.read <- inswabsandblanks.read[,colnames(inswabsandblanks.read) %in% colnames(inswabsandblanks)]
plot(colSums(inswabsandblanks), colSums(inswabsandblanks.read), xlab = "Abundance", ylab = "Read counts", main = "swabs, blanks and samples") #x and y differ

boxplot(colSums(samplesonly.read), colSums(inblanks.read), colSums(inswabsandblanks.read), colSums(inswabs.read), names = c("Samples", "Blanks", "S and B", "Swabs"), log="y")
boxplot(colSums(samplesonly.read), main = "Samples")
boxplot(colSums(inblanks.read), main = "Blanks")
boxplot(colSums(inswabsandblanks.read), main = "SwabsAndBlanks")
boxplot(colSums(inswabs.read), main = "Swabs")

sum(rowSums(samplesonly.read) <=  5)
unique(colSums(samplesonly.read))



#attempting to make data into frame for ggplot()
a <- as.data.frame(colSums(inswabsandblanks))
b <- as.data.frame(colSums(inswabsandblanks.read))
##need to merge by OTU

SB <- data.frame(colSums(inswabsandblanks), colSums(inswabsandblanks.read), rep("S+B", 371), colnames(inswabsandblanks))
S <- data.frame(colSums(inswabs), colSums(inswabs.read), rep("S+B", 83), colnames(inswabs))
B <- data.frame(colSums(inblanks), colSums(inblanks.read), rep("S+B", 294), colnames(inblanks))
Samp <- data.frame(colSums(samplesonly), colSums(samplesonly.read), rep("Reindeer", 5305), colnames(samplesonly))

colnames(SB) <- c("Abundance", "Reads", "Present", "OTU")
rownames(SB) <- SB$OTU
colnames(S) <- c("Abundance", "Reads", "Present", "OTU")
rownames(S) <- S$OTU
colnames(B) <- c("Abundance", "Reads", "Present", "OTU")
rownames(B) <- B$OTU
colnames(Samp) <- c("Abundance", "Reads", "Present", "OTU")
rownames(Samp) <- Samp$OTU

a <- bind_rows(SB, S)
b <- bind_rows(a, B)
all <- bind_rows(b, Samp)
rownames(all) <- all$OTU #there are duplicate values ahhh -> something when setting rownames of the read and abundance taxa // not actually equal