##Data Quality, Read Counts##

library(dplyr)
library(tidyverse)
library(phyloseq)


#reads <- read.delim("./data/kraken2_otu_table_merged_210203-reads.fungi.txt",na.strings = c("","NA"), stringsAsFactors=FALSE) %>% 
  select(which(colSums(.) > 0)) %>%   # remove empty taxa
  filter(rowSums(.) > 0) # remove empty samples

colnames(reads) <- str_replace(colnames(reads), "X", "")
readcounts <- otu_table(reads, taxa_are_rows = FALSE)
readphydata <- phyloseq(readcounts, sampledata, taxotable)

##fitlered to aove 5% threshold

sample_names(readphydata)
sample_names(phydata.filt)
sample_names(phydata)
sum(str_count(sample_names(readphydata), "Rt"))
sum(str_count(sample_names(phydata), "Rt"))
sum(str_count(sample_names(readphydata), "B"))
sum(str_count(sample_names(phydata), "B"))

ntaxa(phydata)
#6492
ntaxa(readphydata)
#6795
nrow(both.contaminants)
#6492


prunedreadphy <- prune_taxa(taxa_names(phydata), readphydata)

ntaxa(prunedreadphy)
#[1] 5491

prunedphy <- prune_taxa(taxa_names(prunedreadphy), phydata)
ntaxa(prunedphy)
#[1] 5491


# transform to relative abundance
phydata.ra  <- transform_sample_counts(prunedphy, function(x) x / sum(x) )

# only OTUs greater than 5% relative abundance are kept.
thresh<-as.numeric(quantile(mean(microbiome::abundances(phydata.ra),na.rm = T),probs = 0.05))
#phydata.filt <- filter_taxa(phydata.ra, function(x) sum(x) >= thresh, TRUE)


### from here prune taxa - > remake phylo with absolute #'s not relative

phydata.filt <- filter_taxa(prune_taxa(taxa_names(phydata.ra), phydata), function(x) sum(x) >= thresh, TRUE)


##### CONTAMINANT FILTERING #####
# create single variable for control samples
controls<-phydata.filt@sam_data$Sample.R_cat %in% c("ExtBlank","LibBlank","Swab")

# both prevalance and frequency based filtering
both.contaminants <- isContaminant(phydata.filt, method="either",neg = controls,conc = "Seq.copies.in.pool",threshold = c(0.1,0.5))
table(both.contaminants$contaminant)

# use phyloseq to prune taxa instead
phywocont <- prune_taxa(both.contaminants$contaminant==FALSE,phydata.filt)

decontaminantsreadphy <- prune_taxa(both.contaminants$contaminant==FALSE, prunedreadphy)

sample <- data.frame(sample_data(phywocont))
readsample <- data.frame(sample_data(decontaminantsreadphy))

##only reindeer
meta_data <- data.frame(sample_data(phywocont)) %>% rownames_to_column("SampleID")
rt.samples<-meta_data %>% filter(grepl("Reindeer",Sample.R_cat))
phywocont.rt <- prune_samples(rt.samples$SampleID,phywocont)
##remove taxa without any reads - not required so none are blank only
##which(colSums(otu.rt) == 0) ## all with reads

##just extraction blanks - whats not there?
meta_data$Sample.R_cat #54 onwards are blanks
blanksamples <- meta_data[54:96,]
phywocont.blanks <- prune_samples(blanksamples$SampleID,phywocont)

otu.blanks <- data.frame(otu_table(phywocont.blanks)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
colnames(otu.blanks) <- str_replace(colnames(otu.blanks), "X", "")
blank <- colnames(otu.blanks)
phyblanks <- prune_taxa(blank, phywocont.blanks)

##just swabs
swabsamples <- meta_data[97:98,]
phywocont.swabs <- prune_samples(swabsamples$SampleID, phywocont)
otu.swabs <- data.frame(otu_table(phywocont.swabs)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
colnames(otu.swabs) <- str_replace(colnames(otu.swabs), "X", "")
swab <- colnames(otu.swabs)
physwabs <- prune_taxa(swab, phywocont.swabs)

##how many overlap between swabs and blanks?
overlapphy <- prune_taxa(taxa_names(physwabs), phyblanks) #400 taxa


##so now need to do all of above but with read data
##only reindeer
meta_data_read <- data.frame(sample_data(decontaminantsreadphy)) %>% rownames_to_column("SampleID")
rt.samples.read <-meta_data_read %>% filter(grepl("Reindeer",Sample.R_cat))
phywocont.rt.read <- prune_samples(rt.samples.read$SampleID,decontaminantsreadphy)
##remove taxa without any reads - not required so none are blank only
##which(colSums(otu.rt) == 0) ## all with reads

##just extraction blanks - whats not there?
meta_data_read$Sample.R_cat #before 48 are blanks
blanksamples.read <- meta_data_read[1:43,]
phywocont.blanks.read <- prune_samples(blanksamples.read$SampleID,decontaminantsreadphy)

otu.blanks.read <- data.frame(otu_table(phywocont.blanks.read)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
colnames(otu.blanks.read) <- str_replace(colnames(otu.blanks.read), "X", "")
blank.read <- colnames(otu.blanks.read)
phyblanks.read <- prune_taxa(blank.read, phywocont.blanks.read)

##just swabs
swabsamples.read <- meta_data_read[44:45,]
phywocont.swabs.read <- prune_samples(swabsamples.read$SampleID, decontaminantsreadphy)
otu.swabs.read <- data.frame(otu_table(phywocont.swabs.read)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
colnames(otu.swabs.read) <- str_replace(colnames(otu.swabs.read), "X", "")
swab.read <- colnames(otu.swabs.read)
physwabs.read <- prune_taxa(swab.read, phywocont.swabs.read)

##how many overlap between swabs and blanks?
overlapphy.read <- prune_taxa(taxa_names(physwabs.read), phyblanks.read) #393 taxa

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


##blanks and reindeer ?? all but one is same extraction
rtsample <- data.frame(sample_data(phywocont.rt))
sample <- data.frame(sample_data(phywocont))
blanks <-sample$Ext.batch %in% rtsample$Ext.batch









##Pull out blanks = counts for pre decontam for balnks + samples
##Decontam? merge phyloseq objects
##Only reindeer blanks
##Plot TAXA not samples -> reads on y, abundances on x, colour by only in blank / blank and samples / sample only
