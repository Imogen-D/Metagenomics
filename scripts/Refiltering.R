#script for re-filtering/swab sorting

phydata <- readRDS("./data/phyloseq-otu-base.rds")
readphydata <- readRDS("./data/phyloseq-read-base.rds")
sampledata <- sample_data(phydata)
read_sampledata <- sample_data(readphydata)
location <- read.csv("data/reindeer_location_summaries.csv", row.names=1)

#### contaminant filtering workflow
# 1. Abundance (base) filtering
# 2. Swab prevalence (taxa found just in swabs not in blanks)
# 3. decontam with library and extraction blanks

#### ABUNDANCE Fitlering
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

## run abundance filtering of taxa on otu data
phy.filt <- abundance_filter_f(as.data.frame(otu_table(phydata)), 0.0005) %>%
  select(which(colSums(.) > 0))

phydata.filt <- phyloseq(otu_table(phy.filt, taxa_are_rows = FALSE), sampledata, tax_table(phydata))

#how many filtered out and how many still in blanks
filteredout <- which(!taxa_names(phydata) %in% taxa_names(phydata.filt))
length(filteredout)

# pruning, keeping only blanks
phyabundance.blanks <- prune_samples(grepl(c("BE|BL|Bk"),sample_names(phydata.filt)),phydata.filt)

# creating otu table, trimming empty taxa
abundanceblanks <- data.frame(otu_table(phyabundance.blanks)) %>%
  select(which(colSums(.) > 0))

# create a phyloseq of the taxa found in blanks
phyabundanceblanks <- prune_taxa(colnames(abundanceblanks), phyabundance.blanks) #282 taxa

## filter out taxa that are just in swabs and samples (i.e. not in blanks)
phy.swabs <- prune_samples(grepl("^BS",sample_names(phydata.filt)), phydata.filt)
otu.swabs <- data.frame(otu_table(phy.swabs)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
physwabs <- prune_taxa(colnames(otu.swabs), phy.swabs) #348 taxa

## how many overlap between swabs and blanks?
# take the swab taxa that  are not found in blanks
swab.excl<- !taxa_names(phy.swabs) %in% taxa_names(phyabundanceblanks)
phy.onlyswabs <- prune_taxa(swab.excl, phy.swabs) #282 taxa

# only keep the taxa that aren't found in swabs
# also remove swabs from phyloseq 
phy.woswab.excl <- !taxa_names(phydata.filt) %in% taxa_names(phy.onlyswabs)
phywoswabs <- prune_taxa(phy.woswab.excl, prune_samples(!grepl("^BS",sample_names(phydata.filt)),phydata.filt))

# how many taxa removed?
print(ntaxa(phydata.filt) - ntaxa(phywoswabs))

# what type of taxa are removed? (maybe more important question to ask for blanks)

##### CONTAMINANT FILTERING #####
### decontam on OTU
### taxa will be filtered out if identified by *either* of the following methods:
### frequency - Contaminants are identified by frequency that varies inversely with sample DNA concentration.
### prevalence - Contaminants are identified by increased prevalence in negative controls.

# create single variable for control samples
controls<-phywoswabs@sam_data$Sample.R_cat %in% c("ExtBlank","LibBlank")

# filter contaminants that are in blanks or are significantly associated with the amount of DNA input for sequencing
both.contaminants <- isContaminant(phywoswabs, method="either",neg = controls,conc = "Seq.copies.in.pool",threshold = c(0.1,0.5),normalize = T)

# how many contaminants were filtered out?
table(both.contaminants$contaminant) #none rip

# keep just the non-contaminant taxa
phywocont <- prune_taxa(both.contaminants$contaminant==FALSE,phywoswabs)
# keep just the contaminant taxa
##EMPTY phyWcont <- prune_taxa(both.contaminants$contaminant==TRUE,phywoswabs)

## Just keep reindeer
rt.samples <- sample_names(phywocont)[which(grepl("^Rt",sample_names(phywocont)))]
phywocont.rt <- prune_samples(rt.samples,phywocont)
otu.rt <- data.frame(otu_table(phywocont.rt)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
rt <- colnames(otu.rt)
phyrt <- prune_taxa(rt, phywocont.rt) #280 taxa, 56 samples

##remove taxa without any reads - 2 blank only????

##just extraction blanks - whats not there?
blanksamples <- sample_names(phywocont)[which(grepl(c("BE|BL|Bk"),sample_names(phywocont)))]
phywocont.blanks <- prune_samples(blanksamples,phywocont) #282 txa

#otu.blanks <- data.frame(otu_table(phywocont.blanks)) %>%
#  select(which(colSums(.) > 0))  # remove empty taxa
#blank <- colnames(otu.blanks)
#phyblanks <- prune_taxa(blank, phywocont.blanks) #still 282 taxa here, none are empty!!

##just swabs
swabsamples <- sample_names(phywocont)[which(grepl("^BS",sample_names(phywocont)))]
phywocont.swabs <- prune_samples(swabsamples, phywocont)
otu.swabs <- data.frame(otu_table(phywocont.swabs)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
swab <- colnames(otu.swabs)
physwabs <- prune_taxa(swab, phywocont.swabs) #54 are empty i.e. found in samples and not swabs, 228 total

