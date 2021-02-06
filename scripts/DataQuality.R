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

decontaminantsreadphy <- prune_taxa(both.contaminants$contaminant==FALSE, prunedreadphy)
##blerg not same ntaxa()

##Pull out blanks = counts for pre decontam for balnks + samples
##Decontam? merge phyloseq objects
##Only reindeer blanks
##Plot TAXA not samples -> reads on y, abundances on x, colour by only in blank / blank and samples / sample only