#script for re-filtering/swab sorting

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
library(reshape2)

phywocont <- readRDS("./data/phyloseqwithoutENVcontaminants.rds")
sampledata <- sample_data(phywocont)
location <- read.csv("data/reindeer_location_summaries.csv", row.names=1)

## Just keep reindeer
rt.samples <- sample_names(phywocont)[which(grepl("^Rt",sample_names(phywocont)))]
phywocont.rt <- prune_samples(rt.samples,phywocont)
otu.rt <- data.frame(otu_table(phywocont.rt)) %>%
  select(which(colSums(.) > 0))  # remove empty taxa
rt <- colnames(otu.rt)
rt <- str_remove(rt, "X")
phyrt <- prune_taxa(rt, phywocont.rt) #1199 taxa, 35 samples

otunort <- data.frame(otu_table(phywocont.rt)) %>%
  select(which(colSums(.) == 0))
taxanotrt <- colnames(otunort) #5

##just extraction blanks - whats not there?
blanksamples <- sample_names(phywocont)[which(grepl(c("BE|BL|Bk"),sample_names(phywocont)))]
phywocont.blanks <- prune_samples(blanksamples,phywocont) #1205 taxa
# what type of taxa are removed? (maybe more important question to ask for blanks)


pdf(file = "./images/TOPTAXA1304.pdf", height = 5, width = 12)
plot_taxa_heatmap(phyrt, subset.top = 10, transformation = "clr",
                  taxonomic.level = "Genus", border_color = "grey60",
                  VariableA = "Sample.R_cat")
dev.off()
tax <- data.frame(tax_table(phyrt))
generainrt <- unique(tax$Genus) #530


##ordination to see if any patterns present
rt <- prune_samples(names(which(sample_sums(phyrt)>0)), phyrt) #3 that are empty - check this out

sample_data(rt) <- sample_data(location)
transformedrt <- microbiome::transform(rt, "compositional")

data.ord <- phyloseq::ordinate(transformedrt, method = "NMDS", distance = "bray") #incomplete dataset
p1 = plot_ordination(transformedrt, data.ord, color = "Spec.nation.metadata")+
  stat_ellipse()

p1
