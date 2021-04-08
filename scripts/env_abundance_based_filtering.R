#load packages
library(tidyverse)
library(phyloseq)
library(funrar)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(dplyr)

abundance_filter_f <- function(count.mat, cutoff) {
  count.filt <- count.mat
  prop.mat <- as.data.frame(prop.table(as.matrix(count.filt), margin = 1))
  prop.mat$Sample <- row.names(prop.mat)
  prop.mat.m <- reshape2::melt(prop.mat, by="Sample", value.name = "Proportion", variable.name = "Taxon")
  samples <- as.character(prop.mat$Sample)
  for (s in samples) {
    exclude.taxa <- as.vector(subset(prop.mat.m, Sample==s & Proportion < cutoff)$Taxon)
    count.filt[s,exclude.taxa] = 0
  }
  return(count.filt)
}

#Community-level taxonomic analysis - Script 5
#Decontamination method based on the relative abundance 
#of taxa insamples and environmental controls
spe_data<-readRDS("DC-Fungi/Metagenomics/data/phyloseqwithoutcontaminants.rds")
spe_data_filt <- data.frame(t(otu_table(spe_data))) %>% 
  rownames_to_column("Tax_ID")

#Transform abundances into relative abundances (using otu table after 0.01% filtering)
species_table_filt_ra <- spe_data_filt %>% 
  select("Tax_ID",starts_with("Rt")) %>% 
  as.data.frame()

## Env controls
reindeer_swabs_filt <- spe_data_filt %>% as.data.frame() %>% dplyr::select("Tax_ID",starts_with("BS"))
  
#Filter for abundance using the same threshold

# add in env swabs from my own data
env_swabs<-spe_data_filt %>% select("Tax_ID","BS001","BS002") %>% 
  mutate(across(where(is.numeric),~./sum(.)))

#Add swab BS005 to the dataset (skull swab with the highest richness) after calculating relative abundance
species_table_filt_ra<-species_table_filt_ra %>% 
  left_join(env_swabs,by="Tax_ID") %>% 
  column_to_rownames("Tax_ID")

# NAs are generate where there was no match. Since this means that the species was absent, turn these to 0.
species_table_filt_ra<-apply(species_table_filt_ra,MARGIN = 2,FUN =  function(x) replace_na(x,0)) %>% 
  as.data.frame()

#### Calculate relative abundance ratios ####
#Define function for calculating relative abundance ratios
rel.abund.ratios <- function(taxa_table_ra, sample_type_vector) {
  #Calculate rel. abundance ratios between each sample and one environmental control
  samplenames <- colnames(taxa_table_ra)[which(sample_type_vector=="sample")]
  controls <- colnames(taxa_table_ra)[which(sample_type_vector=="control")]
  
  #Get long table with columns containing: taxon, sample, rel. abund in sample, rel. abund in control1, rel. abund in control 2 etc...
  taxa_table_ra_m <- reshape(taxa_table_ra, idvar = "taxon", ids=row.names(taxa_table_ra),
                             times=names(taxa_table_ra)[-c(which(sample_type_vector=="control"))],
                             timevar = "sample", varying = list(names(taxa_table_ra)[-c(which(sample_type_vector=="control"))]),
                             direction = "long")
  k <- ncol(taxa_table_ra_m)
  colnames(taxa_table_ra_m)[k-1] <- "rel.abund"
  #Reorder of columns
  taxa_table_ra_m <- taxa_table_ra_m[,c(k,k-2,k-1,1:(k-3))]
  
  #Make long again
  taxa_table_ra_m <- reshape(taxa_table_ra_m, idvar="taxon.sample", ids=rownames(taxa_table_ra_m),
                             times=names(taxa_table_ra_m)[4:ncol(taxa_table_ra_m)], timevar = "control",
                             varying = list(names(taxa_table_ra_m)[4:ncol(taxa_table_ra_m)]),
                             direction = "long")
  
  #Tidy up table a bit
  taxa_table_ra_m[,5] <- as.numeric(taxa_table_ra_m[,5])
  row.names(taxa_table_ra_m) <- NULL
  colnames(taxa_table_ra_m) <- c("taxon", "sample", "rel.abund", "control", "rel.abund.in.control")
  taxa_table_ra_m <- taxa_table_ra_m[,-6]
  
  
  #Remove rows where the control abundances is 0 (will produce infinite values when calculating the ratios)
  taxa_table_ra_m <- taxa_table_ra_m[which(taxa_table_ra_m$rel.abund.in.control>0),]
  
  #Create a matrix with the relative abundance ratio for each sample and control combination
  ra_ratio <- taxa_table_ra_m[,c("taxon","sample","control")]
  ra_ratio$ratio <- taxa_table_ra_m$rel.abund/taxa_table_ra_m$rel.abund.in.control
  return(ra_ratio)
}

#Create vector that indicates in a sample is control or not
control.or.not <- ifelse(grepl("BS", colnames(species_table_filt_ra)),"control", "sample")

#Calculate relative abundance ratios
ra_ratio <- rel.abund.ratios(species_table_filt_ra, control.or.not)

#### Scatterplots ####
#Plot separately for each sample
# ra_ratio_plots <- list()
# 
# samplenames <- colnames(species_table_filt_ra)[which(control.or.not=="sample")]
# 
# for (i in 1:length(samplenames)) {
#   #Get a subset per sample
#   subset <- ra_ratio[ra_ratio$sample==samplenames[i],]
#   #Get the abundance rank of every taxon in the sample
#   order <- species_table_filt_ra[order(species_table_filt_ra[,samplenames[i]], decreasing=TRUE),]
#   order <- cbind(rownames(order), 1:nrow(order))
#   #Add rank info to subset
#   subset$rank <- as.numeric(order[match(subset$taxon, order[,1]),2])
#   #Remove taxa with zero abundance in samples
#   subset <- subset[which(subset$ratio > 0),]
#   #Log transform
#   subset$ratio <- log(subset$ratio)
#   #Plot
#   plot <- ggscatter(subset, "rank", "ratio", color = "control") +
#     annotate("segment", size = 0.5, y = 0, yend = 0, x=0, xend=max(subset$rank), colour = "black") +
#     xlab("Abundance rank in sample") + ylab(colnames(subset)[i])
#   ra_ratio_plots[[i]] <- plot
# }
# 
# #Plot grid
# n <- length(ra_ratio_plots)
# nCol <- floor(sqrt(n))
# ra_plot_grid <- grid.arrange(grobs = ra_ratio_plots, ncol = nCol)
# 
# #Save plot
# ggsave("relabund_ratio_per_sample.png",
#   plot = ra_plot_grid,
#   device = "png",
#   path = "C://Users//MARKELLA//OneDrive - Uppsala universitet//Degree Project//Bioinformatics//T3_community-level",
#   width = 100, height = 100, units = "cm")

##### Get "environmental" taxa ####

#Create a table that shows if a taxon is less abundant in the sample (for every sample-control comparison)
less_in_samples <- ra_ratio[,c(1,2,3)]
less_in_samples$less.abundant.in.sample <- sapply(ra_ratio$ratio, function(x) {x<1})

#Create a table that shows is a taxon passes the criterion 
#of being less abundant in all samples than at least one control
#If at least one sample per group is FALSE (more abundant in the sample), then the whole group is FALSE
always_less <- less_in_samples %>% group_by(taxon, control) %>%
  dplyr::summarize(always.more.in.env = all(less.abundant.in.sample == "TRUE")) %>%
  #if at least one of the comparisons yield TRUE the taxon gets marked as TRUE
  group_by(taxon) %>% dplyr::summarise(environmental = any(always.more.in.env == "TRUE"))
table(always_less$environmental)
#270

#Get species names of the taxa identified as environmental
env_taxa <- always_less$taxon[which(always_less$environmental=="TRUE")]
env_taxa <- as.data.frame(env_taxa)
taxonomy_species<-read.csv("DC/data/taxa-names.csv")
env_taxa$species <- taxonomy_species[match(env_taxa$env_taxa, taxonomy_species$tax_id),"species"]
colnames(env_taxa) <- c("TaxID", "Species")

#Are any of these oral?

#### Lists of oral taxa
# James's list of hominid core oral microbiome taxa.
core_micr <- read.csv("DC/rank-abundance-filtering/data/james-taxa-wID.csv",header = T) %>%
  # filter(!is.na(tax_id) & tax_id %in% unique(diff.abund$Tax_ID)) %>% # only consider the taxa that are found in environmental controls
  mutate(tax_id=as.character(tax_id))

# HOMD
homd<-read.table("DC/data/homd_taxonomy_table.txt",header=T,sep="\t",fill=T) %>%
  # filter(!is.na(NCBI_taxon_id) & NCBI_taxon_id %in% unique(diff.abund$Tax_ID)) %>% # only consider the taxa that are found in environmental controls
  mutate(NCBI_taxon_id=as.character(NCBI_taxon_id)) %>%
  distinct(NCBI_taxon_id,.keep_all=TRUE)

#core hominid microbiome
env_taxa[which(env_taxa$Species %in% core_micr$Taxon),]
# 0

#HOMD
env_taxa[which(env_taxa$TaxID %in% homd$NCBI_taxon_id),]
# 0

#### Final dataset ####
#remove environmental controls
spe_data_final <- subset_samples(spe_data, !(grepl("BS", colnames(spe_data_filt)[-1])))

#keep only present taxa
spe_data_final <- prune_taxa(taxa_sums(spe_data_final)>0,spe_data_final)
spe_data_final <- prune_taxa(!taxa_names(spe_data_final) %in% env_taxa$TaxID,spe_data_final)

spe_data_final <- phyloseq(otu_table(spe_data_final, taxa_are_rows = TRUE),
                           tax_table(spe_data_final),
                           sample_data(spe_data_final))

saveRDS(spe_data_final,file = "DC-Fungi/Metagenomics/data/phyloseqwithoutENVcontaminants.rds")

#How many are oral?
#core hominid microbiome
sum(taxa_names(spe_data_final) %in% core_micr$tax_id)
#111

#HOMD
sum(taxa_names(spe_data_final) %in% homd$NCBI_taxon_id)
#196

#Extract a list of "exogenous" taxa -- based on this list, the taxa that have been removed in this analysis
#will also be removed in the read level using Kraken-tools
exogenous_id_list <- taxa_names(prune_taxa(!(taxa_names(spe_data) %in% taxa_names(spe_data_final)),spe_data))
write.table(exogenous_id_list, file="DC-Fungi/Metagenomics/data/env-exogenous-id-list.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)