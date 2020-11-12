library(tidyverse)
meta <- read.csv("Metagenomics/reindeer_sample_metadate_merged.txt",sep="\t",header=T)

euk<-read.table("Metagenomics/kraken2_otu_table_merged_201112.txt",sep="\t",header=T,skip = 1,comment.char = "") %>% 
  rename(taxa="X.OTU.ID") %>%
  mutate(taxa=paste0("X",taxa)) %>% 
  pivot_longer(cols=-taxa) %>% 
  pivot_wider(names_from = "taxa") %>% 
  mutate(name=gsub("_kraken2_report","",name)) %>% 
  column_to_rownames("name")

source("DC2/scripts/grab-taxa.R")
grab_taxa(euk,"DC-Fungi/reindeer_kraken2_otu_table_merged_201112-taxa.csv")

taxa_names<-read.csv("DC-Fungi/reindeer_kraken2_otu_table_merged_201112-taxa.csv",header=T,na.strings = c("NA","")) %>% filter(kingdom=="Fungi")

taxa_names$genus.species<-ifelse(is.na(taxa_names$species),
                                 ifelse(is.na(taxa_names$genus),
                                        ifelse(is.na(taxa_names$family),as.character(taxa_names$order),as.character(taxa_names$family)),
                                        as.character(taxa_names$genus)),
                                 as.character(taxa_names$species))

rt.euk<-euk %>% 
  rownames_to_column("samples") %>% 
  # filter(grepl("^Rt",x = samples)) %>% # select only reindeer
  mutate(samples=gsub("\\_.*","",samples)) %>%  # clean up column names
  select(samples,which(colnames(.) %in% taxa_names$tax_id)) %>% 
  column_to_rownames("samples")

# colnames(rt.euk) <- ifelse(taxa_names$tax_id %in% colnames(rt.euk),taxa_names$species,taxa_names$tax_id)

write.table(rt.euk,file = "Metagenomics/reindeer_kraken2_otu_table_merged_201112-otu.fungi.txt",quote = F,row.names = T,col.names = T,sep="\t")
