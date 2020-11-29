library(tidyverse)
library(rentrez)
library(taxize)

meta <- read.csv("data/reindeer_sample_metadate_merged.txt",sep="\t",header=T)

euk<-read.table("data/kraken2_otu_table_merged_201129.txt",sep="\t",header=T,skip = 1,comment.char = "") %>% 
  rename(taxa="X.OTU.ID") %>%
  mutate(taxa=paste0("X",taxa)) %>% 
  pivot_longer(cols=-taxa) %>% 
  pivot_wider(names_from = "taxa") %>% 
  mutate(name=gsub("_kraken2_report","",name)) %>% 
  column_to_rownames("name")

# making TaxonomyTable
set_entrez_key("ee1b29805250345f705302e643b1bfc4e007")
OTUtaxa <- classification(gsub("X","",colnames(full_otu)), db = "ncbi")

# source("../../DC2/scripts/grab-taxa.R")
# grab_taxa(euk,"/home/adrian/Lab-Notes/DC-Fungi/Metagenomics/data/reindeer_kraken2_otu_table_merged_201129-taxa.csv")

bound1<-bind_rows(as_tibble(cbind(OTUtaxa)))
rownames(bound1)<-names(OTUtaxa)

write.csv(filter(bound1,kingdom=="Fungi"), "data/reindeer_kraken2_otu_table_merged_201129-fungi-taxa.csv")

taxa_names<-read.csv("data/reindeer_kraken2_otu_table_merged_201129-fungi-taxa.csv",header=T,na.strings = c("NA",""))

taxa_names$genus.species<-ifelse(is.na(taxa_names$species),
                                 ifelse(is.na(taxa_names$genus),
                                        ifelse(is.na(taxa_names$family),as.character(taxa_names$order),as.character(taxa_names$family)),
                                        as.character(taxa_names$genus)),
                                 as.character(taxa_names$species))

rt.euk<-euk %>% 
  rownames_to_column("samples") %>% 
  filter(!grepl("bracken",samples) & 
           !grepl("_201111",samples) & 
           !grepl("_201112",samples)) %>% 
  # filter(grepl("^Rt",x = samples)) %>% # select only reindeer
  mutate(samples=gsub("\\_.*","",samples),
         samples=gsub("_bracken","",samples)) %>%  # clean up column names
  select(samples,which(colnames(.) %in% as.character(paste0("X",taxa_names$X)))) %>% 
  column_to_rownames("samples")

# colnames(rt.euk) <- ifelse(taxa_names$tax_id %in% colnames(rt.euk),taxa_names$species,taxa_names$tax_id)

write.table(rt.euk,file = "data/reindeer_kraken2_otu_table_merged_201129-otu.fungi.txt",quote = F,row.names = T,col.names = T,sep="\t")
