#Analysing OTU and metagenomic data for all reindeer samples; extracting fungal taxa
#phyloseq
#12112020 OTU without fungi so new OTU table


library(usethis)
use_git_config(user.name = "Imogen-D", user.email = "imogen.dumville@gmail.com")
library(dplyr)
library(phyloseq)
library(rentrez)
library(taxize)
library(stringr)
library(tidyr)
library(data.table)


set_entrez_key("ee1b29805250345f705302e643b1bfc4e007")
full_otu <- read.delim("~/MEME/Uppsala_Katja_Project/Metagenomics/data/reindeer_kraken2_otu_table_merged_201112-otu.fungi.txt", row.names=1, stringsAsFactors=FALSE)
metadata <- read.delim("~/MEME/Uppsala_Katja_Project/Metagenomics/data/reindeer_sample_metadate_merged.txt", stringsAsFactors=FALSE)
nodupmeta <- metadata[-c(which(duplicated(metadata$Seq.label))),]
is.na(nodupmeta$Seq.label) #14
nodupmeta <- nodupmeta[-14,]
rownames(nodupmeta) <- nodupmeta$Seq.label

colnames(full_otu) <- str_replace(colnames(full_otu), "X", "")


#maybe tidy up metadata?? 
OTU <- otu_table(full_otu, taxa_are_rows = FALSE)

#attempting to make TaxonomyTable
OTUtaxa <- classification(colnames(full_otu), db = "ncbi")
bound1 <- bind_rows(OTUtaxa, .id = "column_label")

bound1 <- bound1[,-"id"]

wide <- pivot_wider(bound1, names_from = rank, values_from = name, id_cols = column_label)
wide <- wide[,-c(2,4)]
wide <- as.data.frame(wide, row.names = column_label) #so I do this to change from tibble to df
rownames(wide) <- wide$column_label
#rownames(wide) <- as.character(wide$column_label)
colnames(wide)
wide <- wide[,-1]
charwide[] <- lapply(wide, as.character)
#charwide <- data.frame(lapply(wide, as.character), stringsAsFactors=FALSE)

charwide <- as.matrix(charwide) #matrix required for tax table



#is.character(wide[1,1])
#type.convert(wide)

#rownames(wide)
#class(wide)

#colnames(wide)
taxotable <- tax_table(charwide)

View(taxotable)

sampledata <- sample_data(nodupmeta)
phylo
phydata <- phyloseq(OTU, sampledata)
plot_bar(phydata, fill="Genus")

#next = make OTU on x axis & colour by sample/metadata. Go through metadata for interests
colnames(nodupmeta)


