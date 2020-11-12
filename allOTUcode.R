#Analysing OTU and metagenomic data for all reindeer samples; extracting fungal taxa
#phyloseq
#12112020 OTU without fungi so new OTU table


library(usethis)
use_git_config(user.name = "Imogen-D", user.email = "imogen.dumville@gmail.com")

library(phyloseq)
library(rentrez)
library(taxize)
library(stringr)

set_entrez_key("ee1b29805250345f705302e643b1bfc4e007")
full_otu <- read.delim("~/MEME/Uppsala_Katja_Project/Metagenomics/reindeer_kraken2_otu_table_merged_201112-otu.fungi.txt", row.names=1, stringsAsFactors=FALSE)
metadata <- read.delim("~/MEME/Uppsala_Katja_Project/Metagenomics/reindeer_sample_metadate_merged.txt", stringsAsFactors=FALSE)
nodupmeta <- metadata[-c(which(duplicated(metadata$Seq.label))),]
is.na(nodupmeta$Seq.label) #14
nodupmeta <- nodupmeta[-14,]
rownames(nodupmeta) <- nodupmeta$Seq.label


#maybe tidy up metadata?? 
OTU <- otu_table(full_otu, taxa_are_rows = TRUE)
sampledata <- sample_data(nodupmeta)
phylo
phydata <- phyloseq(OTU, sampledata)
plot_bar(phydata, x = "Sample.no")

#next = make OTU on x axis & colour by sample/metadata. Go through metadata for interests
colnames(nodupmeta)
