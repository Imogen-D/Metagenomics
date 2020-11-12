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
transposed <- transpose(bound)
bound <- bind_rows(OTUtaxa, .id = "rank")

#bind list of dataframes with different lengths
#unique rank names as colnames


#lst <- list(a = 1:3, b = 1:4, name = c("A", "X"))  #a list for the example

#n <- max(unlist(lapply(lst, length)))              #check maximum length

#lstnew <- lapply(OTUtaxa, function(x) {ans <- rep(NA, length=16); 
#ans[1:length(x)]<- x; 
#return(ans)})

#tax_table(OTUmatrix)

dt <- setDT(bound1)
dc2 <- dcast.data.table(dt, column_label~rank, value.var = 'name')
dc <-dcast(dt, column_label~rank, value.var='name') #aggregate('name')
melt <- "melt"(bound1, "column_label", "name")

reshape(dt, direction = "wide", idvar=c("name", "column_label"), timevar="rank")

View(dt)

sampledata <- sample_data(nodupmeta)
phylo
phydata <- phyloseq(OTU, sampledata)
plot_bar(phydata, fill="Genus")

#next = make OTU on x axis & colour by sample/metadata. Go through metadata for interests
colnames(nodupmeta)

#classification("71600", db = "ncbi")
#test <- ncbi_get_taxon_summary("71600")
