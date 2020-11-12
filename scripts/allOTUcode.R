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
full_otu <- read.delim("~/MEME/Uppsala_Katja_Project/Metagenomics/reindeer_kraken2_otu_table_merged_201106.txt", row.names=1, stringsAsFactors=FALSE)
metadata <- read.csv("~/MEME/Uppsala_Katja_Project/Metagenomics/Sample_processing_masterlist.csv", stringsAsFactors=FALSE)
colnames(full_otu) <- str_remove(colnames(full_otu), "_kraken2_report")

reindeermeta <- metadata[which(str_detect(metadata$Seq.label, "Rt")),]

reindeermeta[38,5] <- "Rt052" #previously Vl001/Rt052

#maybe tidy up metadata?? 

classes <- data.frame()

classes <- classification(rownames(full_otu), db="ncbi")

Fungi <- which(str_detect(classes, pattern = "Fungi"))



for (x in rownames(full_otu)) {
  class <- classification(x, db = "ncbi")
  classes <- rbind.data.frame(classes, class)
}

###old code from other script

OTUtaxa <- rownames(full_otu)

all_kingdoms <- data.frame()


for (x in OTUtaxa) {
  y <- tax_name(x, "kingdom", db = "ncbi")
  Sys.sleep(0.1)  #ugly and slow, could be shorter but required as API only allows certain no of search per seconds
  all_kingdoms <- rbind.data.frame(all_kingdoms, y)
}

funginums <- which(all_kingdoms$kingdom == "Fungi")
fungisamples <- all_kingdoms[funginums,]

fungaltaxa <- DC1.fungal.otu[,fungisamples$query]

#currently have OTU table of 5 rt samples and table of OTU counts of fungi
#now: analyse

otu_table <- otu_table(fungaltaxa, taxa_are_rows = FALSE)
sampledata <- sample_data(samplereindeer)
phydata <- phyloseq(otu_table, sampledata)
plot_bar(phydata, fill="Sample.weight.g")

#next = make OTU on x axis & colour by sample/metadata. Go through metadata for interests

