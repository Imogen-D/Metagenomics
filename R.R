#Analysing OTU and metagenomic data for all reindeer samples; extracting fungal taxa
#phyloseq


library(usethis)
use_git_config(user.name = "Imogen-D", user.email = "imogen.dumville@gmail.com")

library(phyloseq)
library(rentrez)
library(taxize)

set_entrez_key("ee1b29805250345f705302e643b1bfc4e007")
DC1.fungal.otu <- read.delim("~/MEME/Uppsala_Katja_Project/Metagenomics/DC1-fungal-otu.txt")
Sample_processing_masterlist <- read.csv("~/MEME/Uppsala_Katja_Project/Metagenomics/Sample_processing_masterlist.csv", stringsAsFactors=FALSE)
sample <- c("Rt11", "Rt13", "Rt1", "Rt5", "Rt7") #manually from table
rownames(DC1.fungal.otu) <- sample
row.names(samplereindeer) <- samplereindeer$Seq.label
which(Sample_processing_masterlist$Seq.label == c("Rt11", "Rt13", "Rt1", "Rt5", "Rt7")) #41, 38, 24, 2, 15
samplereindeer <- Sample_processing_masterlist[c(41, 38, 24, 2, 15),]
samplereindeer <- samplereindeer[ , colSums(is.na(samplereindeer)) == 0]
taxa <- colnames(DC1.fungal.otu)

all_kingdoms <- data.frame()

for (x in taxa) {
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

         