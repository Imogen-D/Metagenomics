#script to produce heatmap for top 20 fungal families in reindeer samples WITH CONTAMINANTS REMOVED

library(usethis)
use_git_config(user.name = "Imogen-D", user.email = "imogen.dumville@gmail.com")
library(dplyr)
library(phyloseq)
library(rentrez)
library(taxize)
library(stringr)
library(tidyr)
library(data.table)
library(ggplot2)
library(microbiomeutilities)
library(decontam)
library(tibble)
source("scripts/ancom_v2.1.R")
library(vegan)
library(plyr)


setwd("~/MEME/Uppsala_Katja_Project/Metagenomics") #for local script

full_otu <- read.delim("./data/reindeer_kraken2_otu_table_merged_201129-otu.fungi.txt",na.strings = c("","NA"), stringsAsFactors=FALSE) %>% 
  select(which(colSums(.) > 0))  # remove empty taxa
  # filter(rowSums(.) > 0) # remove empty samples

colnames(full_otu) <- str_replace(colnames(full_otu), "X", "")
OTU <- otu_table(full_otu, taxa_are_rows = FALSE)

#reading metadata table
metadata <- read.delim("./data/Sample_processing_masterlist.txt", stringsAsFactors=FALSE) %>% 
  select(-starts_with("X")) %>% 
  distinct(Seq.label,.keep_all=T) %>%  # remove duplicates, the tidyr way
  filter(!is.na(Seq.label)) %>% 
  mutate(Reindeer.ecotype=ifelse(Sample.R_cat == "Reindeer_pre","pre-historic", # fill in some of the blank ecotypes
                                 ifelse(Sample.R_cat %in% c("ExtBlank","LibBlank","Swab"),"blank",
                                        ifelse(Sample.R_cat == "" & is.na(Sample.R_cat),"unknown",Reindeer.ecotype))))
rownames(metadata)<-metadata$Seq.label # add rownames

#reading and formatting taxonomy table, made with classification
OTUtaxonomyformatted <- read.csv("./data/OTUtaxonomyformattedwcont.csv", row.names=1, stringsAsFactors=FALSE) %>% # read in taxa table saved from taxize 
  rename_all(str_to_title)  # make the column names into title case
taxotable <- tax_table(as.matrix(OTUtaxonomyformatted))

sampledata <- sample_data(metadata[sample_names(OTU),]) # only take the samples that are present in the OTU table

# making full phyloseq data format
phydata <- phyloseq(OTU, sampledata,taxotable)

##### CONTAMINANT FILTERING #####
# create single variable for control samples
controls<-phydata@sam_data$Sample.R_cat %in% c("ExtBlank","LibBlank","Swab")

# prevalence based filtering
prev.contaminants <- isContaminant(phydata, method="prevalence",neg=controls, threshold = 0.5)
table(prev.contaminants$contaminant)

# frequency based filtering
freq.contaminants <- isContaminant(phydata, method="frequency",conc= "Seq.copies.in.pool", threshold = 0.1)
table(freq.contaminants$contaminant)

# both prevalance and frequency based filtering
both.contaminants <- isContaminant(phydata, method="either",neg = controls,conc = "Seq.copies.in.pool",threshold = c(0.1,0.5))
table(both.contaminants$contaminant)

# overlap between frequency and both
length(intersect(rownames(freq.contaminants[freq.contaminants$contaminant==TRUE,]),
          rownames(both.contaminants[both.contaminants$contaminant==TRUE,])))

# overlap between prevalance and both
length(intersect(rownames(prev.contaminants[prev.contaminants$contaminant==TRUE,]),
rownames(both.contaminants[both.contaminants$contaminant==TRUE,])))

# each of these methods exclude a non-overlapping set of contaminants (?)
length(intersect(intersect(rownames(prev.contaminants[prev.contaminants$contaminant==TRUE,]),
          rownames(freq.contaminants[freq.contaminants$contaminant==TRUE,])),
          rownames(both.contaminants[both.contaminants$contaminant==TRUE,])))

# use phyloseq to prune taxa instead
phywocont <- prune_taxa(both.contaminants$contaminant==FALSE,phydata)

# only pulled the reindeer samples for the graphic without contaminants, 
# need to do the same for the contaminant only graphic

# could maybe just subset out ecotype samples again, but then only 8??
# noeco <- (which(metadata$Reindeer.ecotype == "")) #39, some with weird names
# ecotypemeta <- sample_data(metadata[-c(noeco),])

keep_samples<-rownames(sample_data(phywocont)[grepl(sample_data(phywocont)$Seq.label, pattern = "Rt")])
ecophycont <- prune_samples(keep_samples,phywocont)

#aggregating at family level
OTUfamcont <- microbiome::aggregate_taxa(ecophycont, "Family")

#creating heatmap with clr transformation, top 20 families
pdf(file = "./images/heatmap20genuscont.pdf", height = 5, width = 12)
plot_taxa_heatmap(ecophycont, subset.top = 20, transformation = "clr",
                  taxonomic.level = "Genus", border_color = "grey60",
                  VariableA = "Reindeer.ecotype")
dev.off()

#looking only at contamination taxa
phyWcont <- prune_taxa(both.contaminants$contaminant==TRUE,phydata)

OTUfamcont <- microbiome::aggregate_taxa(phyWcont, "Family")

#creating heatmap with clr transformation, top 20 families OF CONTAMINATION

pdf(file = "./images/TOPCONTAMINANTS_genus.pdf", height = 5, width = 12)
plot_taxa_heatmap(phyWcont, subset.top = 20, transformation = "clr",
                  taxonomic.level = "Genus", border_color = "grey60",
                  VariableA = "Sample.R_cat")
dev.off()

### Gut fungi
# list of genera from here:  https://doi.org/10.1371/journal.pone.0151220
# consider adding taxa from Rumen Fungi book chapter: 10.1007/978-81-322-2401-3_7 = same genera
gut.fungi<-c("Neocallimastix","Anaeromyces","Caecomyces","Cyllamyces","Orpinomyces","Piromyces")

# check if any of these are in the contaminant taxa
data.frame(tax_table(phyWcont)) %>% filter(Genus %in% gut.fungi)

# good, they are not excluded, so how does the abundance of these taxa differ?
# filter by matching genus, then grab taxa id
gut.fungi.taxa<-data.frame(tax_table(phywocont)) %>% filter(Genus %in% gut.fungi) %>% rownames(.)
# prune the OTU table
phygut<-prune_taxa(gut.fungi.taxa,phywocont)
# plot a heatmap
plot_taxa_heatmap(phygut,subset.top = 20,taxanomic.level="Genus",VariableA = "Sample.R_cat",transformation = "clr")
# even though the Piromyces are highly abundant in samples & blanks, they were not excluded by decontam, which is good (?)


plot_taxa_heatmap(phygut,subset.top = 20,taxanomic.level="Genus",VariableA = "Spec.coll.year",transformation = "clr")
#plot_taxa_heatmap(phygut,subset.top = 20,taxanomic.level="Genus",VariableA = "Sample.coll.date",transformation = "clr")
plot_taxa_heatmap(phygut,subset.top = 20,taxanomic.level="Genus",VariableA = "Ext.date",transformation = "clr")
plot_taxa_heatmap(phygut,subset.top = 20,taxanomic.level="Genus",VariableA = "Spec.province",transformation = "clr")
plot_taxa_heatmap(phygut,subset.top = 20,taxanomic.level="Genus",VariableA = "Spec.district",transformation = "clr")
plot_taxa_heatmap(phygut,subset.top = 20,taxanomic.level="Genus",VariableA = "Spec.locality",transformation = "clr")

meta_data <- as.data.frame(sample_data(phygut)) #I haven't used the ANCOM filtering method for structural zeros
feature_table <- t(otu_table(phygut))

out <- ANCOM(feature_table, meta_data, main_var = "Sample.R_cat")
write.table(out$out, file = "./images/ANCOM/sampletype.txt")
pdf(file = "./images/ANCOM/sampletype.pdf", height = 5, width = 12)
out$fig
dev.off()

out <- ANCOM(feature_table, meta_data, main_var = "Reindeer.ecotype")
write.table(out$out, file = "./images/ANCOM/ecotype.txt")
pdf(file = "./images/ANCOM/ecotype.pdf", height = 5, width = 12)
out$fig
dev.off()

out <- ANCOM(feature_table, meta_data, main_var = "Ext.date")
write.table(out$out, file = "./images/ANCOM/extractiondate.txt")
pdf(file = "./images/ANCOM/extractiondate.pdf", height = 5, width = 12)
out$fig
dev.off()

out <- ANCOM(feature_table, meta_data, main_var = "Spec.province")
write.table(out$out, file = "./images/ANCOM/province.txt")
pdf(file = "./images/ANCOM/province.pdf", height = 5, width = 12)
out$fig
dev.off()

out <- ANCOM(feature_table, meta_data, main_var = "Spec.district")
write.table(out$out, file = "./images/ANCOM/district.txt")
pdf(file = "./images/ANCOM/district.pdf", height = 5, width = 12)
out$fig
dev.off()


#human oral fungi from https://doi.org/10.1080/21505594.2016.1252015
#combined with rumen fungi, above, no overlap?? Seems odd
humanoralfungi <- c("Agaricus", "Alternaria", "Aspergillus", "Aureobasidium", "Bipolaris",
  "Bullera", "Candida", "Cladosporium", "Coprinus", "Cryptococcus", "Curvularia",
  "Cyberlindnera", "Cystofilobasidium", "Cytospora", "Debaryomyces", "Didymella",
  "Dioszegia", "Epicoccum", "Erythrobasidium", "Exophiala", "Filobasidium", "Fusarium",
  "Glomus", "Hanseniaspora", "Irpex", "Kluyveromyces", "Lenzites", "Leptosphaerulina",
  "Malassezia", "Mrakia", "Naganishia", "Penicillium", "Phaeosphaeria", "Phoma", "Pichia",
  "Pisolithus", "Pyrenochaetopsis", "Ramularia", "Rhizocarpon", "Rhizopus", "Rhodosporidiobolus",
  "Rhodotorula", "Saccharomyces", "Sarcinomyces", "Scedosporium", "Sporobolomyces", "Talaromyces",
  "Taphrina", "Tausonia", "Teratosphaeria", "Torulaspora", "Trametes", "Trichoderma", "Trichosporon",
  "Wallemia", "Neocallimastix","Anaeromyces","Caecomyces","Cyllamyces","Orpinomyces","Piromyces")

# check if any of these are in the contaminant taxa
data.frame(tax_table(phyWcont)) %>% filter(Genus %in% humanoralfungi)
#22 excluded, but remember @ genus level, so some species of genus may still be incl


# filter by matching genus, then grab taxa id
humangut.fungi.taxa<-data.frame(tax_table(phywocont)) %>% filter(Genus %in% humanoralfungi) %>% rownames(.)
# prune the OTU table
humanphygut<-prune_taxa(humangut.fungi.taxa,phywocont)
# plot a heatmap
plot_taxa_heatmap(humanphygut,subset.top = 20,taxanomic.level="Genus",VariableA = "Sample.R_cat",transformation = "clr")
#Two leftmost taxa very abundant? Rt1 and 009
pdf(file = "./images/humanoral.pdf", height = 5, width = 12)
plot_taxa_heatmap(humanphygut,subset.top = 20,taxanomic.level="Genus",VariableA = "Sample.R_cat",transformation = "clr")
dev.off()
#many aspergillus

all_meta_data <- data.frame(sample_data(humanphygut)) #I haven't used the ANCOM filtering method for structural zeros
all_feature_table <- t(otu_table(humanphygut))

##this function here, how do I incoporate a new colour variable??
out <- ANCOM(all_feature_table, all_meta_data, main_var = "Reindeer.ecotype", color = "Sample.R_cat")
write.table(out$out, file = "./images/ANCOM/ecowcolour.txt")
pdf(file = "./images/ANCOM/ecowcolour.pdf", height = 5, width = 12)
out$fig
dev.off()

#this one didn't run, too much data??
#with all data, jsut without contaminants
nocontmetadata <- data.frame(sample_data(phywocont))
nocontfeaturetable <- t(otu_table(phywocont))
nocontaminantsout <- ANCOM(nocontfeaturetable, nocontmetadata, main_var = "Reindeer.ecotype", color = "Sample.R_cat")

write.table(nocontaminantsout$out, file = "./images/ANCOM/alltaxa.txt")
pdf(file = "./images/ANCOM/alltaxa.pdf", height = 5, width = 12)
nocontaminantsout$fig
dev.off()



#ordination


#junk from here down
complete.cases(all_meta_data$Sample.R_cat)
k <- complete.cases()



#PERMANOVA
permanova <- adonis(t(all_feature_table) ~ Sample.R_cat, data = all_meta_data, permutations=99, method = "bray")
#print(as.data.frame(permanova$aov.tab)["Sample.R_cat", "Pr(>F)"])


#ran but missing data - also adonis need distance matrix???

permanova <- adonis(t(otu) ~ group,
                    data = meta, permutations=99, method = "bray")

# P-value
#print(as.data.frame(permanova$aov.tab)["group", "Pr(>F)"])
