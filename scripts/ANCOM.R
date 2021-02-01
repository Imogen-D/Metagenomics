#ANCOM thang
library(phyloseq)
source("scripts/ancom_v2.1.R")

phywocont <- readRDS(file = "phyloseqwithoutcontaminants.rds")

# subset to only reindeer
meta_data <- data.frame(sample_data(phywocont)) %>% rownames_to_column("SampleID")
rt.samples<-meta_data %>% filter(grepl("Reindeer",Sample.R_cat))
phywocont.rt<-prune_samples(rt.samples$SampleID,phywocont)

##### ANCOM STRUCTURAL ZEROS #####

prepro = feature_table_pre_process(feature_table = data.frame(t(otu_table(phywocont.rt)))[1:100,], 
                                   rt.samples,
                                   sample_var = "SampleID",
                                   group_var = "Sample.R_cat",
                                   out_cut = 0.05,
                                   zero_cut = 0.90,
                                   neg_lb = FALSE,
                                   lib_cut = 0)

feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info


#copied from WithContaminants.R
out <- ANCOM(feature_table, meta_data, main_var = "Sample.R_cat", struc_zero = struc_zero)
write.table(out$out, file = "./images/ANCOM/ecowcolour.txt")
pdf(file = "./images/ANCOM/ecowcolour.pdf", height = 5, width = 12)
out$fig + geom_hline(yintercept = quantile(out$out$W,probs = 0.7))
dev.off()

#adjustment for age (I have run this and appears same as above- from using same structural zero table??)

#outadj <- ANCOM(feature_table, meta_data, main_var = "Reindeer.ecotype", struc_zero = struc_zero, adj_formula = "Spec.coll.year", color = "Sample.R_cat")
#write.table(outadj$out, file = "./images/ANCOM/ecowageadj.txt")
#pdf(file = "./images/ANCOM/ecowageadj.pdf", height = 5, width = 12)
#outadj$fig
#dev.off()


#okay so now wanting to visually look at the taxa that are signifcant @0.7
#have gone back to phyloseq object for this, but should I?

#open txt as table
ecowcolour <- out$out
rownames(ecowcolour) <- ecowcolour$taxa_id
significant <- ecowcolour[which(ecowcolour$detected_0.7 == "TRUE"),]
sigotu <- otu_table(full_otu[,which(colnames(full_otu) %in% (significant$taxa_id))], taxa_are_rows = FALSE)
rt.samples<-meta_data %>% filter(grepl("Reindeer",Sample.R_cat))
physig <- merge_phyloseq(sigotu, tax_table(phywocont), sample_data(rt.samples))

plot_taxa_heatmap(physig, subset.top = 20, taxanomic.level="Genus", VariableA = "Reindeer.ecotype", transformation = "clr")

#will also do as barplot
#this is a bit nuts lol, need to subset some of these - how to do top taxa?
plot_bar(physig, x="Reindeer.ecotype", fill="Family")


##### ANCOM STRUCTURAL ZEROS - now for AGE #####

preproage = feature_table_pre_process(feature_table = data.frame(t(otu_table(phywocont.age))), 
                                   rt.samples,
                                   sample_var = "Seq.label",
                                   group_var = "Spec.coll.year",
                                   out_cut = 0.05,
                                   zero_cut = 0.90,
                                   neg_lb = FALSE,
                                   lib_cut = 0)

feature_tableage = preproage$feature_table # Preprocessed feature table
meta_dataage = preproage$meta_data # Preprocessed metadata
struc_zeroage = preproage$structure_zeros # Structural zero info


##Reducing number of dimensions: just Svalbard vs hsitroic vs pre?

##### ANCOM STRUCTURAL ZEROS #####
#Using phywocont.rt, but adding colomn

sample <- data.frame(sample_data(phywocont.rt))
View(sample$Reindeer.ecotype)

generaleco <- c("Svalbard", "Land", "Land", "Land", "Land", "Land",
                "Land", "Land", "Land", "Land", "Land", "Svalbard",
                "Land", "Land", "Land", "Land", "Land", "Svalbard",
                "pre-historic", "Land", "pre-historic", "pre-historic",
                "Land", "Svalbard", "pre-historic", "Land", "Land",
                "pre-historic", "Land", "pre-historic", "Land", "Land",
                "Svalbard", "Land", "Svalbard")
sample$generaleco <- generaleco
sample <- sample_data(sample)
phyloweco <- merge_phyloseq(phywocont.rt, sample)


preproeco = feature_table_pre_process(feature_table = data.frame(t(otu_table(phyloweco))), 
                                   sample,
                                   sample_var = "Seq.label",
                                   group_var = "generaleco",
                                   out_cut = 0.05,
                                   zero_cut = 0.90,
                                   neg_lb = FALSE,
                                   lib_cut = 0)

feature_tableeco = preproeco$feature_table # Preprocessed feature table
meta_dataeco = preproeco$meta_data # Preprocessed metadata
struc_zeroeco = preproeco$structure_zeros # Structural zero info

outeco <- ANCOM(feature_tableeco, meta_dataeco, main_var = "generaleco", struc_zero = struc_zeroeco)
write.table(outeco$out, file = "./images/ANCOM/generaleco.txt")
pdf(file = "./images/ANCOM/generaleco.pdf", height = 5, width = 12)
outeco$fig + geom_hline(yintercept = quantile(out$out$W,probs = 0.7))
dev.off()


#oridnation with smaller eco dataset
reindeerwblanks <- prune_samples(allsample$Ext.batch %in% sample$Ext.batch,phywocont)

reindeersample <- data.frame(sample_data(reindeerwblanks))
blankeco <- c(rep("Blank", 18), generaleco)
reindeersample$generaleco <- blankeco
reindeersample <- sample_data(reindeersample)
reindeerwblanks <- merge_phyloseq(reindeerwblanks, reindeersample)

#sum(rowSums(otu_table(reindeerwblanks))==0)
#non.zero <- prune_samples(names(which(sample_sums(reindeerwblanks)>0)), reindeerwblanks)
#completephy <- microbiome::transform(phywocont.rt, "compositional")
completephy <- microbiome::transform(reindeerwblanks, "compositional")

data.ord <- phyloseq::ordinate(completephy, method = "NMDS", distance = "bray") #incomplete dataset
p1 = plot_ordination(completephy, data.ord,color = "generaleco")+
  stat_ellipse()

#pariwise with more general ecotypes
y <- pairwise.adonis(x=otu_table(completephy), factors=reindeersample$generaleco)
write.csv(y, file = "./images/pairwiseecotype3.csv")
