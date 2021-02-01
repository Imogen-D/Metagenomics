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


#copied from WithContaminants.R
#error atm: 
###Error in alr_data[, -(1:i), drop = FALSE] : 
###incorrect number of dimensions
outage <- ANCOM(feature_tableage, meta_dataage, struc_zero = struc_zeroage, main_var = "Spec.coll.year")
write.table(out$out, file = "./images/ANCOM/ecowcolour.txt")
pdf(file = "./images/ANCOM/ecowcolour.pdf", height = 5, width = 12)
out$fig
dev.off()
