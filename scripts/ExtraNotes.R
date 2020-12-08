#extra crap off of with contam for now
#subsetting for quicker analysis removing those without ecotype information
noeco <- (which(is.na(metadata$Reindeer.ecotype)))
ecotypemeta <- sample_data(metadata[-c(noeco),])
ecophy <- phyloseq(OTU, ecotypemeta,taxotable)

#aggregating at family level
OTUfam <- microbiome::aggregate_taxa(ecophy, "Family")

#creating heatmap with clr transformation, top 20 families
pdf(file = "./images/heatmap20familyclr.pdf", height = 5, width = 10)
plot_taxa_heatmap(OTUfam, subset.top = 20, transformation = "clr",
                  taxonomic.level = "Family", border_color = "grey60",
                  VariableA = "Reindeer.ecotype")
dev.off()


OTUfam <- microbiome::aggregate_taxa(ecophy, "Family")

#creating heatmap with clr transformation, top 20 families
pdf(file = "./images/heatmap20familyclr.pdf", height = 5, width = 10)
plot_taxa_heatmap(OTUfam, subset.top = 20, transformation = "clr",
                  taxonomic.level = "Family", border_color = "grey60",
                  VariableA = "Reindeer.ecotype")
dev.off()