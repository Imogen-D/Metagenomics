#Script for Plots (in combination with allOTUcode.R)

#pruning for top 100 taxa
topN <- 100
most_abundant_taxa <- sort(taxa_sums(ecophy), TRUE)[1:topN]
GP100 <- prune_taxa(names(most_abundant_taxa), ecophy)

plot_heatmap(GP100,sample.label = "Seq.label",taxa.label = "Family", max.label = ntaxa(GP100)) + facet_grid(~ Reindeer.ecotype,scales="free",space="free_x")

#pdf(file = "./images/100eco.pdf", height = 50, width = 50)
#plot_heatmap(GP100, sample.label = "Reindeer.ecotype", taxa.label = "species") #SampleType could be metadata - bind 1 + 2 and other data
#dev.off()

#3000 OTU
top3000 <- 3000
most_abundant_taxa_3000 <- sort(taxa_sums(ecophy), TRUE)[1:top3000]
GP3000 <- prune_taxa(names(most_abundant_taxa_3000), ecophy)
plot_heatmap(GP3000, sample.label = "Seq.label", taxa.label = "Family", max.label = ntaxa(GP3000)) + facet_grid(~ Reindeer.ecotype, scales="free",space="free_x")


pdf(file = "./images/heatmap100wide.pdf", height = 25, width = 50)
plot_taxa_heatmap(ecophy,
                  subset.top = 100, taxonomic.level = "Family", 
                  transformation = "log10", VariableA = "Reindeer.ecotype")
dev.off()



pdf(file = "./images/heatmap1000family.pdf", height = 25, width = 50)
plot_taxa_heatmap(OTUfam, subset.top = 1000, transformation = "log10",
                  taxonomic.level = "Family",
                  VariableA = "Reindeer.ecotype")
dev.off()