##### ordination ##### 
# how many empty rows?
meta_data <- data.frame(sample_data(phywocont))
rt.samples<-meta_data %>% filter(grepl("Reindeer",Sample.R_cat))
phywocont.rt<-prune_samples(rownames(rt.samples),phywocont)
phywocont.rt <- prune_taxa(names(which(taxa_sums(phywocont.rt) > 0)),phywocont.rt)
View(sample_data(phywocont))

sum(rowSums(otu_table(phywocont.rt))==0)
non.zero<- prune_samples(names(which(sample_sums(phywocont.rt)>0)),phywocont.rt)
completephy <- microbiome::transform(phywocont.rt, "compositional")

data.ord <- phyloseq::ordinate(completephy, method = "NMDS", distance = "bray") #incomplete dataset
p1 = plot_ordination(completephy, data.ord,color = "Reindeer.ecotype")+
  stat_ellipse()

p2 = plot_ordination(completephy, data.ord,color = "Spec.coll.year")+
  stat_ellipse()

#when only reindeer V shitty plots

##### PERMANOVA ##### 
completephy.bray<-phyloseq::distance(completephy,method="bray")
permanova <- adonis(completephy.bray ~ Sample.R_cat, data = data.frame(sample_data(completephy)), permutations=99)
#P = 0.06

#print(as.data.frame(permanova$aov.tab)["Sample.R_cat", "Pr(>F)"])

ecotypepermanova <- adonis(completephy.bray ~ Reindeer.ecotype, data = data.frame(sample_data(completephy)), permutations=99)
#P = 0.01


#on AGE / spec.coll.year
age.samples<- meta_data %>% drop_na(Spec.coll.year)
phywocont.age<-prune_samples(rownames(age.samples),phywocont)
completephyage <- microbiome::transform(phywocont.age, "compositional")
completephyage.bray<-phyloseq::distance(completephyage,method="bray")

agepermanova <- adonis(completephyage.bray ~ as.character(Spec.coll.year), data = data.frame(sample_data(completephy)), permutations=99)
#P = 0.01

sample <- data.frame(sample_data(completephy))
y <- pairwise.adonis(x=otu_table(completephy), factors=sample$Reindeer.ecotype)
write.csv(y, file = "./images/pairwiseecotype.csv")
#do I want to use the distance matrix? Or doing bray by default anyway?
  
#completephy.bray ~ Reindeer.ecotype, data = data.frame(sample_data(completephy)))
pairwise