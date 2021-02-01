#Script to produce taxonomy table (in combination with allOTUcode.R)

rentrez::set_entrez_key("ee1b29805250345f705302e643b1bfc4e007")

#making TaxonomyTable
OTUtaxa <- taxize::classification(colnames(full_otu), db = "ncbi")

bound1<-bind_rows(as_tibble(cbind(OTUtaxa))) %>%
  select(kingdom,phylum,class,order,family,genus,species)
rownames(bound1)<-names(OTUtaxa)

write.csv(bound1, "./data/OTUtaxonomyformatted.csv")