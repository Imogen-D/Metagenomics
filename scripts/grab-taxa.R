# get taxa from ncbi database
grab_taxa<-function(taxa.db,otu,taxa.file) {
  # arguments:
  # taxa.db = full path to a zipped archive of the NCBI taxa database. I used ncbitax2lin (https://github.com/zyxue/ncbitax2lin).
  # otu = OTU table with taxa as columns and samples as rows
  # name and path of the resulting subsetted taxa file

  require(tidyverse) # really only dplyr is needed

  #  if column names have an "X" in front, remove it
  #   taxa<-if(grepl("^X",colnames(otu))){
  #   gsub("X","",colnames(otu),ignore.case = F) 
  # }
  
  full_taxa<-read_csv(taxa.db,na = "NA",
                      col_types = cols(tax_id = col_character(),superkingdom = col_character(),phylum = col_character(),class = col_character(),order = col_character(),family = col_character(),genus = col_character(),species = col_character(),biotype = col_character(),clade = col_character(),clade1 = col_character(),
                                       clade10 = col_character(),clade11 = col_character(),clade12 = col_character(),clade13 = col_character(),clade14 = col_character(),clade15 = col_character(),clade16 = col_character(),clade17 = col_character(),clade18 = col_character(),clade19 = col_character(),clade2 = col_character(),
                                       clade3 = col_character(),clade4 = col_character(),clade5 = col_character(),clade6 = col_character(),clade7 = col_character(),clade8 = col_character(),clade9 = col_character(),cohort = col_character(),forma = col_character(),`forma specialis` = col_character(),`forma specialis1` = col_character(),
                                       genotype = col_character(),infraclass = col_character(),infraorder = col_character(),isolate = col_character(),kingdom = col_character(),morph = col_character(),`no rank` = col_character(),`no rank1` = col_character(),`no rank2` = col_character(),`no rank3` = col_character(),`no rank4` = col_character(),
                                       parvorder = col_character(),pathogroup = col_character(),section = col_character(),series = col_character(),serogroup = col_character(),serotype = col_character(),`species group` = col_character(),`species subgroup` = col_character(),strain = col_character(),subclass = col_character(),subcohort = col_character(),
                                       subfamily = col_character(),subgenus = col_character(),subkingdom = col_character(),suborder = col_character(),subphylum = col_character(),subsection = col_character(),subspecies = col_character(),subtribe = col_character(),subvariety = col_character(),superclass = col_character(),superfamily = col_character(),superorder = col_character(),
                                       superphylum = col_character(),tribe = col_character(),varietas = col_character())) %>%
    mutate(tax_id = paste0("X",as.character(tax_id))) %>% # add the 
    filter(tax_id %in% colnames(otu)) # filter taxa that match
    # create a subset df of only the taxa that are matching from the database
    write.csv(x = full_taxa,file = taxa.file,quote = F,row.names = F)  
}
