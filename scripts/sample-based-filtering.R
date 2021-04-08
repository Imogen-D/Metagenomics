library(tidyverse)
library(phyloseq)
library(UpSetR)

## Env controls
# Petrous bone samples: consider this sample as an additional environmental standard
petrous<-read.table("DC/rank-abundance-filtering/data/ERR2503700_m_kraken2_report_bracken_species.tsv",sep="\t",skip = 1,comment.char = "",header = TRUE) %>%
  mutate(X.OTU.ID=as.character(X.OTU.ID)) %>%
  rename(Tax_ID=X.OTU.ID,
         BSPetrous=ERR2503700_m_kraken2_report_bracken_species)

# gorilla skin swab
skin.swab<-read.csv("DC/reproducibility-tests/Gorilla-M/metadata_gorilla_analysis.csv",header=T) %>%
  filter(Seq.label == "ERR2868193") %>% pull(Seq.label)
skin <- read.delim("DC/data/kraken/markella_kraken_table.txt", sep = '\t', skip = 1) %>%
  rename_all(. %>% gsub("\\_.*","",.)) %>%  # clean up column names
  select(X.OTU.ID,all_of(skin.swab)) %>% 
  rename(Tax_ID=X.OTU.ID,
         BSSkin=ERR2868193) %>% 
  mutate(Tax_ID=as.character(Tax_ID))

file.list<-c("full-unfiltered-af-bracken-max-highoral-otu",
             "full-0.01-af-bracken-max-highoral-otu",
             "full-0.03-af-bracken-max-highoral-otu",
             "full-0.05-af-bracken-max-highoral-otu")

for (f in file.list) {
  ps<-filter_taxa(readRDS(paste0("DC/data/RDS/",f,".rds")), function(x) sum(x) > 0, TRUE)
  
  # split into a separate object for sample data
  sam<-data.frame(sample_data(ps))
  
  # make sure to exclude blanks mislabeled as Gorilla
  G.exclude<-sam %>%
    filter(Spec.host=="Gorilla" & grepl("^B",SampleID)) %>% # how did these get in here?
    pull(SampleID)
  
  sam.s<-sam %>%
    filter(!SampleID %in% G.exclude) %>% # include all Gorillas?
    filter(Spec.host2!="Reindeer.pre") %>% # remove non-museum specimens
    select(SampleID,Reads.n.P5,Ext.batch,Spec.host,Sample.R_cat) %>%
    bind_rows(list(SampleID="BSPetrous",Reads.n.P5=NA,Ext.batch="Petrous",Spec.host="control",Sample.R_cat="Swab"),
              list(SampleID="BSSkin",Reads.n.P5=NA,Ext.batch="Skin",Spec.host="control",Sample.R_cat="Swab")) # add in line for petrous environmental sample
  
  # create a data frame of the abundance data
  # generate abundance ranks
  otu<-data.frame(otu_table(ps)) %>%
    rownames_to_column("Tax_ID") %>%
    pivot_longer(cols = c(-Tax_ID),names_to="SampleID",values_to="Abundance") %>% # split by sample type
    left_join(sam.s) %>% # add in the sample data
    filter(Sample.R_cat == "Swab" | Spec.host %in% c("Gorilla","Reindeer","Bear")) %>%  # we only care about swab + sample comparison
    ungroup()
  
  # create a data frame of the abundance data just from environmental samples
  swabs<-otu %>%
    filter(grepl("^BS",SampleID)) %>%
    mutate(Ext.batch=as.character(Ext.batch)) %>%
    rename(swabs.abundance=Abundance) %>%  # rename in prep for joining
    select(Tax_ID,SampleID,Ext.batch,swabs.abundance) %>%
    pivot_wider(id_cols = Tax_ID,names_from = "Ext.batch",values_from = "swabs.abundance",values_fill=NA) %>% # if a certain taxa is not present within one environmental sample, make sure to put an NA
    rename(reindeer.shelf.swab=`18`,reindeer.skull.swab=`19`,bear.shelf.swab=`20`,bear.skull.swab=`21`) %>% 
    full_join(petrous) %>% # add in the petrous dataframe
    full_join(skin) %>% # add in skin dataframe
    rename(Petrous=BSPetrous,
           Skin=BSSkin)

  ##### Differential abundance between samples and environmental controls
  diff.abund<-otu %>% # add abundance data in for the samples we are comparing
    ungroup() %>% 
    rename(samples=Abundance) %>%  # rename in prep for joining
    left_join(swabs,by="Tax_ID") %>% # join with abundance data of environmental controls
    filter(!is.na(samples)) %>% # want to keep the taxa that are not found in samples, but found in controls?
    mutate(across(c(reindeer.shelf.swab,reindeer.skull.swab,bear.shelf.swab,bear.skull.swab,Petrous,Skin),~ samples-.)) %>%  # create individual differential abundances
    pivot_longer(cols = c(ends_with("swab"),"Petrous","Skin"),names_to="env",values_to="abundance.diff") %>%
    filter(!is.na(Spec.host) & 
             !is.na(env) & 
             # !is.na(abundance.diff) &
             !grepl("^BS",SampleID)) # we don't need to compare the swabs back to themselves here
  
  #### Lists of oral taxa
  # James's list of hominid core oral microbiome taxa.
  hominid.df <- read.csv("DC/rank-abundance-filtering/data/james-taxa-wID.csv",header = T) %>%
    filter(!is.na(tax_id) & tax_id %in% diff.abund$Tax_ID) %>% # only consider the taxa that are found in environmental controls
    mutate(tax_id=as.character(tax_id))
  
  # HOMD
  homd.df<-read.table("DC/data/homd_taxonomy_table.txt",header=T,sep="\t",fill=T) %>%
    filter(!is.na(NCBI_taxon_id) & NCBI_taxon_id %in% diff.abund$Tax_ID) %>% # only consider the taxa that are found in environmental controls
    mutate(NCBI_taxon_id=as.character(NCBI_taxon_id)) %>%
    distinct(NCBI_taxon_id,.keep_all=TRUE)
  
  # p1<-diff.abund %>% 
  #   group_by(env,Spec.host,SampleID) %>%
  #   summarise(n_above=n_distinct(Tax_ID[abundance.diff>0])/n_distinct(Tax_ID),
  #             n_below=n_distinct(Tax_ID[abundance.diff<0])/n_distinct(Tax_ID)) %>%
  #   pivot_longer(cols = c(-Spec.host,-SampleID,-env)) %>%
  #   ggplot(aes(x=reorder(SampleID, value),y=value,fill=name))+
  #   geom_bar(stat="identity",color="black")+
  #   facet_grid(env~Spec.host,scales = "free",space = "free_x")+
  #   theme_classic()+
  #   theme(axis.text.x = element_text(angle = 90,hjust = 1,size=6))
  # ggsave(filename = paste0("DC/rank-abundance-filtering/images/",f,"-samples-taxa-above-below.png"),
  #        dpi=300,width = 20,units = "in",plot = p1)
  
  ##### Relative abundance between samples and environmental controls
  rel.abund<-otu %>%
    rename(samples=Abundance) %>%  # rename in prep for joining
    group_by(Spec.host,SampleID) %>% 
    filter(!is.na(samples)) %>% # want to keep the taxa that are not found in samples, but found in controls?
    left_join(swabs,by="Tax_ID") %>% # join with abundance data of environmental controls
    mutate(across(c(reindeer.shelf.swab,reindeer.skull.swab,bear.shelf.swab,bear.skull.swab,Petrous),~ samples/.)) %>%  # create individual differential abundances
    pivot_longer(cols=c(reindeer.shelf.swab,reindeer.skull.swab,bear.shelf.swab,bear.skull.swab,Petrous,Skin),names_to = "env",values_to = "abundance.rel") %>%
    filter(!is.na(Spec.host) & !is.na(env) & 
             !is.na(abundance.rel) &
             !grepl("^BS",SampleID)) # we don't need to compare the swabs back to themselves here
  
  ## this plot should be constructed based on all env samples?
  # p2<-rel.abund %>% group_by(Spec.host,SampleID,env) %>%
  #   arrange(Spec.host,SampleID,env,-samples) %>%
  #   mutate(abundance_rank=row_number()) %>%
  #   filter(env=="Petrous") %>%
  #   ggplot(aes(x=abundance_rank,y=abundance.rel,color=Spec.host))+
  #   geom_point(alpha=0.5)+
  #   facet_wrap(~SampleID)+
  #   geom_hline(yintercept = 1,color="black")+
  #   labs(x="Abundance Rank",y="Log10 Relative Abundance")+
  #   scale_y_log10()+
  #   theme_classic()+
  #   theme(axis.text.x = element_blank())
  # ggsave(filename = paste0("DC/rank-abundance-filtering/images/",f,"-samples-ranked-rel-abundance.png"),
  #        dpi=300,width = 20,height=20,units = "in",plot = p2)
  
  #### Taxa that are always more abundant in all environmental samples
  all.diffs<-diff.abund %>% ungroup() %>% 
    pivot_wider(c(-env,-abundance.diff),names_from=env,values_from=abundance.diff) %>% 
    group_by(SampleID,Spec.host,Tax_ID) %>%
    summarise(museum.env.negative=ifelse(!is.na(reindeer.shelf.swab) & reindeer.shelf.swab < 0 | 
                                           !is.na(reindeer.skull.swab) & reindeer.skull.swab < 0 |
                                           !is.na(bear.shelf.swab) & bear.shelf.swab < 0 |
                                           !is.na(bear.skull.swab) & bear.skull.swab < 0, 1, 0),
              petrous.negative=ifelse(!is.na(Petrous) & Petrous < 0,1,0),
              skin.negative=ifelse(!is.na(Skin) & Skin < 0,1,0)) %>% ungroup()
  
  # p3<-all.diffs %>%
  #   mutate(museum.env.negative.n=museum.env.negative,
  #          petrous.negative.n=museum.env.negative+petrous.negative,
  #          skin.negative.n=museum.env.negative+petrous.negative+skin.negative) %>% 
  #   pivot_longer(cols = ends_with("negative.n"),names_to="env") %>%
  #   group_by(Spec.host,env) %>%
  #   mutate(HOMD_match=ifelse(Tax_ID%in%homd.df$NCBI_taxon_id,1,0),
  #          hominid_match=ifelse(Tax_ID%in%hominid.df$tax_id,1,0)) %>%
  #   summarise(n=sum(value),
  #             n_homd=sum(HOMD_match[value==1]),
  #             n_hominid=sum(hominid_match[value==1])) %>%
  #   mutate(env=factor(x = env,levels = c("museum.env.negative.n", "petrous.negative.n", "skin.negative.n"))) %>% 
  #   ggplot(aes(group=Spec.host,x=env,y=n))+
  #   geom_point()+
  #   geom_line()+
  #   geom_bar(aes(y=n_homd,fill=Spec.host),color="black",stat = "identity")+
  #   geom_bar(aes(y=n_hominid,fill=Spec.host),color="black",stat = "identity")+
  #   facet_grid(.~Spec.host)+theme_bw()+theme(axis.text.x = element_text(angle=90))
  # ggsave(filename = paste0("DC/rank-abundance-filtering/images/",f,"-env-number-breakdown.png"),
  #        dpi=300,width = 20,height=20,units = "in",plot = p3)
  
  # save as a full list
  all.diffs %>% 
    mutate(HOMD_match=ifelse(Tax_ID%in%homd.df$NCBI_taxon_id,1,0),
           hominid_match=ifelse(Tax_ID%in%hominid.df$tax_id,1,0)) %>%      
    write.csv(quote=FALSE,row.names=FALSE,file=paste0("DC/rank-abundance-filtering/data/",f,"-taxa-all-diff-env.csv"))
  
  all.diffs.summary<-all.diffs %>%
    group_by(Spec.host) %>% 
    summarise(museum.env.negative.n=n_distinct(Tax_ID[museum.env.negative==1]),
              petrous.negative.n=n_distinct(Tax_ID[museum.env.negative==1|petrous.negative==1]),
              skin.negative.n=n_distinct(Tax_ID[museum.env.negative==1|petrous.negative==1|skin.negative==1]))
  all.diffs %>% ungroup() %>% 
    summarise(museum.env.negative.n=n_distinct(Tax_ID[museum.env.negative==1]),
              petrous.negative.n=n_distinct(Tax_ID[museum.env.negative==1|petrous.negative==1]),
              skin.negative.n=n_distinct(Tax_ID[museum.env.negative==1|petrous.negative==1|skin.negative==1])) %>% 
    mutate(Spec.host="Combined") %>% 
    bind_rows(all.diffs.summary) %>% 
    write.csv(quote=FALSE,row.names=FALSE,file=paste0("DC/rank-abundance-filtering/data/",f,"-taxa-all-diff-env-summary.csv"))
  
  HOMD.diffs.summary<-all.diffs %>%
    mutate(HOMD_match=ifelse(Tax_ID%in%homd.df$NCBI_taxon_id,1,0)) %>% 
    filter(HOMD_match==1) %>%
    group_by(Spec.host) %>% 
    summarise(starting_HOMD=n_distinct(homd.df$NCBI_taxon_id),
              museum.env.negative.n=n_distinct(Tax_ID[museum.env.negative==1]),
              petrous.negative.n=n_distinct(Tax_ID[museum.env.negative==1|petrous.negative==1]),
              skin.negative.n=n_distinct(Tax_ID[museum.env.negative==1|petrous.negative==1|skin.negative==1]))
  all.diffs %>% ungroup() %>% mutate(HOMD_match=ifelse(Tax_ID%in%homd.df$NCBI_taxon_id,1,0)) %>% 
    filter(HOMD_match==1) %>%
    summarise(starting_HOMD=n_distinct(homd.df$NCBI_taxon_id),
              museum.env.negative.n=n_distinct(Tax_ID[museum.env.negative==1]),
              petrous.negative.n=n_distinct(Tax_ID[museum.env.negative==1|petrous.negative==1]),
              skin.negative.n=n_distinct(Tax_ID[museum.env.negative==1|petrous.negative==1|skin.negative==1])) %>%  
    mutate(Spec.host="Combined") %>% 
    bind_rows(HOMD.diffs.summary) %>% 
    write.csv(quote=FALSE,row.names=FALSE,file=paste0("DC/rank-abundance-filtering/data/",f,"-HOMD-taxa-all-diff-env-summary.csv"))
  
  hominid.diffs.summary<-all.diffs %>%
    mutate(hominid_match=ifelse(Tax_ID%in%hominid.df$tax_id,1,0)) %>% 
    filter(hominid_match==1) %>%
    group_by(Spec.host) %>% 
    summarise(starting_hominid=n_distinct(hominid.df$tax_id),
              museum.env.negative.n=n_distinct(Tax_ID[museum.env.negative==1]),
              petrous.negative.n=n_distinct(Tax_ID[museum.env.negative==1|petrous.negative==1]),
              skin.negative.n=n_distinct(Tax_ID[museum.env.negative==1|petrous.negative==1|skin.negative==1]))
  all.diffs %>% ungroup() %>% mutate(hominid_match=ifelse(Tax_ID%in%hominid.df$tax_id,1,0)) %>% 
    filter(hominid_match==1) %>%
    summarise(starting_hominid=n_distinct(hominid.df$tax_id),
              museum.env.negative.n=n_distinct(Tax_ID[museum.env.negative==1]),
              petrous.negative.n=n_distinct(Tax_ID[museum.env.negative==1|petrous.negative==1]),
              skin.negative.n=n_distinct(Tax_ID[museum.env.negative==1|petrous.negative==1|skin.negative==1])) %>%  
    mutate(Spec.host="Combined") %>% 
    bind_rows(hominid.diffs.summary) %>% 
    write.csv(quote=FALSE,row.names=FALSE,file=paste0("DC/rank-abundance-filtering/data/",f,"-hominid-taxa-all-diff-env-summary.csv"))  
  #### further filtering?
  # decontam?
  # abundance filtering?
  
  # save upset plot to show intersections
  ## at least one env...
  # num.taxa.upset<- all.diffs %>% 
  # filter(!is.na(positive.one)) %>%
  # pivot_wider(id_cols = Tax_ID,values_from = positive.one,values_fill = 0,names_from = Spec.host) %>% as.data.frame()
  
  # png(paste0("DC/rank-abundance-filtering/images/",f,"-taxa-diff-one-positive-upset.png"),res = 300,height = 1400,width = 1400,units = "px")
  # print(upset(num.taxa.upset,nintersects = NA,order.by = "freq",empty.intersections = T))
  # dev.off()
  
  # num.taxa.upset<- all.diffs %>% 
  #   filter(!is.na(negative.one)) %>% 
  #   pivot_wider(id_cols = Tax_ID,values_from = negative.one,values_fill = 0,names_from = Spec.host) %>% as.data.frame()
  # 
  # png(paste0("DC/rank-abundance-filtering/images/",f,"-taxa-diff-one-negative-upset.png"),res = 300,height = 1400,width = 1400,units = "px")
  # print(upset(num.taxa.upset,nintersects = NA,order.by = "freq",empty.intersections = T))
  # dev.off()
  # 
  # # always more abundant
  # num.taxa.upset<- taxa.always.diff %>% 
  #   filter(!is.na(museum.env.negative)) %>% group_by(Spec.host,Tax_ID) %>% summarise(museum.env.negative=max(museum.env.negative,na.rm=T)) %>% 
  #   pivot_wider(id_cols = Tax_ID,values_from = museum.env.negative,values_fill = 0,names_from = Spec.host) %>% as.data.frame()
  # 
  # png(paste0("DC/rank-abundance-filtering/images/",f,"-taxa-always-negative-upset.png"),res = 300,height = 1400,width = 1400,units = "px")
  # print(upset(num.taxa.upset,nintersects = NA,order.by = "freq",empty.intersections = T))
  # dev.off()
  # 
  # num.taxa.upset<- taxa.always.diff %>% 
  #   filter(!is.na(museum.env.positive)) %>% group_by(Spec.host,Tax_ID) %>% summarise(museum.env.positive=max(museum.env.positive,na.rm=T)) %>% 
  #   pivot_wider(id_cols = Tax_ID,values_from = museum.env.positive,values_fill = 0,names_from = Spec.host) %>% as.data.frame()
  # 
  # png(paste0("DC/rank-abundance-filtering/images/",f,"-taxa-always-positive-upset.png"),res = 300,height = 1400,width = 1400,units = "px")
  # print(upset(num.taxa.upset,nintersects = NA,order.by = "freq",empty.intersections = T))
  # dev.off()
}
