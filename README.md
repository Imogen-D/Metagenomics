# Metagenomics

**OVERALL GOALS**
Analysis of differences between:
  Historical vs prehistoric samples (and associated change in distribtuion, enviornment, populaiton size etc)
  Additional comparison of Svlabard samples (with assumption that historic and svalbard samples will 'diverge' from prehistoric)
  General comparisons between locations (with _caution_) = reindeer_location_summaries
  
Assumptions that Svalbard may have differences due to dietary and environmental composition


**Code present:**
  reindeer_merged_kraken2 is Kraken code (tbc by Adrian)
  grab_taxa.R gets taxa ID from NCBI
  ancom2.1 needed for ANCOM function
  
  ANCOM and PERMANOVA are 'scrap' scripts for these processes
  
  MakingRDS makes orignal, uncleaned data (phyloseq-otu-base objects)
  Which are then used in AbundanceAndReadPlots in decontam + abundance filtering and to produce heatmaps of location and top taxa
  Can also use refiltering for decontam/abundance filtering and then also does ordination

  env_abundance_based_filtering (and sample_filtering) produced by Adrian for environmental filtering -> produces object that shouldn't need decontam or abundance filtering
  
Previously checked variaitons between abundances and read counts (hence 'read' data objects)

**Data used:**
    Sample_processing_masterlist (for metadata)
    kraken2_otu_table_merged_210216-otu.fungi.txt (for OTU table)
    TAXONOMY TABLE FROM ADRIAN
    
_all images are replacable - have not been remade since cleaning_


Ideal next steps:
  Using RDS as phywocont object to continue with ordination, heatmaps and ancom. Samples pulled as enviornmental contaminants by Adrian need to be removed (env-endogenous-id-list)
  Observing other research to see specific taxa and/or functional prcoesses in taxa in similar species
      Suggested papers: https://www.nature.com/articles/nbt.4110 ; https://www.nature.com/articles/ismej2017108
      Include compariosons of those fungi found in other ruminants but also associate present taxa with soil presence
  Observing taxa in samples that are dominant or _not_ dominant and expected or _not_ expected to be present
  
Note: there are little previous studies on fungi in dental calculus in ruminants so using rumen studies and assuming overlap due to mastication
