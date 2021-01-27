#ANCOM thang
library(phyloseq)

phywocont <- readRDS(file = "phyloseqwithoutcontaminants.rds")

##### ANCOM STRUCTURAL ZEROS #####
meta_data <- data.frame(sample_data(phywocont))

#trouble is here - as soon as I add the group_var it says "Error in feature_table_pre_process(feature_table, meta_data, sample_var,  : 
#argument "neg_lb" is missing, with no default
#In addition: Warning message:
#  In e[!is.na(group)] <- residuals(f_fit) :
#  number of items to replace is not a multiple of replacement length"
  
feature_table = t(otu_table(phywocont)); sample_var = "Seq.label"; group_var = "Reindeer.ecotype";
out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, out_cut, zero_cut, neg_lb)

feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

View(struc_zero) #without group_var is NULL



#copied from WithContaminants.R
out <- ANCOM(all_feature_table, all_meta_data, main_var = "Reindeer.ecotype", color = "Sample.R_cat")
write.table(out$out, file = "./images/ANCOM/ecowcolour.txt")
pdf(file = "./images/ANCOM/ecowcolour.pdf", height = 5, width = 12)
out$fig
dev.off()

#this one didn't run, too much data??
#with all data, just without contaminants
nocontmetadata <- data.frame(sample_data(phywocont))
nocontfeaturetable <- t(otu_table(phywocont))
nocontaminantsout <- ANCOM(nocontfeaturetable, nocontmetadata, main_var = "Reindeer.ecotype", color = "Sample.R_cat")

write.table(nocontaminantsout$out, file = "./images/ANCOM/alltaxa.txt")
pdf(file = "./images/ANCOM/alltaxa.pdf", height = 5, width = 12)
nocontaminantsout$fig
dev.off()


##OLD_SCRIPT##
#even dropping to no threshold it removes everything? I tried with the normal values as on https://github.com/FrederickHuangLin/ANCOM/blob/master/README.md
feature_table = otu_table(phywocont); sample_var = "Seq.label"; group_var = NULL;
out_cut = 0.00; zero_cut = 1; neg_lb = FALSE; lib_cut = 0
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, neg_lb)
#feature_table = prepro$feature_table # Preprocessed feature table
#meta_data = prepro$meta_data # Preprocessed metadata
#struc_zero = prepro$structure_zeros # Structural zero info

prepro$feature_table

main_var = "Reindeer.ecotype"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL
feature_table = otu_table(phywocont)
feature_table <- t(as.data.frame(feature_table))
meta_data = sample_data(phywocont)
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, main_var, struc_zero = NULL, p_adj_method, alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start # around 30s


#with adjustment for age of sample? (check spec.coll.year before running)
#haven't run this yet, was going to try do it on UPPMAX but struggling to set it up
main_var = "Reindeer.ecotype"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "Spec.coll.year"; rand_formula = NULL
feature_table = otu_table(phywocont)
feature_table <- t(as.data.frame(feature_table))
meta_data = sample_data(phywocont)
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, main_var, struc_zero = NULL, p_adj_method, alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start 

write_csv(res$out, "outputs/res_moving_pics.csv")
