#ANCOM thang

View(otu_table(phywocont))
meta_data = as.data.frame(sample_data(phywocont))
feature_table = as.data.frame(otu_table(phywocont))


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
