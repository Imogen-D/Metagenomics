#ANCOM thang


feature_table = otu_table(phywocont); sample_var = "Seq.label"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, sample_data(phywocont), sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

main_var = "Reindeer.ecotype"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start # around 30s

write_csv(res$out, "outputs/res_moving_pics.csv")