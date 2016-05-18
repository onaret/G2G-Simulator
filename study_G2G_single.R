#rm(list = ls())
cat("\014")  
#library(rjson)
library(dplyr)
library(ggplot2)
library(tidyr)
library(parallel)

source("G2G_full.R")
source("G2G_single.R")
source("GWAS.R")
trace <- TRUE

####Different stratification scenarios
if(FALSE) {
	repnb = 1000
	study_design_5000 = get_study_design(5000,2,2)
	study_design_10000 = get_study_design(10000,2,2)
	
	########C1
	
	###CORRECTED C1
	c1 =  get_G2G_setup(1000, s_stratified = c("P1","P2"), s_biased = c("A","B"), s_partial_bias = "P1", y_stratified = c("A","B"), y_biased = c("P1","P2"), y_partial_strat = "P1")
	res = test_G2G_setup(study_design, c1, fst_pop_strat = 0.2, fst_pop_bias = 0.02, fst_strain_strat = 0.02, fst_strain_bias = 0.2)
	#######
	
	###TOFIX : Limite case, will make error
	c1 =  get_G2G_setup(1000, s_stratified = c("P1","P2"), s_biased = c("A","B"), s_partial_bias = "P1", y_stratified = c("A","B"), y_partial_strat = "P1")
	##########
	c1 =  get_G2G_setup(100, s_biased = c("A","B"), s_partial_bias = "P1", y_stratified = c("A","B"), y_partial_strat = "P1")
	res = test_G2G_setup(study_design_5000, c1, fst_pop_bias = 0.03, fst_strain_strat = 0.2)
	
	c1 =  get_G2G_setup(10, s_stratified = c("P1"," P2"), s_biased = c("A","B"), s_partial_bias = "P1", y_stratified = c("A","B"), y_partial_strat = "P1", y_biased = c("P1","P2"))
	res = test_G2G_setup(study_design_5000, c1, fst_pop_strat = 0.03, fst_pop_bias = 0.03, fst_strain_strat = 0.2, fst_strain_bias = 0.2)
	
	c1_eq = get_G2G_setup(repnb,s_biased = c("A","B"), s_partial_bias = "P1", y_stratified = c("A","B"), y_partial_strat = "P1")
	res = test_G2G_setup(study_design_5000, c1_eq, fst_pop_bias = 0.2, fst_strain_strat = 0.2)
	
	c1_without_human_bias = get_G2G_setup(repnb, s_stratified = c("P1"," P2"), y_stratified = c("A","B"), y_partial_strat = "P1")
	res = test_G2G_setup(study_design_5000, c1_without_human_bias, fst_pop_strat = 0.2, fst_strain_strat = 0.2)
	
	c1_with_full_bias_full_strat = get_G2G_setup(repnb, s_stratified =  c("P1"," P2"), s_biased = c("A","B"), y_stratified = c("A","B"), y_biased =  c("P1"," P2"))
	res = test_G2G_setup(study_design_5000, c1_with_full_bias_full_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.03, fst_pop_bias = 0.03)
	
	c1_with_full_bias_full_strat_eq = get_G2G_setup(repnb, s_stratified =  c("P1"," P2"), s_biased =c("A","B"), y_stratified = c("A","B"), y_biased =  c("P1"," P2"))
	res = test_G2G_setup(study_design_5000, c1_with_full_bias_full_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.2, fst_pop_bias = 0.2)
	
	
	########FS1
	fs1 = get_G2G_setup(repnb, s_stratified = c("P1", "P2"), y_biased = c("P1", "P2"))
	res = test_G2G_setup(study_design_5000, fs1, fst_pop_strat = 0.2, fst_strain_bias = 0.02)
	
	fs1_full_strat = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full")
	res = test_G2G_setup(study_design_5000, fs1_full_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.2)
	
	######fs2_ss1 
	fs2_ss1 = get_G2G_setup(repnb, s_stratified = "full", associated_strains = "full", associated_populations = "full",beta = 0.1)
	res = test_G2G_setup(study_design_5000, fs2_ss1, fst_pop_strat = 0.2, tag = "beta_0.1")
	fs2_ss1 = get_G2G_setup(repnb, s_stratified = "full", associated_strains = "full", associated_populations = "full",beta = 0.5)
	res = test_G2G_setup(study_design_5000, fs2_ss1, fst_pop_strat = 0.2, tag = "beta_0.5")
	fs2_ss1 = get_G2G_setup(repnb, s_stratified = "full", associated_strains = "full", associated_populations = "full",beta = 0.25)
	res = test_G2G_setup(study_design_5000, fs2_ss1, fst_pop_strat = 0.2, tag = "beta_0.25")
	
	fs2_ss1_without_strat = get_G2G_setup(repnb, associated_strains = "full", associated_populations = "full", beta = 0.1)
	res = test_G2G_setup(study_design_5000, fs2_ss1_without_strat, tag = "beta_0.1")
	fs2_ss1_without_strat = get_G2G_setup(repnb, associated_strains = "full", associated_populations = "full", beta = 0.5)
	res = test_G2G_setup(study_design_5000, fs2_ss1_without_strat, tag = "beta_0.5")
	fs2_ss1_without_strat = get_G2G_setup(repnb, associated_strains = "full", associated_populations = "full", beta = 0.25)
	res = test_G2G_setup(study_design_5000, fs2_ss1_without_strat, tag = "beta_0.25")
	######
	
	
	######ss2
	ss2 = get_G2G_setup(repnb, s_stratified = "full", associated_strains = "full", associated_populations = "P1", beta = 0.25)
	res = test_G2G_setup(study_design_5000, ss2, fst_pop_strat = 0.2, tag = "sample_size_5000")
	res = test_G2G_setup(study_design_10000, ss2, fst_pop_strat = 0.2, tag = "sample_size_10000")
	
	ss2_with_vir_strat = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "full", associated_populations = "half", beta = 0.25)
	res = test_G2G_setup(study_design_5000, ss2_with_vir_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_5000")
	res = test_G2G_setup(study_design_10000, ss2_with_vir_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_10000")
	
	ss2_with_vir_strat_asymetric = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "full", associated_populations = "P1", beta = 0.25)
	res = test_G2G_setup(study_design_10000, ss2_with_vir_strat_asymetric, fst_pop_strat = 0.2, fst_strain_strat = 0.05)
	
	
	ss2_with_vir_strat_eq = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "full", associated_populations = "half", beta = 0.25)
	res = test_G2G_setup(study_design_5000, ss2_with_vir_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000")
	res = test_G2G_setup(study_design_10000, ss2_with_vir_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000")
	
	ss2_with_vir_strat_asymetric_eq = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "full", associated_populations = "P1", beta = 0.25)
	res = test_G2G_setup(study_design_10000, ss2_with_vir_strat_asymetric_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2)
	
	
	####fs3_ss3
	fs3_ss3 = get_G2G_setup(repnb, s_stratified = "full", associated_strains = "half", associated_populations = "full", beta = 0.25)
	
	res = test_G2G_setup(study_design_5000, fs3_ss3, fst_pop_strat = 0.2, tag = "sample_size_5000_beta_0.25")
	res = test_G2G_setup(study_design_10000, fs3_ss3, fst_pop_strat = 0.2, tag = "sample_size_10000_beta_5")
	fs3_ss3 = get_G2G_setup(repnb, s_stratified = "full", associated_strains = "half", associated_populations = "full", beta = 0.5)
	res = test_G2G_setup(study_design_5000, fs3_ss3, fst_pop_strat = 0.2, tag = "sample_size_5000_beta_0.5")
	res = test_G2G_setup(study_design_10000, fs3_ss3, fst_pop_strat = 0.2, tag = "sample_size_10000_beta_0.5")
	
	
	fs3_ss3_strat = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "full", beta = 0.25)
	res = test_G2G_setup(study_design_5000, fs3_ss3_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_5000")
	res = test_G2G_setup(study_design_10000, fs3_ss3_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_10000")
	
	fs3_ss3_strat = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "full", beta = 0.5)
	res = test_G2G_setup(study_design_5000, fs3_ss3_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_5000_beta_0.5")
	res = test_G2G_setup(study_design_10000, fs3_ss3_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_10000_beta_0.5")
	
	fs3_ss3_strat_eq = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "full", beta = 0.25)
	res = test_G2G_setup(study_design_5000, fs3_ss3_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000")
	res = test_G2G_setup(study_design_10000, fs3_ss3_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000")
	
	fs3_ss3_strat_eq = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "full", beta = 0.5)
	res = test_G2G_setup(study_design_5000, fs3_ss3_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.5")
	res = test_G2G_setup(study_design_10000, fs3_ss3_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000_beta_0.5")
	
	
	#####fs4
	fs4 = get_G2G_setup(repnb, s_stratified = "full",y_biased = "full", associated_strains = "half", associated_populations = "full", beta = 0.25)
	res = test_G2G_setup(study_design_5000, fs4, fst_pop_strat = 0.2, fst_strain_bias = 0.2, tag = "sample_size_5000_beta_0.25")
	res = test_G2G_setup(study_design_10000, fs4, fst_pop_strat = 0.2, fst_strain_bias = 0.2, tag = "sample_size_10000_beta_0.25")
	
	fs4 = get_G2G_setup(repnb, s_stratified = "full",y_biased = "full", associated_strains = "half", associated_populations = "full", beta = 0.5)
	res = test_G2G_setup(study_design_5000, fs4, fst_pop_strat = 0.2, fst_strain_bias = 0.2, tag = "sample_size_5000_beta_0.5")
	res = test_G2G_setup(study_design_10000, fs4, fst_pop_strat = 0.2, fst_strain_bias = 0.2, tag = "sample_size_10000_beta_0.5")
	
	fs4_strat = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full",y_biased = "full", associated_strains = "half", associated_populations = "full", beta = 0.25)
	res = test_G2G_setup(study_design_5000, fs4_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, fst_strain_bias = 0.05, tag = "sample_size_5000_beta_0.25")
	res = test_G2G_setup(study_design_10000, fs4_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, fst_strain_bias = 0.05, tag = "sample_size_10000_beta_0.25")
	
	fs4_strat = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full",y_biased = "full", associated_strains = "half", associated_populations = "full", beta = 0.5)
	res = test_G2G_setup(study_design_5000, fs4_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, fst_strain_bias = 0.2, tag = "sample_size_5000_beta_0.5")
	res = test_G2G_setup(study_design_10000, fs4_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, fst_strain_bias = 0.2, tag = "sample_size_10000_beta_0.5")
	
	fs4_strat_eq = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full",y_biased = "full", associated_strains = "half", associated_populations = "full", beta = 0.25)
	res = test_G2G_setup(study_design_5000, fs4_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.05, tag = "sample_size_5000_beta_0.25")
	res = test_G2G_setup(study_design_10000, fs4_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.05, tag = "sample_size_10000_beta_0.25")
	
	fs4_strat_eq = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full",y_biased = "full", associated_strains = "half", associated_populations = "full", beta = 0.5)
	res = test_G2G_setup(study_design_5000, fs4_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.2, tag = "sample_size_5000_beta_0.5")
	res = test_G2G_setup(study_design_10000, fs4_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.2, tag = "sample_size_10000_beta_0.5")
	
	
	####fs4 unbiased
	#fs4_unbiased = get_G2G_setup(repnb, s_stratified = "full",y_biased = "full", associated_strains = "half", associated_populations = "full", Ordered_Bias = FALSE)
	
	#res = test_G2G_setup(study_design_5000, fs4_unbiased, fst_pop_strat = 0.2, fst_strain_bias = 0.2, beta = 0.25)
	#res = test_G2G_setup(study_design_10000, fs4_unbiased, fst_pop_strat = 0.2, fst_strain_bias = 0.2, beta = 0.25)
	
	####ss4
	ss4 = get_G2G_setup(repnb, s_stratified = "full", associated_strains = "half", associated_populations = "half", beta = 0.25)
	res = test_G2G_setup(study_design_5000, ss4, fst_pop_strat = 0.2, tag = "sample_size_5000_beta_0.25")
	res = test_G2G_setup(study_design_10000, ss4, fst_pop_strat = 0.2, tag = "sample_size_10000_beta_0.25")
	
	ss4 = get_G2G_setup(repnb, s_stratified = "full", associated_strains = "half", associated_populations = "half", beta = 0.5)
	res = test_G2G_setup(study_design_5000, ss4, fst_pop_strat = 0.2, tag = "sample_size_5000_beta_0.5")
	res = test_G2G_setup(study_design_10000, ss4, fst_pop_strat = 0.2, tag = "sample_size_10000_beta_0.5")
	
	ss4_with_strat = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "half", beta = 0.25)
	res = test_G2G_setup(study_design_5000, ss4_with_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_5000_beta_0.25")
	res = test_G2G_setup(study_design_10000, ss4_with_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_10000_beta_0.25")
	
	ss4_with_strat = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "half", beta = 0.5)
	res = test_G2G_setup(study_design_5000, ss4_with_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_5000_beta_0.5")
	res = test_G2G_setup(study_design_10000, ss4_with_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_10000_beta_0.5")
	
	ss4_with_strat_eq = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "half", beta = 0.25)
	res = test_G2G_setup(study_design_5000, ss4_with_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.25")
	res = test_G2G_setup(study_design_10000, ss4_with_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000_beta_0.25")
	
	ss4_with_strat_eq = get_G2G_setup(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "half", beta = 0.5)
	res = test_G2G_setup(study_design_5000, ss4_with_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.5")
	res = test_G2G_setup(study_design_10000, ss4_with_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000_beta_0.5")
	
	######SC1
	sc1 = get_G2G_setup(repnb, s_stratified = c("P1","P2"), y_stratified = c("A","B"),y_biased = c("P1","P2"))
	res = test_G2G_setup(study_design_5000, sc1, fst_pop_strat = 0.2, fst_strain_strat = 0.02, fst_strain_bias = 0.2, tag = "*")
	#with vir>pop
	res = test_G2G_setup(study_design_5000, sc1, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.02, tag = "with_vir_pop")
	
	res = test_G2G_setup(study_design_5000, sc1, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.2, tag = "_eq")
	
	
	sc1_inv = get_G2G_setup(repnb, s_stratified = "P2", y_stratified = "full",y_biased = "full")
	res = test_G2G_setup(study_design_5000, sc1_inv, fst_pop_strat = 0.2, fst_strain_strat = 0.02, fst_strain_bias = 0.2, tag = "*")
	
	
	####sc1b
	sc1b = get_G2G_setup(repnb, s_stratified = c("P1","P2"), y_stratified = c("A","B"), y_biased = c("P1","P2"),  associated_strains = "A", associated_populations = "full", beta = 0.25)
	res = test_G2G_setup(study_design_5000, sc1b, fst_pop_strat = 0.2, fst_strain_strat = 0.02, fst_strain_bias = 0.2)
	res = test_G2G_setup(study_design_5000, sc1b, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.02, tag = "with_vir_pop")
	res = test_G2G_setup(study_design_5000, sc1b, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.2, tag = "_eq")
	
	
	####sc2
	study_design_sc2_5000 = to_study_design(as.data.frame(list(`P1` = c(`A` = 1000, `B` = 1500), `P2` = c(`A` = 1000, `B`  = 1500))))
	study_design_sc2_10000 = to_study_design(as.data.frame(list(`P1` = c(`A` = 2000, `B` = 3000), `P2` = c(`A` = 2000, `B`  = 3000))))
	
	sc2 = get_G2G_setup(repnb, s_stratified = c("P1","P2"), y_stratified = c("A","B"),  associated_strains = "A", associated_populations = "full", beta = 0.25)
	res = test_G2G_setup(study_design_sc2_5000, sc2, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.25")
	res = test_G2G_setup(study_design_sc2_10000, sc2, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000_beta_0.25")
	
	sc2 = get_G2G_setup(repnb, s_stratified = c("P1","P2"), y_stratified = c("A","B"),  associated_strains = "A", associated_populations = "full")
	res = test_G2G_setup(study_design_sc2_5000, sc2, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.5")
	res = test_G2G_setup(study_design_sc2_10000, sc2, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000_beta_0.5")
	
	####Demonstration
	demo_no_kill = get_G2G_setup(repnb, s_stratified = "full")
	res = test_G2G_setup(study_design_5000, demo_no_kill, fst_pop_strat = 0.2)
	
	demo_power_gain3 = get_G2G_setup(repnb, s_stratified = c("P1","P2"), associated_populations = c("P1","P2"))
	###0.01,0.2
	res = test_G2G_setup(study_design_5000, demo_power_gain3, fst_pop_strat = 0.2, fst_strain_bias = 0.2,beta = 0.25, get_viral = get_viral_output_sp)
	
	####sc2_ret
	sc2_ret = get_G2G_setup(repnb, s_stratified = c("P1","P2"), y_stratified = c("A","B"),  associated_strains = "A", associated_populations = "full", beta = 0.25)
	res = test_G2G_setup(study_design_sc2_5000, sc2_ret, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.25", get_viral = get_viral_output_sp3)
	res = test_G2G_setup(study_design_sc2_10000, sc2_ret, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000_beta_0.25", get_viral = get_viral_output_sp3)
	
	sc2_ret = get_G2G_setup(repnb, s_stratified = c("P1","P2"), y_stratified = c("A","B"),  associated_strains = "A", associated_populations = "full", beta = 0.5)
	res = test_G2G_setup(study_design_sc2_5000, sc2_ret, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.5", get_viral = get_viral_output_sp3)
	res = test_G2G_setup(study_design_sc2_10000, sc2_ret, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000_beta_0.5", get_viral = get_viral_output_sp3)
	
	sc2_ret_double_ass = get_G2G_setup(repnb, s_stratified = c("P1","P2"), y_stratified = c("A","B"), associated_strains = "full", associated_populations = "full", beta = 0.25)
	res = test_G2G_setup(study_design_sc2_5000, sc2_ret_double_ass, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.25", get_viral = get_viral_output_sp3)
	
	################Experimenting on parallel
}
####NOP
study_design_Ab = to_study_design(as.data.frame(list(`P1` = c(`A` = 1500, `B` = 1000), `P2` = c(`A` = 1000, `B`  = 1500))))
c1_cor =  get_G2G_setup(1000, s_stratified = c("P1","P2"), y_stratified = c("A","B"))
res = test_G2G_setup(study_design_Ab, c1_cor, fst_pop_strat = 0.2, fst_strain_strat = 0.2)

###imroved 2 
study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 1250, `B` = 1250), `P2` = c(`A` = 1, `B`  = 2499))))
c1_improved_2 = get_G2G_setup(1000, s_stratified = c("P1","P2"), y_stratified = c("A","B"))
res = test_G2G_setup(study_design, c1_improved_2, fst_pop_strat = 0.2, fst_pop_bias = 0.2, fst_strain_strat = 0.02)

study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 1250, `B` = 1250), `P2` = c(`A` = 1250, `B`  = 1250))))
c1_improved =  get_G2G_setup(1000, s_stratified = c("P1","P2"), y_stratified = c("A","B"), y_biased = c("P1","P2"), y_partial_strat = "P1")
res = test_G2G_setup(study_design, c1_improved, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.2)


study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 1250, `B` = 1250), `P2` = c(`A` = 1250, `B`  = 1250))))
fs2_ss1_cor_cm = get_G2G_setup(1000, s_stratified = "full", associated_strains = "full", associated_populations = "full",beta = 0.25)
res1 = test_G2G_setup(study_design, fs2_ss1_cor_cm, fst_pop_strat = 0.2, tag = "beta_0.25", get_viral = generate_AAs)
fs2_ss1_cor2_cm = get_G2G_setup(1000, associated_strains = "full", associated_populations = "full",beta = 0.25)
res2 = test_G2G_setup(study_design, fs2_ss1_cor2_cm, tag = "beta_0.25", get_viral = generate_AAs)



####FOR MASTER THESIS
##Neutral
study_design_Ab = to_study_design(as.data.frame(list(`P1` = c(`A` = 1500, `B` = 1000), `P2` = c(`A` = 1000, `B`  = 1500))))
study_design_P1A_P2B = to_study_design(as.data.frame(list(`P1` = c(`A` = 2500, `B` = 1), `P2` = c(`A` = 1, `B`  = 2500))))
study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 1250, `B` = 1250), `P2` = c(`A` = 1250, `B`  = 1250))))

study_design_P1AB_P2B = to_study_design(as.data.frame(list(`P1` = c(`A` = 2500, `B` = 0), `P2` = c(`A` = 0, `B`  = 2500))))
fs1_cor = get_G2G_setup(1000, s_stratified = c("P1", "P2"), y_stratified = c("A", "B"))
res = test_G2G_setup(study_design_P1AB_P2B, fs1_cor, fst_pop_strat = 0.2, fst_strain_strat = 0.2)

study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 1500, `B` = 1500), `P2` = c(`A` = 1500, `B`  = 1500))))
classic_strat =  get_G2G_setup(1000, s_stratified = c("P1","P2"), y_stratified = c("A","B"))
res = test_G2G_setup(study_design, classic_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.2)

study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 1125, `B` = 625), `P2` = c(`A` = 1125, `B`  = 2150))))
classic_strat_st2 =  get_G2G_setup(1000, s_stratified = c("P1","P2"), y_stratified = c("A","B"))
res = test_G2G_setup(study_design, classic_strat_st2, fst_pop_strat = 0.2, fst_strain_strat = 0.2)


study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 1250, `B` = 1250), `P2` = c(`A` = 1, `B`  = 2499))))
c1_ref_cor_2 = get_G2G_setup(1000, s_biased = c("A","B"), s_stratified = c("P1","P2"), s_partial_bias = "P1", y_stratified = c("A","B"), y_partial_strat = "P1")
res = test_G2G_setup(study_design, c1_ref_cor_2, fst_pop_strat = 0.2, fst_pop_bias = 0.2, fst_strain_strat = 0.02)

##Classic

##Associated

study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 1250, `B` = 1250), `P2` = c(`A` = 1250, `B`  = 1250))))
fs2_ss1_cor_cl = get_G2G_setup(1000, s_stratified = "full", associated_strains = "full", associated_populations = "full",beta = 0.25)
res1 = test_G2G_setup(study_design, fs2_ss1_cor_cl, fst_pop_strat = 0.2, tag = "beta_0.25", get_viral = generate_AAs_classic)
fs2_ss1_cor2_cl = get_G2G_setup(1000, associated_strains = "full", associated_populations = "full",beta = 0.25)
res2 = test_G2G_setup(study_design, fs2_ss1_cor2_cl, tag = "beta_0.25", get_viral = generate_AAs_classic)

study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 2500, `B` = 0), `P2` = c(`A` = 0, `B`  = 2500))))
fs2_ss1_cor3_cm = get_G2G_setup(1000, s_stratified = c("P1","P2"), y_bias = c("A","B"), associated_strains = "full", associated_populations = "full",beta = 0.25)
res2 = test_G2G_setup(study_design, fs2_ss1_cor3_cm, tag = "beta_0.25", fst_pop_strat = 0.2, fst_strain_bias = 0.2, get_viral = generate_AAs)

###
study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 1500, `B` = 1000), `P2` = c(`A` = 1500, `B`  = 1000))))
sc2_ret = get_G2G_setup(1000, s_stratified = c("P1","P2"), y_stratified = c("A","B"),  associated_strains = "A", associated_populations = "full", beta = 0.25)
#res = test_G2G_setup(study_design, sc2_ret, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.5_power_gain_old_work_ps", get_viral = generate_AAs_power_gain)
#res = test_G2G_setup(study_design, sc2_ret, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.5_power_gain_old_do_not_work_ps", get_viral = generate_AAs_power_gain_do_not_works_like_old)
res = test_G2G_setup(study_design, sc2_ret, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "ps", get_viral = generate_AAs_pg_new)


study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 1250, `B` = 1250), `P2` = c(`A` = 1250, `B`  = 1250))))
sc2_ret = get_G2G_setup(1000, s_stratified = c("P1","P2"), y_stratified = c("A","B"),  associated_strains = "A", associated_populations = "full", beta = 0.25)
#res = test_G2G_setup(study_design, sc2_ret, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.5_power_gain_old_work", get_viral = generate_AAs_power_gain)
#res = test_G2G_setup(study_design, sc2_ret, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.5_power_gain_old_do_not_work", get_viral = generate_AAs_power_gain_do_not_works_like_old)
res = test_G2G_setup(study_design, sc2_ret, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "s", get_viral = generate_AAs_pg_new)

###Study design do not change results


study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 1500, `B` = 1000), `P2` = c(`A` = 1500, `B`  = 1000))))
sc2_ret_double_ass = get_G2G_setup(1000, s_stratified = c("P1","P2"), y_stratified = c("A","B"), associated_strains = "full", associated_populations = "full", beta = 0.25)
res = test_G2G_setup(study_design, sc2_ret_double_ass, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "ps", get_viral = generate_AAs_pg_new_da)

study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 1250, `B` = 1250), `P2` = c(`A` = 1250, `B`  = 1250))))
sc2_ret_double_ass = get_G2G_setup(1000, s_stratified = c("P1","P2"), y_stratified = c("A","B"), associated_strains = "full", associated_populations = "full", beta = 0.25)
res = test_G2G_setup(study_design, sc2_ret_double_ass, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "s", get_viral = generate_AAs_pg_new_da)






res = replicate(30, {
	SNP.freq = get_frequencies(SNP.scenario, populations, strains)
	SNP.data = generate_SNPs_for_G2G(SNP.freq, study_design)
	sum(SNP.data)})