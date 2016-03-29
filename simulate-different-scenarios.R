#rm(list = ls())
cat("\014")  
library(rjson)
library(parallel)
library(dplyr)
library(ggplot2)

source("generate-full-model.R")
source("summarize.R")

#######With virus
trace <- FALSE

###Input
threshold = 5e-05
repnb = 1000
study_design_5000 = get_study_design(5000,2,2)
study_design_10000 = get_study_design(10000,2,2)

########C1
c1 =  get_scenario(100, s_biased = c("A","B"), s_partial_bias = "P1", y_stratified = c("A","B"), y_partial_strat = "P1")
res = test_scenario(study_design_5000, c1, fst_pop_bias = 0.03, fst_strain_strat = 0.2)

c1 =  get_scenario(10, s_stratified = c("P1"," P2"), s_biased = c("A","B"), s_partial_bias = "P1", y_stratified = c("A","B"), y_partial_strat = "P1", y_biased = c("P1","P2"))
res = test_scenario(study_design_5000, c1, fst_pop_strat = 0.03, fst_pop_bias = 0.03, fst_strain_strat = 0.2, fst_strain_bias = 0.2)


c1_eq = get_scenario(repnb,s_biased = c("A","B"), s_partial_bias = "P1", y_stratified = c("A","B"), y_partial_strat = "P1")
res = test_scenario(study_design_5000, c1_eq, fst_pop_bias = 0.2, fst_strain_strat = 0.2)

c1_without_human_bias = get_scenario(repnb, s_stratified = c("P1"," P2"), y_stratified = c("A","B"), y_partial_strat = "P1")
res = test_scenario(study_design_5000, c1_without_human_bias, fst_pop_strat = 0.2, fst_strain_strat = 0.2)

c1_with_full_bias_full_strat = get_scenario(repnb, s_stratified =  c("P1"," P2"), s_biased = c("A","B"), y_stratified = c("A","B"), y_biased =  c("P1"," P2"))
res = test_scenario(study_design_5000, c1_with_full_bias_full_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.03, fst_pop_bias = 0.03)

c1_with_full_bias_full_strat_eq = get_scenario(repnb, s_stratified =  c("P1"," P2"), s_biased =c("A","B"), y_stratified = c("A","B"), y_biased =  c("P1"," P2"))
res = test_scenario(study_design_5000, c1_with_full_bias_full_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.2, fst_pop_bias = 0.2)


########FS1
fs1 = get_scenario(repnb, s_stratified = c("P1", "P2"), y_biased = c("P1", "P2"))
res = test_scenario(study_design_5000, fs1, fst_pop_strat = 0.2, fst_strain_bias = 0.02)

fs1_full_strat = get_scenario(repnb, s_stratified = "full",y_stratified = "full")
res = test_scenario(study_design_5000, fs1_full_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.2)

######fs2_ss1 
fs2_ss1 = get_scenario(repnb, s_stratified = "full", associated_strains = "full", associated_populations = "full",beta = 0.1)
res = test_scenario(study_design_5000, fs2_ss1, fst_pop_strat = 0.2, tag = "beta_0.1")
fs2_ss1 = get_scenario(repnb, s_stratified = "full", associated_strains = "full", associated_populations = "full",beta = 0.5)
res = test_scenario(study_design_5000, fs2_ss1, fst_pop_strat = 0.2, tag = "beta_0.5")
fs2_ss1 = get_scenario(repnb, s_stratified = "full", associated_strains = "full", associated_populations = "full",beta = 0.25)
res = test_scenario(study_design_5000, fs2_ss1, fst_pop_strat = 0.2, tag = "beta_0.25")

fs2_ss1_without_strat = get_scenario(repnb, associated_strains = "full", associated_populations = "full", beta = 0.1)
res = test_scenario(study_design_5000, fs2_ss1_without_strat, tag = "beta_0.1")
fs2_ss1_without_strat = get_scenario(repnb, associated_strains = "full", associated_populations = "full", beta = 0.5)
res = test_scenario(study_design_5000, fs2_ss1_without_strat, tag = "beta_0.5")
fs2_ss1_without_strat = get_scenario(repnb, associated_strains = "full", associated_populations = "full", beta = 0.25)
res = test_scenario(study_design_5000, fs2_ss1_without_strat, tag = "beta_0.25")
######


######ss2
ss2 = get_scenario(repnb, s_stratified = "full", associated_strains = "full", associated_populations = "half", beta = 0.25)
res = test_scenario(study_design_5000, ss2, fst_pop_strat = 0.2, tag = "sample_size_5000")
res = test_scenario(study_design_10000, ss2, fst_pop_strat = 0.2, tag = "sample_size_10000")

ss2_with_vir_strat = get_scenario(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "full", associated_populations = "half", beta = 0.25)
res = test_scenario(study_design_5000, ss2_with_vir_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_5000")
res = test_scenario(study_design_10000, ss2_with_vir_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_10000")

ss2_with_vir_strat_asymetric = get_scenario(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "full", associated_populations = "P1", beta = 0.25)
res = test_scenario(study_design_10000, ss2_with_vir_strat_asymetric, fst_pop_strat = 0.2, fst_strain_strat = 0.05)


ss2_with_vir_strat_eq = get_scenario(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "full", associated_populations = "half", beta = 0.25)
res = test_scenario(study_design_5000, ss2_with_vir_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000")
res = test_scenario(study_design_10000, ss2_with_vir_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000")

ss2_with_vir_strat_asymetric_eq = get_scenario(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "full", associated_populations = "P1", beta = 0.25)
res = test_scenario(study_design_10000, ss2_with_vir_strat_asymetric_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2)


####fs3_ss3
fs3_ss3 = get_scenario(repnb, s_stratified = "full", associated_strains = "half", associated_populations = "full", beta = 0.25)

res = test_scenario(study_design_5000, fs3_ss3, fst_pop_strat = 0.2, tag = "sample_size_5000_beta_0.25")
res = test_scenario(study_design_10000, fs3_ss3, fst_pop_strat = 0.2, tag = "sample_size_10000_beta_5")
fs3_ss3 = get_scenario(repnb, s_stratified = "full", associated_strains = "half", associated_populations = "full", beta = 0.5)
res = test_scenario(study_design_5000, fs3_ss3, fst_pop_strat = 0.2, tag = "sample_size_5000_beta_0.5")
res = test_scenario(study_design_10000, fs3_ss3, fst_pop_strat = 0.2, tag = "sample_size_10000_beta_0.5")


fs3_ss3_strat = get_scenario(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "full", beta = 0.25)
res = test_scenario(study_design_5000, fs3_ss3_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_5000")
res = test_scenario(study_design_10000, fs3_ss3_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_10000")

fs3_ss3_strat = get_scenario(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "full", beta = 0.5)
res = test_scenario(study_design_5000, fs3_ss3_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_5000_beta_0.5")
res = test_scenario(study_design_10000, fs3_ss3_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_10000_beta_0.5")

fs3_ss3_strat_eq = get_scenario(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "full", beta = 0.25)
res = test_scenario(study_design_5000, fs3_ss3_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000")
res = test_scenario(study_design_10000, fs3_ss3_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000")

fs3_ss3_strat_eq = get_scenario(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "full", beta = 0.5)
res = test_scenario(study_design_5000, fs3_ss3_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.5")
res = test_scenario(study_design_10000, fs3_ss3_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000_beta_0.5")


#####fs4
fs4 = get_scenario(repnb, s_stratified = "full",y_biased = "full", associated_strains = "half", associated_populations = "full", beta = 0.25)
res = test_scenario(study_design_5000, fs4, fst_pop_strat = 0.2, fst_strain_bias = 0.2, tag = "sample_size_5000_beta_0.25")
res = test_scenario(study_design_10000, fs4, fst_pop_strat = 0.2, fst_strain_bias = 0.2, tag = "sample_size_10000_beta_0.25")

fs4 = get_scenario(repnb, s_stratified = "full",y_biased = "full", associated_strains = "half", associated_populations = "full", beta = 0.5)
res = test_scenario(study_design_5000, fs4, fst_pop_strat = 0.2, fst_strain_bias = 0.2, tag = "sample_size_5000_beta_0.5")
res = test_scenario(study_design_10000, fs4, fst_pop_strat = 0.2, fst_strain_bias = 0.2, tag = "sample_size_10000_beta_0.5")

fs4_strat = get_scenario(repnb, s_stratified = "full",y_stratified = "full",y_biased = "full", associated_strains = "half", associated_populations = "full", beta = 0.25)
res = test_scenario(study_design_5000, fs4_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, fst_strain_bias = 0.05, tag = "sample_size_5000_beta_0.25")
res = test_scenario(study_design_10000, fs4_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, fst_strain_bias = 0.05, tag = "sample_size_10000_beta_0.25")

fs4_strat = get_scenario(repnb, s_stratified = "full",y_stratified = "full",y_biased = "full", associated_strains = "half", associated_populations = "full", beta = 0.5)
res = test_scenario(study_design_5000, fs4_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, fst_strain_bias = 0.2, tag = "sample_size_5000_beta_0.5")
res = test_scenario(study_design_10000, fs4_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, fst_strain_bias = 0.2, tag = "sample_size_10000_beta_0.5")

fs4_strat_eq = get_scenario(repnb, s_stratified = "full",y_stratified = "full",y_biased = "full", associated_strains = "half", associated_populations = "full", beta = 0.25)
res = test_scenario(study_design_5000, fs4_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.05, tag = "sample_size_5000_beta_0.25")
res = test_scenario(study_design_10000, fs4_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.05, tag = "sample_size_10000_beta_0.25")

fs4_strat_eq = get_scenario(repnb, s_stratified = "full",y_stratified = "full",y_biased = "full", associated_strains = "half", associated_populations = "full", beta = 0.5)
res = test_scenario(study_design_5000, fs4_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.2, tag = "sample_size_5000_beta_0.5")
res = test_scenario(study_design_10000, fs4_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.2, tag = "sample_size_10000_beta_0.5")


####fs4 unbiased
#fs4_unbiased = get_scenario(repnb, s_stratified = "full",y_biased = "full", associated_strains = "half", associated_populations = "full", Ordered_Bias = FALSE)

#res = test_scenario(study_design_5000, fs4_unbiased, fst_pop_strat = 0.2, fst_strain_bias = 0.2, beta = 0.25)
#res = test_scenario(study_design_10000, fs4_unbiased, fst_pop_strat = 0.2, fst_strain_bias = 0.2, beta = 0.25)

####ss4
ss4 = get_scenario(repnb, s_stratified = "full", associated_strains = "half", associated_populations = "half", beta = 0.25)
res = test_scenario(study_design_5000, ss4, fst_pop_strat = 0.2, tag = "sample_size_5000_beta_0.25")
res = test_scenario(study_design_10000, ss4, fst_pop_strat = 0.2, tag = "sample_size_10000_beta_0.25")

ss4 = get_scenario(repnb, s_stratified = "full", associated_strains = "half", associated_populations = "half", beta = 0.5)
res = test_scenario(study_design_5000, ss4, fst_pop_strat = 0.2, tag = "sample_size_5000_beta_0.5")
res = test_scenario(study_design_10000, ss4, fst_pop_strat = 0.2, tag = "sample_size_10000_beta_0.5")

ss4_with_strat = get_scenario(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "half", beta = 0.25)
res = test_scenario(study_design_5000, ss4_with_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_5000_beta_0.25")
res = test_scenario(study_design_10000, ss4_with_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_10000_beta_0.25")

ss4_with_strat = get_scenario(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "half", beta = 0.5)
res = test_scenario(study_design_5000, ss4_with_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_5000_beta_0.5")
res = test_scenario(study_design_10000, ss4_with_strat, fst_pop_strat = 0.2, fst_strain_strat = 0.05, tag = "sample_size_10000_beta_0.5")

ss4_with_strat_eq = get_scenario(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "half", beta = 0.25)
res = test_scenario(study_design_5000, ss4_with_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.25")
res = test_scenario(study_design_10000, ss4_with_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000_beta_0.25")

ss4_with_strat_eq = get_scenario(repnb, s_stratified = "full",y_stratified = "full", associated_strains = "half", associated_populations = "half", beta = 0.5)
res = test_scenario(study_design_5000, ss4_with_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.5")
res = test_scenario(study_design_10000, ss4_with_strat_eq, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000_beta_0.5")

######SC1
sc1 = get_scenario(repnb, s_stratified = c("P1","P2"), y_stratified = c("A","B"),y_biased = c("P1","P2"))
res = test_scenario(study_design_5000, sc1, fst_pop_strat = 0.2, fst_strain_strat = 0.02, fst_strain_bias = 0.2, tag = "*")
#with vir>pop
res = test_scenario(study_design_5000, sc1, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.02, tag = "with_vir_pop")

res = test_scenario(study_design_5000, sc1, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.2, tag = "_eq")


sc1_inv = get_scenario(repnb, s_stratified = "P2", y_stratified = "full",y_biased = "full")
res = test_scenario(study_design_5000, sc1_inv, fst_pop_strat = 0.2, fst_strain_strat = 0.02, fst_strain_bias = 0.2, tag = "*")


####sc1b
sc1b = get_scenario(repnb, s_stratified = c("P1","P2"), y_stratified = c("A","B"), y_biased = c("P1","P2"),  associated_strains = "A", associated_populations = "full", beta = 0.25)
res = test_scenario(study_design_5000, sc1b, fst_pop_strat = 0.2, fst_strain_strat = 0.02, fst_strain_bias = 0.2)
res = test_scenario(study_design_5000, sc1b, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.02, tag = "with_vir_pop")
res = test_scenario(study_design_5000, sc1b, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.2, tag = "_eq")


####sc2
study_design_sc2_5000 = get_study_design(5000, 2, c(2, 0.35, 0.65))
study_design_sc2_10000 = get_study_design(10000, 2, c(2, 0.35, 0.65))

sc2 = get_scenario(repnb, s_stratified = c("P1","P2"), y_stratified = c("A","B"),  associated_strains = "A", associated_populations = "full", beta = 0.25)
res = test_scenario(study_design_sc2_5000, sc2, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.25")
res = test_scenario(study_design_sc2_10000, sc2, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000_beta_0.25")

sc2 = get_scenario(repnb, s_stratified = c("P1","P2"), y_stratified = c("A","B"),  associated_strains = "A", associated_populations = "full")
res = test_scenario(study_design_sc2_5000, sc2, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.5")
res = test_scenario(study_design_sc2_10000, sc2, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000_beta_0.5")

####Demonstration
demo_no_kill = get_scenario(repnb, s_stratified = "full")
res = test_scenario(study_design_5000, demo_no_kill, fst_pop_strat = 0.2)

demo_power_gain3 = get_scenario(repnb, s_stratified = c("P1","P2"), associated_populations = c("P1","P2"))
###0.01,0.2
res = test_scenario(study_design_5000, demo_power_gain3, fst_pop_strat = 0.2, fst_strain_bias = 0.2,beta = 0.25, get_viral = get_viral_output_sp)

####sc2_ret
sc2_ret = get_scenario(repnb, s_stratified = c("P1","P2"), y_stratified = c("A","B"),  associated_strains = "A", associated_populations = "full", beta = 0.25)
res = test_scenario(study_design_sc2_5000, sc2_ret, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.25", get_viral = get_viral_output_sp3)
res = test_scenario(study_design_sc2_10000, sc2_ret, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000_beta_0.25", get_viral = get_viral_output_sp3)

sc2_ret = get_scenario(repnb, s_stratified = c("P1","P2"), y_stratified = c("A","B"),  associated_strains = "A", associated_populations = "full", beta = 0.5)
res = test_scenario(study_design_sc2_5000, sc2_ret, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.5", get_viral = get_viral_output_sp3)
res = test_scenario(study_design_sc2_10000, sc2_ret, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_10000_beta_0.5", get_viral = get_viral_output_sp3)

sc2_ret_double_ass = get_scenario(repnb, s_stratified = c("P1","P2"), y_stratified = c("A","B"), associated_strains = "full", associated_populations = "full", beta = 0.25)
res = test_scenario(study_design_sc2_5000, sc2_ret_double_ass, fst_pop_strat = 0.2, fst_strain_strat = 0.2, tag = "sample_size_5000_beta_0.25", get_viral = get_viral_output_sp3)

################Experimenting on parallel

cl = makeCluster(30, type = "FORK")
#parLapply(cl, seq(200, 5000, length.out = 10), function(sample_size) {
#  res = res = test_scenario(sample_size, nb_strain, nb_pop, causal_S =  causal_S, causal_NS = causal_NS , viral_aa = viral_aa,beta =  beta, fst_pop_strat =  fst_pop_strat, fst_strain_bias = fst_strain_bias) })

#size = seq(200, 5000, length.out = 10)
#fst_pop_strat = seq(0.01, 0.3, length.out = 10)
#beta = seq(0.1, 2, length.out = 10)


sc2_simple_implem = get_scenario(repnb, s_stratified = "P1", y_stratified = "A", y_biased = "hB",  associated_strains = "A", associated_populations = "full", size = 100)
res = res = test_scenario(sample_size = 100,  nb_strains = c(2, 0.35, 0.65), nb_pop, sc2_simple_implem, fst_pop_strat = 0.2, fst_strain_strat = 0.2, fst_strain_bias = 0.05, beta = 0.25)


res = res = test_scenario(study_design_5000, nb_strain, nb_pop, fs3_ss3, fst_pop_strat = 0.2, beta = 0.5, tag = "sample_size_5000_beta_0.5")


scale_y_continuous(limits=c(0,100), breaks=seq(0,100,10), expand = c(0, 0))

parLapply(beta, seq(0.1, 2, length.out = 30), function(sample_size) {
  res = test_scenario(sample_size, nb_strain, nb_pop, SNPs, fst_pop_strat, fst_strain_strat, fst_strain_bias, fst_pop_bias, beta) })

res = parLapply(cl, seq(0.01, 0.2, length.out = 30), function(param) {
  test = select(as.data.frame(test_scenario(sample_size, nb_strain, nb_pop, fs1_full_strat, fst_pop_strat = param, fst_strain_strat = param, tag = "fs1_strat_impact")$pval ), ends_with("pval"))
  medianpval = apply(pval, 2, function(p) median(p))
  rm(test,pval)
  c(param, medianpval)
})

res = parLapply(cl, seq(0.01, 1, length.out = 100), function(param) {
  res = test_scenario(sample_size, nb_strain, nb_pop, fs1, fst_pop_strat = param, fst_strain_bias = param, tag = "fs1_strat_impact_fp")
  test = unlist(res$summary_sim["FP_ratio",])----------------------------------------
    rm(res)
  c(param, test)
})


test_fp = unlist(res$summary_sim["FP_ratio",])

write.table(res,file = "result_on_scenario_fs1_bias")


res = test_scenario(sample_size, nb_strain, nb_pop, fs1, fst_pop_strat = param, fst_strain_bias = param, tag = "fs1_strat_impact_fp")

parLapply(beta, seq(0.1, 2, length.out = 30), function(sample_size) {
  res = test_scenario(sample_size, nb_strain, nb_pop, SNPs, fst_pop_strat, fst_strain_strat, fst_strain_bias, fst_pop_bias, beta) })

parLapply(cl, seq(200, 5000, length.out = 10), function(sample_size) {
  for(fst_pop_strat in seq(0.01, 0.3, length.out = 4)) {
    for(fst_strain_strat in seq(0.01, 0.3, length.out = 4)) {
      for(beta in seq(0.1, 2, length.out = 4)) {
        for(fst_pop_bias in seq(0.01, 0.3, length.out = 4)) {
          for(fst_strain_bias in seq(0.01, 0.3, length.out = 4)) {
            res = test_scenario(sample_size, nb_strain, nb_pop, SNPs, fst_pop_strat, fst_strain_strat, fst_strain_bias, fst_pop_bias, beta) }}}}}})