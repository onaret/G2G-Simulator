#rm(list = ls())
cat("\014")  
#library(rjson)
library(dplyr)
library(ggplot2)
library(tidyr)
library(parallel)

source("G2G.R")
source("GWAS.R")
source("summarize.R")
trace <- TRUE

####Full scenario
if(FALSE) {
  study_design = get_study_design(sample_size=5000, nb_strain=2, nb_pop=2)
  
  nb_cpu = 30
  
  #c1 = get_aa(stratified = "full", associated_strains = "half", associated_populations = "full", beta = 0.25, associated_SNPs = SNP_input, size = 10, study_design = study_design, pop_structure = pop_structure)
  
  ########FS1
  aa1= get_aa(study_design,25, fst_bias=0.02, biased = "full")
  aa2= get_aa(study_design,25, fst_bias=0.2, biased = "full")
  aa3= get_aa(study_design,25, fst_strat=0.02, stratified = "full")
  aa4 = get_aa(study_design,25, fst_strat=0.2, stratified = "full")
  
  aa5 = get_aa(study_design,25, fst_bias=0.02, fst_strat=0.2, stratified = "full", biased = "full")
  aa6 = get_aa(study_design,25, fst_bias=0.2, fst_strat=0.02, stratified = "full", biased = "full")
  
  aa7 = get_aa(study_design,25, fst_bias=0.02, fst_strat=0.2, stratified = "P1", biased = "A")
  aa8 = get_aa(study_design,25, fst_bias=0.2, fst_strat=0.02, stratified = "P2", biased = "A")
  
  SNP1 = get_SNP(study_design,10, stratified = "full", fst_strat=0.2)
  aa9 = get_aa(study_design,10, associated_strains = "full", associated_populations = "full", beta = 0.1, associated_SNPs = SNP1)
  aa10 = get_aa(study_design,10, associated_strains = "full", associated_populations = "full", beta = 0.25, associated_SNPs = SNP1)
  aa11 = get_aa(study_design,10, associated_strains = "full", associated_populations = "full", beta = 0.5, associated_SNPs = SNP1)
  
  SNP2 = get_SNP(study_design,10, stratified = "full", fst_strat=0.02)
  aa12 = get_aa(study_design,10, associated_strains = "full", associated_populations = "full", beta = 0.1, associated_SNPs = SNP2)
  aa13 = get_aa(study_design,10, associated_strains = "full", associated_populations = "full", beta = 0.25, associated_SNPs = SNP2)
  aa14 = get_aa(study_design,10, associated_strains = "full", associated_populations = "full", beta = 0.5, associated_SNPs = SNP2)
  
  SNP3 = get_SNP(study_design,10, stratified = "full", fst_strat=0.2)
  aa15 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, associated_strains = "full", associated_populations = "full", beta = 0.1, associated_SNPs = SNP3)
  aa16 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, associated_strains = "full", associated_populations = "full", beta = 0.25, associated_SNPs = SNP3)
  aa17 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, associated_strains = "full", associated_populations = "full", beta = 0.5, associated_SNPs = SNP3)
  
  SNP4 = get_SNP(study_design,10, stratified = "full", fst_strat=0.02)
  aa18 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, associated_strains = "full", associated_populations = "full", beta = 0.1, associated_SNPs = SNP4)
  aa19 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, associated_strains = "full", associated_populations = "full", beta = 0.25, associated_SNPs = SNP4)
  aa20 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, associated_strains = "full", associated_populations = "full", beta = 0.5, associated_SNPs = SNP4)
  
  SNP5 = get_SNP(study_design,10, stratified = "full", fst_strat=0.2)
  aa21 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, associated_strains = "A", associated_populations = "P2", beta = 0.1, associated_SNPs = SNP5)
  aa22 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, associated_strains = "B", associated_populations = "P1", beta = 0.25, associated_SNPs = SNP5)
  aa23 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, associated_strains = "half", associated_populations = "half", beta = 0.5, associated_SNPs = SNP5)
  
  SNP6 = get_SNP(study_design,10, stratified = "full", fst_strat=0.02)
  aa24 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, associated_strains = "half", associated_populations = "half", beta = 0.1, associated_SNPs = SNP6)
  aa25 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, associated_strains = "A", associated_populations = "P1", beta = 0.25, associated_SNPs = SNP6)
  aa26 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, associated_strains = "B", associated_populations = "P2", beta = 0.5, associated_SNPs = SNP6)
  
  SNP7 = get_SNP(study_design,10, stratified = "full", fst_strat=0.2)
  aa27 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, biased = "full", fst_bias=0.2, associated_strains = "A", associated_populations = "P2", beta = 0.1, associated_SNPs = SNP7)
  aa28 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, biased = "full", fst_bias=0.2, associated_strains = "B", associated_populations = "P1", beta = 0.25, associated_SNPs = SNP7)
  aa29 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, biased = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "half", beta = 0.5, associated_SNPs = SNP7)
  
  SNP8 = get_SNP(study_design,10, stratified = "full", fst_strat=0.02)
  aa30 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, biased = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "half", beta = 0.1, associated_SNPs = SNP8)
  aa31 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, biased = "full", fst_bias=0.2, associated_strains = "A", associated_populations = "P1", beta = 0.25, associated_SNPs = SNP8)
  aa32 = get_aa(study_design,10, stratified = "full", fst_strat=0.2, biased = "full", fst_bias=0.2, associated_strains = "B", associated_populations = "P2", beta = 0.5, associated_SNPs = SNP8)
  
  SNP9 = get_SNP(study_design,10, stratified = "full", fst_strat=0.2)
  aa33 = get_aa(study_design,10, stratified = "A", fst_strat=0.2, biased = "P1", fst_bias=0.2, associated_strains = "A", associated_populations = "P2", beta = 0.1, associated_SNPs = SNP9)
  aa34 = get_aa(study_design,10, stratified = "B", fst_strat=0.2, biased = "full", fst_bias=0.2, associated_strains = "B", associated_populations = "P1", beta = 0.25, associated_SNPs = SNP9)
  aa35 = get_aa(study_design,10, stratified = "A", fst_strat=0.2, biased = "P2", fst_bias=0.2, associated_strains = "half", associated_populations = "half", beta = 0.5, associated_SNPs = SNP9)
  
  SNP10 = get_SNP(study_design,10, stratified = "full", fst_strat=0.02)
  aa36 = get_aa(study_design,10, stratified = "B", fst_strat=0.2, biased = "P2", fst_bias=0.2, associated_strains = "half", associated_populations = "half", beta = 0.1, associated_SNPs = SNP10)
  aa37 = get_aa(study_design,10, stratified = "A", fst_strat=0.2, biased = "full", fst_bias=0.2, associated_strains = "A", associated_populations = "P1", beta = 0.25, associated_SNPs = SNP10)
  aa38 = get_aa(study_design,10, stratified = "B", fst_strat=0.2, biased = "P1", fst_bias=0.2, associated_strains = "B", associated_populations = "P2", beta = 0.5, associated_SNPs = SNP10)
  
  ######SC1
  aa39 = get_aa(study_design,10, biased = c("P1","P2"), fst_bias=0.2, stratified = c("A","B"), fst_strat=0.02)
  aa40 = get_aa(study_design,10, biased = c("P1","P2"), fst_bias=0.02, stratified = c("A","B"), fst_strat=0.2)
  aa41 = get_aa(study_design,10, biased = c("P1","P2"), fst_bias=0.2, stratified = c("A","B"), fst_strat=0.2)
  
  ####sc1b
  SNP12 = get_SNP(study_design,10, stratified = c("P1","P2"), fst_strat=0.2)
  aa42 = get_aa(study_design,10, biased = c("P1","P2"), fst_bias=0.2, stratified = c("A","B"), fst_strat=0.02, associated_strains = "A", associated_populations = "full", beta = 0.25, associated_SNPs = SNP12)
  aa43 = get_aa(study_design,10, biased = c("P1","P2"), fst_bias=0.02, stratified = c("A","B"), fst_strat=0.2, associated_strains = "A", associated_populations = "full", beta = 0.25, associated_SNPs = SNP12)
  aa44 = get_aa(study_design,10, biased = c("P1","P2"), fst_bias=0.2, stratified = c("A","B"), fst_strat=0.2, associated_strains = "A", associated_populations = "full", beta = 0.25, associated_SNPs = SNP12)
  
  SNPstrat = get_SNP(study_design,2500, stratified = "full", fst_strat=0.1)
  SNPunstrat = get_SNP(study_design,50000)
  aaunstrat = get_aa(study_design,1000)
  
  data = parse_G2G_config(aa1, aa2, aa3, aa4, aa5, aa6, aa7, aa8, SNP1, aa9, aa10, aa11, SNP2, aa12, aa13, aa14, SNP3, aa15, aa16, aa17, SNP4, aa18, aa19, aa20, SNP5, aa21, aa22, aa23, SNP6, aa24, aa25, aa26, SNP7, aa27, aa28, aa29, SNP8, aa30, aa31, aa32, SNP9, aa33, aa34, aa35, SNP10, aa36, aa37, aa38, aa39, aa40, aa41, SNP12, aa42, aa43, aa44, SNPstrat, SNPunstrat, aaunstrat)
  rm(aa1, aa2, aa3, aa4, aa5, aa6, aa7, aa8, SNP1, aa9, aa10, aa11, SNP2, aa12, aa13, aa14, SNP3, aa15, aa16, aa17, SNP4, aa18, aa19, aa20, SNP5, aa21, aa22, aa23, SNP6, aa24, aa25, aa26, SNP7, aa27, aa28, aa29, SNP8, aa30, aa31, aa32, SNP9, aa33, aa34, aa35, SNP10, aa36, aa37, aa38, aa39, aa40, aa41, SNP12, aa42, aa43, aa44, SNPstrat, SNPunstrat, aaunstrat)
  
  AA = data$AA
  SNP = data$SNP
  AA.scenarios = data$AA.scenarios
  SNP.scenarios = data$SNP.scenarios
  associations = data$associations
  rm(data)
  gc()
  res = analyse_G2G(SNP, AA, study_design, nb_cpu)
  plot_G2G(res, data.associations, data.AA.scenarios, data.SNP.scenarios)
}

####Different stratification scenarios
if(FALSE) {
  repnb = 1000
  study_design_5000 = get_study_design(5000,2,2)
  study_design_10000 = get_study_design(10000,2,2)
  
  ########C1
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
  ss2 = get_G2G_setup(repnb, s_stratified = "full", associated_strains = "full", associated_populations = "half", beta = 0.25)
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
  study_design_sc2_5000 = get_study_design(5000, 2, c(2, 0.35, 0.65))
  study_design_sc2_10000 = get_study_design(10000, 2, c(2, 0.35, 0.65))
  
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

####Analyse viral PC
if(FALSE) {
  ###Establish a scenario where viral stratification make false positive
  #study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 1500, `B` = 1000), `P2` = c(`A` = 1000, `B`  = 1500))))
  study_design = get_study_design(sample_size = 5000,nb_pop = 2,nb_strain = 2)
  #aa_low_f = get_aa(study_design,50, fst_strat=0.02, stratified = "full")
  aa_high_f = get_aa(study_design,50, fst_strat=0.2, stratified = "full")
  aa_unstrat = get_aa(study_design,950)
  SNP = get_SNP(study_design,fst_bias = 0.2, biased = "full", 1000)
  
  data = parse_config(aa_high_f, SNP, aa_unstrat)
  res = analyse_full(data$SNP, data$AA, study_design, nb_cpu = 10,WO_correction = T,W_strain_group = T, W_strain_PC = T)
  plot_G2G_exp(res, associations = data$associations, AA.scenarios = data$AA.scenarios, SNP.scenarios = data$SNP.scenarios)
  
  if(FALSE) {
    pc_low_f = prcomp(cbind(aa_unstrat$AA.data, aa_low_f$AA.data))
    pc_high_f = prcomp(cbind(aa_unstrat$AA.data, aa_high_f$AA.data))
    
    summary(glm(aa_low_f$AA.data~pc_low_f$x[,1:5]))
    summary(glm(aa_low_f$AA.data~study_design[,"Population"]))
    
    
    aa_high_f = get_aa(study_design,150, fst_strat=0.2, stratified = "full")
    aa_unstrat = get_aa(study_design,850)
    pc_high_f = prcomp(cbind(aa_unstrat$AA.data, aa_high_f$AA.data))
    summary(pc_high_f)
    
    
    SNP_strat = get_SNP(study_design, 5000,stratified = "full", fst_strat = 0.02)
    SNP_unstrat = get_SNP(study_design, 95000)
    pc_high_f = prcomp(cbind(SNP_strat$SNP.data, SNP_unstrat$SNP.data))
    head(summary(pc_high_f))
    
    SNP_strat1 = get_SNP(study_design, 5000,stratified = "full", fst_strat = 0.2)
    SNP_unstrat1 = get_SNP(study_design, 95000)
    pc_high_f = prcomp(cbind(SNP_strat1$SNP.data, SNP_unstra1t$SNP.data))
    head(summary(pc_high_f))
  }
  tester =  get_G2G_setup(100, y_stratified = c("A","B"), s_biased = c("A","B"))
  res = test_G2G_setup(study_design, tester, fst_strain_strat = 0.2, fst_pop_bias = 0.2)
}

###Generate SNPs
if(FALSE) {
  
  ####Constants
  C0 = generate_population_structure_for_CC(list(`P1` = c(`case` = 1250, `control` = 1250), `P2` = c(`case` = 1250, `control`  = 1250)))
  C1 = generate_population_structure_for_CC(list(`P1` = c(`case` = 200, `control` = 400), `P2` = c(`case` = 400, `control`  = 200)))
  C2 = generate_population_structure_for_CC(list(`P1` = c(`case` = 400, `control` = 200), `P2` = c(`case` = 200, `control` = 400)))
  C3 = generate_population_structure_for_CC(list(`P1` = c(`case` = 300, `control` = 0), `P2` = c(`case` = 300, `control` = 600)))
  C4 = generate_population_structure_for_CC(list(`P1` = c(`case` = 300, `control` = 200), `P2` = c(`case` = 200, `control` = 100), `P3` = c(`case` = 100, `control` = 300)))
  C5 = generate_population_structure_for_CC(list(`P1` = c(`case` = 200, `control` = 0), `P2` = c(`case` = 400, `control` = 200), `P3` = c(`case` = 0, `control` = 400)))
  AllPop = generate_population_structure_for_CC(list(`C1` = C1, `C2`= C2, `C3` = C3, `C4` = C4, `C5` = C5))
  
  fcoeff  = 0.2 ##### Wright's coefficient for inbreeding
  trace = TRUE
  
  #Para here
  load("G_Study_SNP_strat_correction.RData")
  cl = makeCluster(20, type = "FORK", outfile='outcluster.log')
  res = parLapply(cl, seq(0, 0.1,length.out = 100), function (strat_rate){
    GWAS_scenario(C1, 100000, strat_rate, causal_S = NULL, seq(1,2, by = 0.05), 0.01)
  })
  t3 = do.call(rbind, lapply(res, function(r) data.frame(`S_rate` = r$params$neutral_S_rate, `fst_strat` = r$params$fst_strat, r$summary_sim["lambda",], `FP` = r$summary_sim["FP_sum",])))  
  plot_GWAS(res[[306]]$pvalues, res[[306]]$SNP_params, title = paste("neutral rate: ", res[306]$params$neutral_S_rate)) 
  plot_GWAS(res[[76]]$pvalues, res[[76]]$SNP_params, title = paste("neutral rate: ", res[76]$params$neutral_S_rate))
  
  rest3 = filter(t3, fst_strat < 0.103 & fst_strat > 0.101)
  
  ###High strat, few SNP 2
  load("G_Study_SNP_strat_correction2.RData")
  cl = makeCluster(25, type = "FORK", outfile='outcluster.log')
  res = parLapply(cl, seq(0, 0.01,length.out = 30), function (strat_rate){
    GWAS_scenario(C1, 50000, strat_rate, causal_S = NULL, seq(1,2, by = 0.05), 0.2)
  })
  
  t4 = do.call(rbind, lapply(res, function(r) data.frame(`S_rate` = r$params$neutral_S_rate, `fst_strat` = r$params$fst_strat, r$summary_sim["lambda",], `FP` = r$summary_sim["FP_sum",])))  
  t4$strat_nb = t4$S_rate * 50000
  threshold <- 0.5/50000
  plot_GWAS(res[[1]]$pvalues, res[[1]]$SNP_params, title = paste("neutral rate: ", res[1]$params$neutral_S_rate)) 
  plot_GWAS(res[[2]]$pvalues, res[[2]]$SNP_params, title = paste("neutral rate: ", res[2]$params$neutral_S_rate)) 
  plot_GWAS(res[[3]]$pvalues, res[[3]]$SNP_params, title = paste("neutral rate: ", res[3]$params$neutral_S_rate)) 
  plot_GWAS(res[[4]]$pvalues, res[[4]]$SNP_params, title = paste("neutral rate: ", res[4]$params$neutral_S_rate)) 
  
  
  ###High strat, fewer SNP 3
  load("G_Study_SNP_strat_correction3.RData")
  cl = makeCluster(25, type = "FORK", outfile='outcluster.log')
  res = parLapply(cl, seq(0, 0.001,length.out = 30), function (strat_rate){
    GWAS_scenario(C1, 50000, strat_rate, causal_S = NULL, seq(1,2, by = 0.05), 0.2)
  })
  t5 = do.call(rbind, lapply(res, function(r) data.frame(`S_rate` = r$params$neutral_S_rate, `fst_strat` = r$params$fst_strat, r$summary_sim["lambda",], `FP` = r$summary_sim["FP_sum",])))  
  t5$strat_nb = t5$S_rate * 50000
  
  
  load( "G_Study_SNP_strat_correction4(best).RData")
  cl = makeCluster(25, type = "FORK", outfile='outcluster.log')
  res = parLapply(cl, seq(0, 0.002,length.out = 100), function (strat_rate){
    GWAS_scenario(C1, 50000, strat_rate, causal_S = NULL, seq(1,2, by = 0.05), 0.2)
  })
  nb_pc = 5
  analyse_FP_in_function_of_s_rate(res)
  
  
  
  load("G_Study_SNP_strat_correction4(best+PC).RData")
  nb_SNP = 50000
  nb_pc = 10
  cl = makeCluster(25, type = "FORK", outfile='outcluster.log')
  res = parLapply(cl, seq(0, 0.002,length.out = 100), function (strat_rate){
    GWAS_scenario(C1, nb_SNP, strat_rate, causal_S = NULL, seq(1,2, by = 0.05), 0.2,nb_pc = nb_pc)
  })
  
  analyse_FP_in_function_of_s_rate(res)
  
  analyse_FP_in_function_of_s_rate <- function(res) {
    t7 = do.call(rbind, lapply(res, function(r) data.frame(`S_rate` = r$params$neutral_S_rate, `FP` = r$summary_sim["FP_sum",])))  
    t7 = t7  %>% gather(S_rate, FP)
    colnames(t7) <- c("S_rate", "Correction", "FP")
    t7$FP <- t7$FP/nb_SNP
    p <- ggplot(data = t7, aes(S_rate, FP, colour=Correction))
    p + geom_smooth() +  
      labs(title = paste("FP in function of stratified SNPs number with", res[[1]]$params$fst_strat, "Fst and", nb_pc,"PCs"), y = "False positive per SNP", x=paste("Stratification rate in percentage for", nb_SNP,"SNP")) }
  
  nb_SNP = 50000
  nb_pc = 10
  cl = makeCluster(25, type = "FORK", outfile='outcluster.log')
  res = parLapply(cl, seq(0, 0.01,length.out = 500), function (strat_rate){
    GWAS_scenario(C1, nb_SNP, strat_rate, causal_S = NULL, seq(1,2, by = 0.05), 0.2,nb_pc = nb_pc)
  })
  analyse_FP_in_function_of_s_rate(res)
  
  res = GWAS_scenario(C1, nb_SNP, 0.01, causal_S = NULL, causal_NS = NULL, 0.2,nb_pc = nb_pc)
  
}

#####ANalyse full scenario

load("gen-data/big_first_res.RData")