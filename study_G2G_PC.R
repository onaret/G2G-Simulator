#rm(list = ls())
cat("\014")  
#library(rjson)
library(dplyr)
library(ggplot2)
library(tidyr)
library(parallel)

source("G2G.R")
source("G2G_full.R")
source("G2G_single.R")
trace <- TRUE

####Analyse viral PC
if(FALSE) {
	###Establish a scenario where viral stratification make false positive
	#study_design = to_study_design(as.data.frame(list(`P1` = c(`A` = 1500, `B` = 1000), `P2` = c(`A` = 1000, `B`  = 1500))))
	study_design = get_study_design(sample_size = 5000,nb_pop = 2,nb_strain = 2)
	#aa_low_f = get_AA(study_design,50, fst_strat=0.02, stratified = "full")
	aa_high_f = get_AA(study_design,50, fst_strat=0.2, stratified = "full")
	aa_unstrat = get_AA(study_design,950)
	SNP = get_SNP(study_design,fst_bias = 0.2, biased = "full", 1000)
	
	data = parse_G2G_config(aa_high_f, SNP, aa_unstrat)
	res = analyse_G2G(data$SNP, data$AA, study_design, nb_cpu = 10,WO_correction = T,W_strain_group = T, W_strain_PC = T)
	
	plot_G2G_exp(res, associations = data$associations, AA.scenarios = data$AA.scenarios, SNP.scenarios = data$SNP.scenarios)
	
	library(homals)
	res = analyse_G2G(data$SNP, data$AA, study_design, nb_cpu = 10,WO_correction = T,W_strain_group = T, W_strain_PC = T)
	
	if(FALSE) {
		pc_low_f = prcomp(cbind(aa_unstrat$AA.data, aa_low_f$AA.data))
		pc_high_f = prcomp(cbind(aa_unstrat$AA.data, aa_high_f$AA.data))
		
		summary(glm(aa_low_f$AA.data~pc_low_f$x[,1:5]))
		summary(glm(aa_low_f$AA.data~study_design[,"Population"]))
		
		aa_high_f = get_AA(study_design,150, fst_strat=0.2, stratified = "full")
		aa_unstrat = get_AA(study_design,850)
		pc_high_f = prcomp(cbind(aa_unstrat$AA.data, aa_high_f$AA.data))
		summary(pc_high_f)
		
		SNP_strat = get_SNP(study_design, 5000,stratified = "full", fst_strat = 0.02)
		SNP_unstrat = get_SNP(study_design, 95000)
		pc_high_f = prcomp(cbind(SNP_strat$SNP.data, SNP_unstrat$SNP.data))
		head(summary(pc_high_f))
		
		SNP_strat1 = get_SNP(study_design, 5000,stratified = "full", fst_strat = 0.2)
		SNP_unstrat1 = get_SNP(study_design, 95000)
		pc_high_f = prcomp(cbind(SNP_strat1$SNP.data, SNP_unstra1t$SNP.data))
		head(summary(pc_high_f))}
}


data_base =	parse_G2G_config(
	study_design,
	G2G_conf(
		SNP(1000, fst_bias = 0.2, biased = "full"),
		AA(50,fst_strat = 0.2,stratified = "full"),
		AA(950)))

library(homals)
res = analyse_G2G(data$SNP, data$AA, study_design, nb_cpu = 10, WO_correction = T, W_strain_group = T, W_strain_PC = T)
