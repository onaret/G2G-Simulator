#rm(list = ls())
cat("\014")  
library(rjson)
library(parallel)
library(dplyr)
library(ggplot2)
library(reshape2)
require(gridExtra)

source("generate-full-model.R")
source("summarize.R")

#######With virus
trace <- FALSE

###Input
threshold = 5e-05
study_design = get_study_design(sample_size=5000, nb_strain=2, nb_pop=2)

if(FALSE) {
  trace = TRUE
  

  SNP_input = get_SNP(study_design, 10, stratified = "full", fst_strat = 0.01, biased = "full", fst_bias = 0.02)
  SNP_input2 = get_SNP(study_design, 10, stratified = c("P2", "P1"), fst_strat = 0.01, biased = "B", fst_bias = 0.02)
  
  aa_input = get_aa(study_design, 10, stratified = "full", fst_strat = 0.03, associated_strains = "half", associated_populations = "full", beta = 0.25, associated_SNPs = SNP_input)
  aa_input2 = get_aa(study_design, 10, biased = "full", fst_bias = 0.04)
  
  res = scenario_full(study_design, aa_input,aa_input2, SNP_input,SNP_input2)
}


nb_cpu = 30

#c1 = get_aa(stratified = "full", associated_strains = "half", associated_populations = "full", beta = 0.25, associated_SNPs = SNP_input, size = 10, study_design = study_design, pop_structure = pop_structure)

########FS1
aa = get_aa(study_design,10, fst_bias=0.02, biased = "full")

SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)

aa = get_aa(study_design,size, fst_strat=0.2, stratified = "full")
SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)

######fs2_ss1 
SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design,size, associated_strains = "full", associated_populations = "full", beta = 0.1, associated_SNPs = SNP)
aa = get_aa(study_design,size, associated_strains = "full", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
aa = get_aa(study_design,size, associated_strains = "full", associated_populations = "full", beta = 0.5, associated_SNPs = SNP)
SNP = get_SNP(study_design,size)
aa = get_aa(study_design,size, associated_strains = "full", associated_populations = "full", beta = 0.1, associated_SNPs = SNP)
aa = get_aa(study_design,size, associated_strains = "full", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
aa = get_aa(study_design,size, associated_strains = "full", associated_populations = "full", beta = 0.5, associated_SNPs = SNP)

######ss2
SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design,size, associated_strains = "full", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
res = scenario_full(aa, SNP, study_design=study_design)unlist(lapply(res, function(r) -log10(head(sort(r$W_human_groups.pval),n=1))))

SNP = get_SNP(study_design_10000,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design_10000,size, associated_strains = "full", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
SNP = get_SNP(study_design_10000,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design_10000,size, associated_strains = "full", associated_populations = "full", beta = 0.5, associated_SNPs = SNP)
SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design,size, stratified ="full", fst_strat=0.05, associated_strains = "full", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
SNP = get_SNP(study_design_10000,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design_10000,size, stratified ="full", fst_strat=0.05, associated_strains = "full", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)

SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design,size, stratified ="full", fst_strat=0.05, associated_strains = "full", associated_populations = "P1", beta = 0.25, associated_SNPs = SNP)
SNP = get_SNP(study_design_10000,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design_10000,size, stratified ="full", fst_strat=0.05, associated_strains = "full", associated_populations = "P1", beta = 0.25, associated_SNPs = SNP)

SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design,size, stratified ="full", fst_strat=0.2, associated_strains = "full", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
SNP = get_SNP(study_design_10000,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design_10000,size, stratified ="full", fst_strat=0.2, associated_strains = "full", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)

SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design,size, stratified ="full", fst_strat=0.2, associated_strains = "full", associated_populations = "P1", beta = 0.25, associated_SNPs = SNP)
SNP = get_SNP(study_design_10000,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design_10000,size, stratified ="full", fst_strat=0.2, associated_strains = "full", associated_populations = "P1", beta = 0.25, associated_SNPs = SNP)

####fs3_ss3
SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design,size, associated_strains = "half", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
SNP = get_SNP(study_design_10000,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design_10000,size, associated_strains = "half", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design,size, associated_strains = "half", associated_populations = "full", beta = 0.5, associated_SNPs = SNP)
SNP = get_SNP(study_design_10000,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design_10000,size, associated_strains = "half", associated_populations = "full", beta = 0.5, associated_SNPs = SNP)

SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design,size, stratified = "full", fst_strat=0.2, associated_strains = "half", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
SNP = get_SNP(study_design_10000,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design_10000,size, stratified = "full", fst_strat=0.2, associated_strains = "half", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design,size, stratified = "full", fst_strat=0.2, associated_strains = "half", associated_populations = "full", beta = 0.5, associated_SNPs = SNP)
SNP = get_SNP(study_design_10000,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design_10000,size, stratified = "full", fst_strat=0.2, associated_strains = "half", associated_populations = "full", beta = 0.5, associated_SNPs = SNP)
#####fs4
SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design,size, biased = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
aa = get_aa(study_design,size, biased = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "full", beta = 0.5, associated_SNPs = SNP)
SNP = get_SNP(study_design_10000,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design_10000,size, biased = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
aa = get_aa(study_design_10000,size, biased = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "full", beta = 0.5, associated_SNPs = SNP)

SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design,size, biased = "full", stratified = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
aa = get_aa(study_design,size, biased = "full", stratified = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "full", beta = 0.5, associated_SNPs = SNP)
SNP = get_SNP(study_design_10000,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design_10000,size, biased = "full", stratified = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
aa = get_aa(study_design_10000,size, biased = "full", stratified = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "full", beta = 0.5, associated_SNPs = SNP)

####ss4
SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design,size, biased = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "half", beta = 0.25, associated_SNPs = SNP)
aa = get_aa(study_design,size, biased = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "half", beta = 0.5, associated_SNPs = SNP)
SNP = get_SNP(study_design_10000,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design_10000,size, biased = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "half", beta = 0.25, associated_SNPs = SNP)
aa = get_aa(study_design_10000,size, biased = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "half", beta = 0.5, associated_SNPs = SNP)

SNP = get_SNP(study_design,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design,size, biased = "full", stratified = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "half", beta = 0.25, associated_SNPs = SNP)
aa = get_aa(study_design,size, biased = "full", stratified = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "half", beta = 0.5, associated_SNPs = SNP)
SNP = get_SNP(study_design_10000,size, stratified = "full", fst_strat=0.2)
aa = get_aa(study_design_10000,size, biased = "full", stratified = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "half", beta = 0.25, associated_SNPs = SNP)
aa = get_aa(study_design_10000,size, biased = "full", stratified = "full", fst_bias=0.2, associated_strains = "half", associated_populations = "half", beta = 0.5, associated_SNPs = SNP)

######SC1
SNP = get_SNP(study_design,size, stratified = c("P1","P2"), fst_strat=0.2)

aa = get_aa(study_design,size, biased = c("P1","P2"), fst_bias=0.2, stratified = c("A","B"), fst_strat=0.02)
aa = get_aa(study_design,size, biased = c("P1","P2"), fst_bias=0.02, stratified = c("A","B"), fst_strat=0.2)
aa = get_aa(study_design,size, biased = c("P1","P2"), fst_bias=0.2, stratified = c("A","B"), fst_strat=0.2)
####sc1b
SNP = get_SNP(study_design,size, stratified = c("P1","P2"), fst_strat=0.2)

aa = get_aa(study_design,size, biased = c("P1","P2"), fst_bias=0.2, stratified = c("A","B"), fst_strat=0.02, associated_strains = "A", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
aa = get_aa(study_design,size, biased = c("P1","P2"), fst_bias=0.02, stratified = c("A","B"), fst_strat=0.2, associated_strains = "A", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
aa = get_aa(study_design,size, biased = c("P1","P2"), fst_bias=0.2, stratified = c("A","B"), fst_strat=0.2, associated_strains = "A", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)

####sc2, to set new imp
SNP = get_SNP(study_design,size, stratified = c("P1","P2"), fst_strat=0.2)
aa = get_aa(study_design,size, stratified = c("A","B"), fst_strat=0.2, associated_strains = "A", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
aa = get_aa(study_design,size, stratified = c("A","B"), fst_strat=0.2, associated_strains = "A", associated_populations = "full", beta = 0.5, associated_SNPs = SNP)
SNP = get_SNP(study_design_10000,size, stratified = c("P1","P2"), fst_strat=0.2)
aa = get_aa(study_design_10000,size, stratified = c("A","B"), fst_strat=0.2, associated_strains = "A", associated_populations = "full", beta = 0.25, associated_SNPs = SNP)
aa = get_aa(study_design_10000,size, stratified = c("A","B"), fst_strat=0.2, associated_strains = "A", associated_populations = "full", beta = 0.5, associated_SNPs = SNP)