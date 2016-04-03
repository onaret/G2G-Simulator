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

SNPstrat = get_SNP(study_design,5000, stratified = "full", fst_strat=0.1)
SNPunstrat = get_SNP(study_design,90000)
aaunstrat = get_aa(study_design,1000)


scenario_full(study_design, aa1, aa2, aa3, aa4, aa5, aa6, aa7, aa8, SNP1, aa9, aa10, aa11, SNP2, aa12, aa13, aa14, SNP3, aa15, aa16, aa17, SNP4, aa18, aa19, aa20, SNP5, aa21, aa22, aa23, SNP6, aa24, aa25, aa26, SNP7, aa27, aa28, aa29, SNP8, aa30, aa31, aa32, SNP9, aa33, aa34, aa35, SNP10, aa36, aa37, aa38, aa39, aa40, aa41, SNP12, aa42, aa43, aa44, SNPstrat, SNPunstrat, aaunstrat)

scenario_full(study_design, aa7, aa8, SNP1, aa9, aa10, aa11, SNP2, aa12, aa13, aa14)


####sc2, to set new imp
SNP13 = get_SNP(study_design,size, stratified = c("P1","P2"), fst_strat=0.2)
aa45 = get_aa(study_design,size, stratified = c("A","B"), fst_strat=0.2, associated_strains = "A", associated_populations = "full", beta = 0.25, associated_SNPs = SNP13)
aa46 = get_aa(study_design,size, stratified = c("A","B"), fst_strat=0.2, associated_strains = "A", associated_populations = "full", beta = 0.5, associated_SNPs = SNP13)
SNP14 = get_SNP(study_design_10000,size, stratified = c("P1","P2"), fst_strat=0.2)
aa47 = get_aa(study_design_10000,size, stratified = c("A","B"), fst_strat=0.2, associated_strains = "A", associated_populations = "full", beta = 0.25, associated_SNPs = SNP14)
aa48 = get_aa(study_design_10000,size, stratified = c("A","B"), fst_strat=0.2, associated_strains = "A", associated_populations = "full", beta = 0.5, associated_SNPs = SNP14)