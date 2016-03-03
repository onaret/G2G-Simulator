rm(list = ls())
cat("\014")  
library(rjson)
library(parallel)
library(dplyr)
#library(compiler)
#library(pryr)
#library(microbenchmark)
#library(ggplot2)

source("generate-SNP-datas.R")
source("generate-viral-datas.R")
#source("analyse-generated-datas.R")
#source("system-tools.R")

#enableJIT(3)

#######With virus

###Input
sequence = 1
threshold = 5e-07
#threshold = 0.05
pop_nb = 2
vir_pop = c("A","B")
fcoeff_pop  = 0.3 ##### Wright's coefficient for inbreeding
fcoeff_vir  = 0.3
vir_bias = 0.02
beta = 0.1
sample_size = 5000
repnb = 100

get_SNP_input <- function(S_Stratified = NA, S_Biased = NA,Y_Stratified = NA, Y_Biased = NA, Associated_Strains = NA, Associated_Populations = NA, size) {
  lvl = c("full","no","half")
  prototype = data.frame(`S_Stratified` = factor(levels = lvl), `S_Biased` = factor(levels = lvl), `Y_Stratified` = factor(levels =lvl), `Y_Biased` = factor(levels = lvl), `Associated_Strains` = character(), `Associated_Populations` = factor(levels = lvl))
  input = data.frame(`S_Stratified` = rep(S_Stratified,size), `S_Biased` = S_Biased, `Y_Stratified` = Y_Stratified, `Y_Biased` = Y_Biased, `Associated_Strains` = Associated_Strains, `Associated_Populations` = Associated_Populations)
  rbind_all(list(prototype,input))
  #input = data.frame(`S_Stratified` = factor(S_Stratified, levels = lvl), `S_Biased` = factor(S_Biased, levels = lvl), `Y_Stratified` = factor(Y_Stratified, levels = lvl), `Y_Biased` = factor(Y_Biased, levels = lvl), `Associated_Strains` = Associated_Strains, `Associated_Populations` = factor(Associated_Populations, levels = lvl)) 
  #as.data.frame(t(replicate(size,input, simplify = "data.frame")))
}

##NB : When half associated, there is stratification
fs1 = get_SNP_input(S_Stratified = "full",Y_Biased = "full", size = repnb)
fs2_ss1 = get_SNP_input(S_Stratified = "full", Associated_Strains = "AB", Associated_Populations = "full", size = repnb)
ss2 = get_SNP_input(S_Stratified = "full", Associated_Strains = "AB", Associated_Populations = "half", size = repnb)
fs3_ss3 = get_SNP_input(S_Stratified = "full", Associated_Strains = "A", Associated_Populations = "full", size = repnb)
fs4 = get_SNP_input(S_Stratified = "full", Y_Biased = "full",  Associated_Strains = "A", Associated_Populations = "full", size = repnb)
ss4 = get_SNP_input(S_Stratified = "full", Associated_Strains = "A", Associated_Populations = "half", size = repnb)
c1 = get_SNP_input(S_Biased = "half",Y_Stratified = "half", size = repnb)
sc1 = get_SNP_input(S_Stratified = "full", Y_Stratified = "full",Y_Biased = "full", size = repnb)
sc1b = get_SNP_input(S_Stratified = "full", Y_Stratified = "full",Y_Biased = "full",  Associated_Strains = "A", Associated_Populations = "full", size = repnb)


SNPs = rbind_all(list(prototype,fs1, fs2_ss1, ss2, fs3_ss3, fs4, ss4, c1, sc1, sc1b))
scenario_viral(sample_size, vir_pop, pop_nb, SNPs, fcoeff_pop = 0.3, fcoeff_vir = 0.1, vir_bias = 0.01, pop_bias = 0.03, beta = beta)


########C1
c1 =  get_SNP_input(S_Biased = "half",Y_Stratified = "half", size = repnb)
scenario_viral(sample_size, vir_pop, pop_nb, c1, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.03, pop_bias = 0.03, beta = beta)

c1_without_human_bias = get_SNP_input(S_Stratified = "full",Y_Stratified = "half", size = repnb)
scenario_viral(sample_size, vir_pop, pop_nb, c1_without_human_bias, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.03, pop_bias = 0.03, beta = beta)

c1_without_human_bias = get_SNP_input(S_Stratified = "full",Y_Stratified = "half", size = repnb)
scenario_viral(sample_size, vir_pop, pop_nb, c1_without_human_bias, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.03, pop_bias = 0.03, beta = beta)

c1_with_full_bias_full_strat = get_SNP_input(S_Stratified = "full", S_Biased = "full", Y_Stratified = "full", Y_Biased = "full", size = repnb)
scenario_viral(sample_size, vir_pop, pop_nb, c1_with_full_bias_full_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.03, pop_bias = 0.03, beta = beta)


########FS1
fs1 = get_SNP_input(S_Stratified = "full",Y_Biased = "full", size = repnb)
scenario_viral(sample_size, vir_pop, pop_nb, fs1, fcoeff_pop = 0.2, fcoeff_vir = 0.01, vir_bias = 0.003, pop_bias = 0.03, beta = beta)

fs1_full_strat = get_SNP_input(S_Stratified = "full",Y_Stratified = "full", size = repnb)
scenario_viral(sample_size, vir_pop, pop_nb, fs1_full_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.001, vir_bias = 0.001, pop_bias = 0.03, beta = beta)

cl = makeCluster(10, type = "FORK")
#parLapply(cl, seq(200, 5000, length.out = 10), function(sample_size) {
#  res = scenario_viral(sample_size, vir_pop, pop_nb, causal_S =  causal_S, causal_NS = causal_NS , viral_aa = viral_aa,beta =  beta, fcoeff_pop =  fcoeff_pop, vir_bias = vir_bias) })

#size = seq(200, 5000, length.out = 10)
#fcoeff_pop = seq(0.01, 0.3, length.out = 10)
#beta = seq(0.1, 2, length.out = 10)

parLapply(cl, seq(200, 5000, length.out = 10), function(sample_size) {
  for(fcoeff_pop in seq(0.01, 0.3, length.out = 4)) {
    for(fcoeff_vir in seq(0.01, 0.3, length.out = 4)) {
      for(beta in seq(0.1, 2, length.out = 4)) {
        for(pop_bias in seq(0.01, 0.3, length.out = 4)) {
          for(vir_bias in seq(0.01, 0.3, length.out = 4)) {
            scenario_viral(sample_size, vir_pop, pop_nb, SNPs, fcoeff_pop, fcoeff_vir, vir_bias, pop_bias, beta) }}}}}})