#rm(list = ls())
cat("\014")  
library(rjson)
library(parallel)
library(dplyr)
library(ggplot2)
library(reshape2)
require(gridExtra)
#library(compiler)
#library(pryr)
#library(microbenchmark)


source("generate-SNP-datas.R")
source("generate-viral-datas.R")
#source("analyse-generated-datas.R")
#source("system-tools.R")

#enableJIT(3)

#######With virus
trace <- FALSE

###Input
threshold = 5e-05
nb_pop = 2
nb_strain = 2
sample_size = 5000
repnb = 1000

get_SNP_input <- function(S_Stratified = NA, S_Biased = NA,Y_Stratified = NA, Y_Biased = NA, Associated_Strains = NA, Associated_Populations = NA, size = size, Ordered_Bias = FALSE) {
  lvl = c("full","no","half", "P1", "P2", "A", "B")
  if(is.na(Associated_Strains) & !is.na(Associated_Populations)) Associated_Strains = "full"
  if(is.na(Associated_Populations) & !is.na(Associated_Strains)) Associated_Populations = "full"
  #  if(is.na(Associated_Strains) & !is.na(Associated_Populations)) throw("Define associated_Populations")
  # if(is.na(Associated_Populations) & !is.na(Associated_Strains)) throw("Define Associated_Strains")
  tbl_df(data.frame(`S_Stratified` = factor(rep(S_Stratified, size), levels = lvl), `S_Biased` = factor(S_Biased, levels = lvl), `Y_Stratified` = factor(Y_Stratified, levels = lvl), 
                    `Y_Biased` = factor(Y_Biased, levels = lvl), `Associated_Strains` = Associated_Strains, `Associated_Populations` = Associated_Populations, `Ordered_Bias` = Ordered_Bias,
                    `size` = size, stringsAsFactors = FALSE) )}

##NB : When half associated, there is stratification
#fs1 = get_SNP_input(S_Stratified = "full",Y_Biased = "full", size = repnb)
#fs2_ss1 = get_SNP_input(S_Stratified = "full", Associated_Strains = "full", Associated_Populations = "full", size = repnb)
#ss2 = get_SNP_input(S_Stratified = "full", Associated_Strains = "full", Associated_Populations = "half", size = repnb)
#fs3_ss3 = get_SNP_input(S_Stratified = "full", Associated_Strains = "A", Associated_Populations = "full", size = repnb)
#fs4 = get_SNP_input(S_Stratified = "full", Y_Biased = "full",  Associated_Strains = "A", Associated_Populations = "full", size = repnb)
#ss4 = get_SNP_input(S_Stratified = "full", Associated_Strains = "A", Associated_Populations = "half", size = repnb)
#c1 = get_SNP_input(S_Biased = "half",Y_Stratified = "half", size = repnb)
#sc1 = get_SNP_input(S_Stratified = "full", Y_Stratified = "full",Y_Biased = "full", size = repnb)
#sc1b = get_SNP_input(S_Stratified = "full", Y_Stratified = "full",Y_Biased = "full",  Associated_Strains = "A", Associated_Populations = "full", size = repnb)


#SNPs = rbind_all(list(prototype,fs1, fs2_ss1, ss2, fs3_ss3, fs4, ss4, c1, sc1, sc1b))
#res = scenario_viral(sample_size, nb_strain, nb_pop, SNPs, fcoeff_pop = 0.3, fcoeff_vir = 0.1, vir_bias = 0.01, pop_bias = 0.03, beta = beta)

########C1
c1 =  get_SNP_input(S_Biased = "half",Y_Stratified = "half", size = repnb, Ordered_Bias = TRUE)
res = scenario_viral(sample_size, nb_strain, nb_pop, c1, pop_bias = 0.03, fcoeff_vir = 0.2)

c1_eq =  get_SNP_input(S_Biased = "half",Y_Stratified = "half", size = repnb, Ordered_Bias = TRUE)
res = scenario_viral(sample_size, nb_strain, nb_pop, c1_eq, pop_bias = 0.2, fcoeff_vir = 0.2)

c1_without_human_bias = get_SNP_input(S_Stratified = "full",Y_Stratified = "half", size = repnb, Ordered_Bias = TRUE)
res = scenario_viral(sample_size, nb_strain, nb_pop, c1_without_human_bias, fcoeff_pop = 0.2, fcoeff_vir = 0.2)

c1_with_full_bias_full_strat = get_SNP_input(S_Stratified = "full", S_Biased = "full", Y_Stratified = "full", Y_Biased = "full", size = repnb, Ordered_Bias = TRUE)
res = scenario_viral(sample_size, nb_strain, nb_pop, c1_with_full_bias_full_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.03, pop_bias = 0.03)

c1_with_full_bias_full_strat_eq = get_SNP_input(S_Stratified = "full", S_Biased = "full", Y_Stratified = "full", Y_Biased = "full", size = repnb, Ordered_Bias = TRUE)
res = scenario_viral(sample_size, nb_strain, nb_pop, c1_with_full_bias_full_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.2, pop_bias = 0.2)


########FS1
fs1 = get_SNP_input(S_Stratified = "full",Y_Biased = "full", size = repnb, Ordered_Bias = TRUE)
res = scenario_viral(sample_size, nb_strain, nb_pop, fs1, fcoeff_pop = 0.2, vir_bias = 0.02)

fs1_full_strat = get_SNP_input(S_Stratified = "full",Y_Stratified = "full", size = repnb)
res = scenario_viral(sample_size, nb_strain, nb_pop, fs1_full_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.2)

######fs2_ss1 
fs2_ss1 = get_SNP_input(S_Stratified = "full", Associated_Strains = "full", Associated_Populations = "full", size = repnb)
res = scenario_viral(sample_size, nb_strain, nb_pop, fs2_ss1, fcoeff_pop = 0.2, beta = 0.1, tag = "beta_0.1")
res = scenario_viral(sample_size, nb_strain, nb_pop, fs2_ss1, fcoeff_pop = 0.2, beta = 0.5, tag = "beta_0.5")
res = scenario_viral(sample_size, nb_strain, nb_pop, fs2_ss1, fcoeff_pop = 0.2, beta = 0.25, tag = "beta_0.25")

fs2_ss1_without_strat = get_SNP_input(Associated_Strains = "full", Associated_Populations = "full", size = repnb)
res = scenario_viral(sample_size, nb_strain, nb_pop, fs2_ss1_without_strat, beta = 0.1, tag = "beta_0.1")
res = scenario_viral(sample_size, nb_strain, nb_pop, fs2_ss1_without_strat, beta = 0.5, tag = "beta_0.5")
res = scenario_viral(sample_size, nb_strain, nb_pop, fs2_ss1_without_strat, beta = 0.25, tag = "beta_0.25")
######


######ss2
ss2 = get_SNP_input(S_Stratified = "full", Associated_Strains = "full", Associated_Populations = "half", size = repnb)
res = scenario_viral(sample_size = 5000, nb_strain, nb_pop, ss2, fcoeff_pop = 0.2, beta = 0.25, tag = "sample_size_5000")
res = scenario_viral(sample_size = 10000, nb_strain, nb_pop, ss2, fcoeff_pop = 0.2, beta = 0.25, tag = "sample_size_10000")

ss2_with_vir_strat = get_SNP_input(S_Stratified = "full",Y_Stratified = "full", Associated_Strains = "full", Associated_Populations = "half", size = repnb)
res = scenario_viral(sample_size = 5000, nb_strain, nb_pop, ss2_with_vir_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.05, beta = 0.25, tag = "sample_size_5000")
res = scenario_viral(sample_size = 10000, nb_strain, nb_pop, ss2_with_vir_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.05, beta = 0.25, tag = "sample_size_10000")

ss2_with_vir_strat_asymetric = get_SNP_input(S_Stratified = "full",Y_Stratified = "full", Associated_Strains = "full", Associated_Populations = "P1", size = repnb)
res = scenario_viral(sample_size = 10000, nb_strain, nb_pop, ss2_with_vir_strat_asymetric, fcoeff_pop = 0.2, fcoeff_vir = 0.05, beta = 0.25)


ss2_with_vir_strat_eq = get_SNP_input(S_Stratified = "full",Y_Stratified = "full", Associated_Strains = "full", Associated_Populations = "half", size = repnb)
res = scenario_viral(sample_size = 5000, nb_strain, nb_pop, ss2_with_vir_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.25, tag = "sample_size_5000")
res = scenario_viral(sample_size = 10000, nb_strain, nb_pop, ss2_with_vir_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.25, tag = "sample_size_10000")

ss2_with_vir_strat_asymetric_eq = get_SNP_input(S_Stratified = "full",Y_Stratified = "full", Associated_Strains = "full", Associated_Populations = "P1", size = repnb)
res = scenario_viral(sample_size = 10000, nb_strain, nb_pop, ss2_with_vir_strat_asymetric_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.25)


####fs3_ss3
fs3_ss3 = get_SNP_input(S_Stratified = "full", Associated_Strains = "half", Associated_Populations = "full", size = repnb)

res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, fs3_ss3, fcoeff_pop = 0.2, beta = 0.25, tag = "sample_size_5000_beta_0.25")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, fs3_ss3, fcoeff_pop = 0.2, beta = 0.25, tag = "sample_size_10000_beta_5")
res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, fs3_ss3, fcoeff_pop = 0.2, beta = 0.5, tag = "sample_size_5000_beta_0.5")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, fs3_ss3, fcoeff_pop = 0.2, beta = 0.5, tag = "sample_size_10000_beta_0.5")


fs3_ss3_strat = get_SNP_input(S_Stratified = "full",Y_Stratified = "full", Associated_Strains = "half", Associated_Populations = "full", size = repnb)

res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, fs3_ss3_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.05, beta = 0.25, tag = "sample_size_5000")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, fs3_ss3_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.05, beta = 0.25, tag = "sample_size_10000")
res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, fs3_ss3_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.05, beta = 0.5, tag = "sample_size_5000_beta_0.5")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, fs3_ss3_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.05, beta = 0.5, tag = "sample_size_10000_beta_0.5")

fs3_ss3_strat_eq = get_SNP_input(S_Stratified = "full",Y_Stratified = "full", Associated_Strains = "half", Associated_Populations = "full", size = repnb)

res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, fs3_ss3_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.25, tag = "sample_size_5000")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, fs3_ss3_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.25, tag = "sample_size_10000")
res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, fs3_ss3_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.5, tag = "sample_size_5000_beta_0.5")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, fs3_ss3_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.5, tag = "sample_size_10000_beta_0.5")


#####fs4
fs4 = get_SNP_input(S_Stratified = "full",Y_Biased = "full", Associated_Strains = "half", Associated_Populations = "full", size = repnb)

res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, fs4, fcoeff_pop = 0.2, vir_bias = 0.2, beta = 0.25, tag = "sample_size_5000_beta_0.25")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, fs4, fcoeff_pop = 0.2, vir_bias = 0.2, beta = 0.25, tag = "sample_size_10000_beta_0.25")
res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, fs4, fcoeff_pop = 0.2, vir_bias = 0.2, beta = 0.5, tag = "sample_size_5000_beta_0.5")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, fs4, fcoeff_pop = 0.2, vir_bias = 0.2, beta = 0.5, tag = "sample_size_10000_beta_0.5")

fs4_strat = get_SNP_input(S_Stratified = "full",Y_Stratified = "full",Y_Biased = "full", Associated_Strains = "half", Associated_Populations = "full", size = repnb)

res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, fs4_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.05, vir_bias = 0.05, beta = 0.25, tag = "sample_size_5000_beta_0.25")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, fs4_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.05, vir_bias = 0.05, beta = 0.25, tag = "sample_size_10000_beta_0.25")
res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, fs4_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.05, vir_bias = 0.2, beta = 0.5, tag = "sample_size_5000_beta_0.5")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, fs4_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.05, vir_bias = 0.2, beta = 0.5, tag = "sample_size_10000_beta_0.5")

fs4_strat_eq = get_SNP_input(S_Stratified = "full",Y_Stratified = "full",Y_Biased = "full", Associated_Strains = "half", Associated_Populations = "full", size = repnb)

res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, fs4_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.05, beta = 0.25, tag = "sample_size_5000_beta_0.25")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, fs4_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.05, beta = 0.25, tag = "sample_size_10000_beta_0.25")
res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, fs4_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.2, beta = 0.5, tag = "sample_size_5000_beta_0.5")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, fs4_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.2, beta = 0.5, tag = "sample_size_10000_beta_0.5")


####fs4 unbiased
#fs4_unbiased = get_SNP_input(S_Stratified = "full",Y_Biased = "full", Associated_Strains = "half", Associated_Populations = "full", Ordered_Bias = FALSE, size = repnb)

#res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, fs4_unbiased, fcoeff_pop = 0.2, vir_bias = 0.2, beta = 0.25)
#res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, fs4_unbiased, fcoeff_pop = 0.2, vir_bias = 0.2, beta = 0.25)

####ss4
ss4 = get_SNP_input(S_Stratified = "full", Associated_Strains = "half", Associated_Populations = "half", size = repnb)

res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, ss4, fcoeff_pop = 0.2, beta = 0.25, tag = "sample_size_5000_beta_0.25")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, ss4, fcoeff_pop = 0.2, beta = 0.25, tag = "sample_size_10000_beta_0.25")
res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, ss4, fcoeff_pop = 0.2, beta = 0.5, tag = "sample_size_5000_beta_0.5")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, ss4, fcoeff_pop = 0.2, beta = 0.5, tag = "sample_size_10000_beta_0.5")

ss4_with_strat = get_SNP_input(S_Stratified = "full",Y_Stratified = "full", Associated_Strains = "half", Associated_Populations = "half", size = repnb)

res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, ss4_with_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.05, beta = 0.25, tag = "sample_size_5000_beta_0.25")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, ss4_with_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.05, beta = 0.25, tag = "sample_size_10000_beta_0.25")
res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, ss4_with_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.05, beta = 0.5, tag = "sample_size_5000_beta_0.5")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, ss4_with_strat, fcoeff_pop = 0.2, fcoeff_vir = 0.05, beta = 0.5, tag = "sample_size_10000_beta_0.5")

ss4_with_strat_eq = get_SNP_input(S_Stratified = "full",Y_Stratified = "full", Associated_Strains = "half", Associated_Populations = "half", size = repnb)

res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, ss4_with_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.25, tag = "sample_size_5000_beta_0.25")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, ss4_with_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.25, tag = "sample_size_10000_beta_0.25")
res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, ss4_with_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.5, tag = "sample_size_5000_beta_0.5")
res = scenario_viral(sample_size= 10000, nb_strain, nb_pop, ss4_with_strat_eq, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.5, tag = "sample_size_10000_beta_0.5")


######SC1
sc1 = get_SNP_input(S_Stratified = "full", Y_Stratified = "full",Y_Biased = "full", size = repnb, Ordered_Bias = TRUE)
res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, sc1, fcoeff_pop = 0.2, fcoeff_vir = 0.02, vir_bias = 0.2, tag = "*")
#with vir>pop
res = scenario_viral(sample_size = 5000, nb_strain, nb_pop, sc1, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.02, tag = "with_vir_pop")

res = scenario_viral(sample_size = 5000, nb_strain, nb_pop, sc1, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.2, tag = "_eq")


sc1_inv = get_SNP_input(S_Stratified = "P2", Y_Stratified = "full",Y_Biased = "full", size = repnb)
res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, sc1_inv, fcoeff_pop = 0.2, fcoeff_vir = 0.02, vir_bias = 0.2, tag = "*")


####sc1b
sc1b = get_SNP_input(S_Stratified = "full", Y_Stratified = "full", Y_Biased = "full",  Associated_Strains = "A", Associated_Populations = "full", size = repnb, Ordered_Bias = TRUE)
res = scenario_viral(sample_size = 5000, nb_strain, nb_pop, sc1b, fcoeff_pop = 0.2, fcoeff_vir = 0.02, vir_bias = 0.2, beta = 0.25)
res = scenario_viral(sample_size = 5000, nb_strain, nb_pop, sc1b, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.02, beta = 0.25, tag = "with_vir_pop")
res = scenario_viral(sample_size = 5000, nb_strain, nb_pop, sc1b, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.2, beta = 0.25, tag = "_eq")


####sc2
sc2 = get_SNP_input(S_Stratified = "full", Y_Stratified = "full",  Associated_Strains = "A", Associated_Populations = "full", size = repnb, Ordered_Bias = TRUE)
res = scenario_viral(sample_size= 5000,  c(2, 0.35, 0.65), nb_pop, sc2, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.25, tag = "sample_size_5000_beta_0.25")
res = scenario_viral(sample_size= 10000,  c(2, 0.35, 0.65), nb_pop, sc2, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.25, tag = "sample_size_10000_beta_0.25")
res = scenario_viral(sample_size= 5000,  c(2, 0.35, 0.65), nb_pop, sc2, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.5, tag = "sample_size_5000_beta_0.5")
res = scenario_viral(sample_size= 10000,  c(2, 0.35, 0.65), nb_pop, sc2, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.5, tag = "sample_size_10000_beta_0.5")

####Demonstration
demo_no_kill = get_SNP_input(S_Stratified = "full", size = repnb)
res = scenario_viral(sample_size= 5000,  nb_strain, nb_pop, demo_no_kill, fcoeff_pop = 0.2)

demo_power_gain3 = get_SNP_input(S_Stratified = "P1", Associated_Populations = "P1", size = 100, Ordered_Bias = FALSE)
###0.01,0.2
res = scenario_viral(sample_size= 5000,  nb_strain, nb_pop, demo_power_gain3, fcoeff_pop = 0.2, vir_bias = 0.2,beta = 0.25, get_viral = get_viral_output_sp)

####sc2_ret
sc2_ret = get_SNP_input(S_Stratified = "full", Y_Stratified = "full",  Associated_Strains = "A", Associated_Populations = "full", size = 100, Ordered_Bias = TRUE)
res = scenario_viral(sample_size= 5000, nb_strains = c(2, 0.35, 0.65), nb_pop, sc2_ret, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.25, tag = "sample_size_5000_beta_0.25", get_viral = get_viral_output_sp2)
res = scenario_viral(sample_size= 10000,  nb_strains = c(2, 0.35, 0.65), nb_pop, sc2_ret, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.25, tag = "sample_size_10000_beta_0.25", get_viral = get_viral_output_sp2)
res = scenario_viral(sample_size= 5000,  nb_strains = c(2, 0.35, 0.65), nb_pop, sc2_ret, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.5, tag = "sample_size_5000_beta_0.5", get_viral = get_viral_output_sp2)
res = scenario_viral(sample_size= 10000,  nb_strains = c(2, 0.35, 0.65), nb_pop, sc2_ret, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.5, tag = "sample_size_10000_beta_0.5", get_viral = get_viral_output_sp2)

res = scenario_viral(sample_size= 10000,  nb_strains = c(2, 0.35, 0.65), nb_pop, sc2_ret, fcoeff_pop = 0.2, fcoeff_vir = 0.2, beta = 0.5, tag = "sp1_sample_size_10000_beta_0.5", get_viral = get_viral_output_sp)

sc2_simple_implem = get_SNP_input(S_Stratified = "P1", Y_Stratified = "A", Y_Biased = "hB",  Associated_Strains = "A", Associated_Populations = "full", size = 100)
res = scenario_viral(sample_size= 10000,  nb_strains = c(2, 0.35, 0.65), nb_pop, sc2_simple_implem, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.05, beta = 0.25, tag = "sample_size_10000_beta_0.5", get_viral = get_viral_output)

################Experimenting on parallel

cl = makeCluster(30, type = "FORK")
#parLapply(cl, seq(200, 5000, length.out = 10), function(sample_size) {
#  res = res = scenario_viral(sample_size, nb_strain, nb_pop, causal_S =  causal_S, causal_NS = causal_NS , viral_aa = viral_aa,beta =  beta, fcoeff_pop =  fcoeff_pop, vir_bias = vir_bias) })

#size = seq(200, 5000, length.out = 10)
#fcoeff_pop = seq(0.01, 0.3, length.out = 10)
#beta = seq(0.1, 2, length.out = 10)


sc2_simple_implem = get_SNP_input(S_Stratified = "P1", Y_Stratified = "A", Y_Biased = "hB",  Associated_Strains = "A", Associated_Populations = "full", size = 100)
res = res = scenario_viral(sample_size= 100,  nb_strains = c(2, 0.35, 0.65), nb_pop, sc2_simple_implem, fcoeff_pop = 0.2, fcoeff_vir = 0.2, vir_bias = 0.05, beta = 0.25)


res = res = scenario_viral(sample_size= 5000, nb_strain, nb_pop, fs3_ss3, fcoeff_pop = 0.2, beta = 0.5, tag = "sample_size_5000_beta_0.5")


scale_y_continuous(limits=c(0,100), breaks=seq(0,100,10), expand = c(0, 0))

parLapply(beta, seq(0.1, 2, length.out = 30), function(sample_size) {
  res = scenario_viral(sample_size, nb_strain, nb_pop, SNPs, fcoeff_pop, fcoeff_vir, vir_bias, pop_bias, beta) })

res = parLapply(cl, seq(0.01, 0.2, length.out = 30), function(param) {
  test = select(as.data.frame(scenario_viral(sample_size, nb_strain, nb_pop, fs1_full_strat, fcoeff_pop = param, fcoeff_vir = param, tag = "fs1_strat_impact")$pval ), ends_with("pval"))
  medianpval = apply(pval, 2, function(p) median(p))
  rm(test,pval)
  c(param, medianpval)
})

res = parLapply(cl, seq(0.01, 1, length.out = 100), function(param) {
  res = scenario_viral(sample_size, nb_strain, nb_pop, fs1, fcoeff_pop = param, vir_bias = param, tag = "fs1_strat_impact_fp")
  test = unlist(res$summary_sim["FP_ratio",])----------------------------------------
    rm(res)
  c(param, test)
})


test_fp = unlist(res$summary_sim["FP_ratio",])

write.table(res,file = "result_on_scenario_fs1_bias")


res = scenario_viral(sample_size, nb_strain, nb_pop, fs1, fcoeff_pop = param, vir_bias = param, tag = "fs1_strat_impact_fp")

parLapply(beta, seq(0.1, 2, length.out = 30), function(sample_size) {
  scenario_viral(sample_size, nb_strain, nb_pop, SNPs, fcoeff_pop, fcoeff_vir, vir_bias, pop_bias, beta) })

parLapply(cl, seq(200, 5000, length.out = 10), function(sample_size) {
  for(fcoeff_pop in seq(0.01, 0.3, length.out = 4)) {
    for(fcoeff_vir in seq(0.01, 0.3, length.out = 4)) {
      for(beta in seq(0.1, 2, length.out = 4)) {
        for(pop_bias in seq(0.01, 0.3, length.out = 4)) {
          for(vir_bias in seq(0.01, 0.3, length.out = 4)) {
            res = scenario_viral(sample_size, nb_strain, nb_pop, SNPs, fcoeff_pop, fcoeff_vir, vir_bias, pop_bias, beta) }}}}}})