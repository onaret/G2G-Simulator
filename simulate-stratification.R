rm(list = ls())
cat("\014")  
library(rjson)
library(parallel)
library(plyr)
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
#threshold = 5e-07
threshold = 0.05
pop_nb = 2
vir_pop = c("A","B")
fcoeff_pop  = 0.3 ##### Wright's coefficient for inbreeding
fcoeff_vir  = 0.2
vir_bias = 0.1
beta = 0.1
sample_size = 5000

part_causal_s_and_y =  data.frame(`S_Stratified` = c(rep(TRUE,40), rep(FALSE,40)),`S_Biased` = c(rep(TRUE, 20), rep(FALSE,20)),`Y_Stratified` = c(rep(TRUE,10), rep(FALSE,10)),`Y_Biased` = c(rep(TRUE, 5), rep(FALSE,5)), `Viral_Association` = "A", `R` = c(0.1))
part_causal_s_and_y = filter(part_causal_s_and_y, !(`S_Stratified` == FALSE & `S_Biased` == TRUE))
causal =  data.frame(`S_Stratified` = c(rep(TRUE,5), rep(FALSE,5)), `Viral_Association` = "AB", `R` = c(0.1))
neutral_param_square_2 = data.frame(`S_Stratified` = c(rep(TRUE,40), rep(FALSE,40)),`S_Biased` = c(rep(TRUE, 20), rep(FALSE,20)),`Y_Stratified` = c(rep(TRUE,10), rep(FALSE,10)),`Y_Biased` = c(rep(TRUE, 5), rep(FALSE,5)))
neutral_param_square_2 = filter(neutral_param_square_2, !(`S_Stratified` == FALSE & `S_Biased` == TRUE))

ucausal = unique(causal)
uneutral_param_square_2 = unique(neutral_param_square_2)
ucpart_causal_s_and_y = unique(part_causal_s_and_y)







s1 = data.frame(`S_Stratified` = rep(TRUE, 100),`Y_Biased` = TRUE)
s2 = data.frame(`S_Stratified` = rep(TRUE, 100),  `Viral_Association` = "AB", `R` = c(0.1))
s3 = data.frame(`S_Stratified` = rep(TRUE, 1), `Viral_Association` = "A", `R` = c(0.1))
s4 = data.frame(`S_Stratified` = rep(TRUE, 1),`Y_Biased` = TRUE,  `Viral_Association` = "A", `R` = c(0.1))

c1 = data.frame(`S_Half_Biased` = rep(TRUE, 100),`Y_Half_Stratified` = TRUE)
sc1 = data.frame(`S_Stratified` = rep(TRUE, 100), `Y_Stratified` = TRUE,`Y_Biased` = TRUE)
sc1b = data.frame(`S_Stratified` = rep(TRUE, 100), `Y_Stratified` = TRUE,`Y_Biased` = TRUE,  `Viral_Association` = "A", `R` = c(0.1))





SNPs = rbind.fill(s1,s2,s3,s4,sc1,sc1b)

scenario_viral(sample_size, vir_pop, pop_nb, SNPs = c1, fcoeff_pop = 0.3, fcoeff_vir = 0.1, vir_bias = 0.01, pop_bias = 0.03, beta = beta)

scenario_viral(sample_size, vir_pop, pop_nb, SNPs = sc1b, fcoeff_pop = 0.3, fcoeff_vir = 0.1, vir_bias = 0.01, pop_bias = 0.03, beta = beta)

SC1 = c(`fcoeff_pop`=0.1, `fp_vir` = 0.05)














#sc = data.frame(`S_Stratified` = rep(TRUE, 10),`S_Biased` = TRUE, `Y_Stratified` = TRUE,`Y_Biased` = TRUE, `Viral_Association` = "A", `R` = c(0.1))


s1 = data.frame(`S_Stratified` = rep(TRUE, 100),`S_Biased` = FALSE, `Y_Stratified` = FALSE,`Y_Biased` = TRUE)
s2 = data.frame(`S_Stratified` = rep(TRUE, 100),`S_Biased` = FALSE,  `Viral_Association` = "AB", `R` = c(0.1))
s3 = data.frame(`S_Stratified` = rep(TRUE, 1),`S_Biased` = FALSE, `Y_Stratified` = FALSE,`Y_Biased` = FALSE,  `Viral_Association` = "A", `R` = c(0.1))
s4 = data.frame(`S_Stratified` = rep(TRUE, 1),`S_Biased` = FALSE, `Y_Stratified` = FALSE,`Y_Biased` = TRUE,  `Viral_Association` = "A", `R` = c(0.1))
#TODO c1 = 
sc1 = data.frame(`S_Stratified` = rep(TRUE, 100),`S_Biased` = FALSE, `Y_Stratified` = TRUE,`Y_Biased` = TRUE)
sc1b = data.frame(`S_Stratified` = rep(TRUE, 100),`S_Biased` = FALSE, `Y_Stratified` = TRUE,`Y_Biased` = TRUE,  `Viral_Association` = "A", `R` = c(0.1))

SNPs = rbind.fill(s1,s2,s3,s4,sc1,sc1b)

scenario_viral(sample_size, vir_pop, pop_nb, SNPs = SNPs, fcoeff_pop = 0.3, fcoeff_vir = 0.1, vir_bias = 0.01, pop_bias = 0.03, beta = beta)

scenario_viral(sample_size, vir_pop, pop_nb, SNPs = sc1b, fcoeff_pop = 0.3, fcoeff_vir = 0.1, vir_bias = 0.01, pop_bias = 0.03, beta = beta)

SC1 = c(`fcoeff_pop`=0.1, `fp_vir` = 0.05)




cl = makeCluster(1, type = "FORK")
parLapply(cl, seq(200, 5000, length.out = 10), function(sample_size) {
  res = scenario_viral(sample_size, vir_pop, pop_nb, causal_S =  causal_S, causal_NS = causal_NS , viral_aa = viral_aa,beta =  beta, fcoeff_pop =  fcoeff_pop, vir_bias = vir_bias)
})

size = seq(200, 5000, length.out = 10)
fcoeff_pop = seq(0.01, 0.3, length.out = 10)
beta = seq(0.1, 2, length.out = 10)

res = c()
for(sample_size in seq(200, 5000, length.out = 4)) {
  for(fcoeff_pop in seq(0.01, 0.3, length.out = 4)) {
    for(beta in seq(0.1, 2, length.out = 4)) {
      scenario_viral(sample_size, vir_pop, pop_nb, SNPs = SNPs, fcoeff_pop = 0.3, fcoeff_vir = 0.1, vir_bias = 0.01, pop_bias = 0.03, beta = beta)
      }}}


####Scenario

###1
#scenario_viral(sample_size, vir_pop, pop_nb,  neutral = 1,  neutral_S_rate = 1)

###2
#scenario_viral(sample_size, vir_pop, pop_nb, causal_S = causal_S_2)

#scenario_viral(sample_size = sample_size, vir_pop = vir_pop, pop_nb =  pop_nb,causal_S = causal_S_3)

##New step: 
# - Viral repartition between group, S1 : --> DONE
# - Y generation linked or not with SNPs, S2 : --> DONE
#   -not linked : --> DONE
#   -linked : --> DONE
# - Virus side stratification S3 and S4
# - PC with, SNPs, human strat, viral strat S1-S4 : --> DONE
# - Take a function to generate Y to allow special things : --> DONE
# - Implement multiple AA strain