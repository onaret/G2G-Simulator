#rm(list = ls())
library(ggplot2)
setwd("~/workspace")
cat("\014")  

####Function

scenario <- function(populations, SNP_association, stratification_rate) {
  if(stratification_rate>1) stratification_rate = stratification_rate/100
  if(stratification_rate>100) stratification_rate = 1
  sort(SNP_association)
  stratified_SNP = min(trunc(stratification_rate * length(SNP_association)), table(SNP_association)["0"])
  
  SNP_names = c()
  print("Generating stratified allele frequency")
  allele_frequencies = if(stratified_SNP>0) {
    SNP_names = c(sapply(1:stratified_SNP, function(snp_num)c(paste0('Snp_S_',snp_num,"_R",SNP_association[snp_num]))))
    sapply(1:length(populations), function(population){
      if(population %% 2 == 0) RP = replicate(stratified_SNP, runif(1,0.7,0.9))
      else RP = replicate(stratified_SNP, runif(1,0.1,0.3))
      get_alternate_allele(RP)
    })
  }
  
  print("Generating unstratified allele frequency")
  allele_frequencies = rbind(allele_frequencies, if((length(SNP_association) - stratified_SNP) >0) {
    SNP_names = c(SNP_names, sapply((stratified_SNP+1):length(SNP_association), function(snp_num) c(paste0('Snp_NS_',snp_num,"_R",SNP_association[snp_num]))))
    RP = replicate(length(SNP_association) - stratified_SNP, runif(1,0.4,0.5))
    sapply(1:length(populations), function(population) get_alternate_allele(RP))
  })
  
  populations = mapply(function(population, pos) {
    population$alternate_allele = allele_frequencies[,pos]
    list(population)}, populations, 1:length(populations))
  
  SNP_struct = list(SNP_association=SNP_association, SNP_names=SNP_names)
  samples_SNPs = generate_SNPs(SNP_struct = SNP_struct, populations = populations)
  
  cases_nb = sum(sapply(populations, function(population) population$case))
  controls_nb = sum(sapply(populations, function(population) population$control))
  study_infos = list(`cases_nb`= cases_nb, controls_nb = `controls_nb`, `SNP_names`= SNP_names)
  list(`samples_SNPs` = samples_SNPs, `population_structure` = populations, `study_infos` = study_infos)
}

get_alternate_allele <- function(reference_alleles) {
  sapply(reference_alleles, function(reference_allele) {
    s1 = reference_allele*(1-fcoeff)/fcoeff
    s2 = (1-reference_allele)*(1-fcoeff)/fcoeff
    rbeta(n = 1, shape1 = s1,shape2 = s2)
  })
}

generate_SNPs <- function(SNP_struct, populations) {
  print("Generating datas")
  data = sapply(1:length(SNP_struct$SNP_association), function(snp_num){
    data = sapply(1:length(populations), function(sample) {
      pCase = get_genotype_probability(populations[[sample]]$alternate_allele[snp_num], SNP_struct$SNP_association[snp_num])
      pControl = get_genotype_probability(populations[[sample]]$alternate_allele[snp_num], 0)
      cases=sample(c(0,1,2), size = populations[[sample]]$case, prob = pCase, replace = TRUE)
      controls=sample(c(0,1,2), size = populations[[sample]]$control, prob = pControl, replace = TRUE)
      c(cases, controls)
    })
    unlist(data)
  })
  
  print("Formatting matrix")
  row_names = c(unlist(sapply(1:length(populations), function(sample) {
    if(populations[[sample]]$case > 0) row_names = c(paste('Case', sample, 1:populations[[sample]]$case,sep='_'))
    if(populations[[sample]]$control > 0) row_names = c(row_names, paste('Control', sample, 1:populations[[sample]]$control,sep='_'))
    row_names
  })))
  
  generated_SNP = matrix(data, ncol = length(SNP_struct$SNP_association), nrow = length(row_names), dimnames = list(row_names, SNP_struct$SNP_names))
  as.matrix(generated_SNP[ order(rownames(generated_SNP)), ])
}

get_genotype_probability <- function(alternate_allele, R) {
  if(R!=0) {
    p_associated = c((1-alternate_allele)^2, 2*R*alternate_allele*(1-alternate_allele), R^2 * alternate_allele^2 )
    p_associated/(sum(p_associated) + (1-alternate_allele^2))
  }
  else c((1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2)
}

analyse <- function(data) {
  simulated_PC = unlist(c(mapply(function(population, length) {rep(length, population$case)}, data$population_structure, 1:length(data$population_structure))
                          , mapply(function(population, length) {rep(length, population$control)}, data$population_structure, 1:length(data$population_structure))))
  
  y = c(rep(1,data$study_infos$cases_nb),rep(0,data$study_infos$controls_nb))
  
  print("Computing PCA")
  ca = prcomp(data$samples_SNPs, scale. = FALSE) 
  print("Computing logistic regression with computed PC")
  summary_glm_w_computed_PC = sapply(1:ncol(data$samples_SNPs), function(SNP) coef(summary(glm(y~data$samples_SNPs[,SNP]+ca$x[,1:5])))[,4][2])
  print("Computing logistic regression with simulated  PC")
  summary_glm_w_simulated_PC = sapply(1:ncol(data$samples_SNPs), function(SNP) coef(summary(glm(y~data$samples_SNPs[,SNP]+simulated_PC)))[,4][2])
  print("Computing logistic regression without PC")
  summary_glm_wo_PC = sapply(1:ncol(data$samples_SNPs), function(SNP) coef(summary(glm(y~data$samples_SNPs[,SNP])))[,4][2])
  print("Done")
  data.frame(summary_glm_wo_PC, summary_glm_w_computed_PC, summary_glm_w_simulated_PC, row.names = data$study_infos$SNP_names)
}

####Constants
C1 = list(`P1` = list(`case` = 201, `control` = 401), `P2` = list(`case` = 399, `control`  = 199))
C2 = list(`P1` = c(`case` = 400, `control` = 200), `P2` = c(`case` = 200, `control` = 400))
C3 = list(`P1` = c(`case` = 300, `control` = 0), `P2` = c(`case` = 300, `control` = 600))
C4 = list(`P1` = c(`case` = 300, `control` = 200), `P2` = c(`case` = 200, `control` = 100), `P3` = c(`case` = 100, `control` = 300))
C5 = list(`P1` = c(`case` = 200, `control` = 0), `P2` = c(`case` = 400, `control` = 200), `P3` = c(`case` = 0, `control` = 400))

fcoeff  = 0.01 ##### Wright's coefficient for inbreeding

######Scenarios
SNPdataNoAssoc = scenario(populations=C1, SNP_association = c(rep(0,100000), 1,2,3,4,5,6,7,8,9), stratification_rate = 0.05)
SNPdataNoAssocAnalyze = analyse(SNPdataNoAssoc)
#summarise(SNPdataNoAssoc)
#SNPdataXXLNoAssoc = scenario(populations=C1, SNP_association = c(rep(0,100000)), stratification_rate = 0.05)
#SNPdataXXLNoAssocAnalyzed = analyse(SNPdataXXLNoAssoc)

par(mfrow=c(3,2))
qqnorm(SNPdataNoAssocAnalyze[,1],ylab="Standardized Residuals", xlab="Normal Scores", main="without_PC")
qqnorm(SNPdataNoAssocAnalyze[,2],ylab="Standardized Residuals", xlab="Normal Scores", main="with_computed_PC")
qqnorm(SNPdataNoAssocAnalyze[,3],ylab="Standardized Residuals", xlab="Normal Scores", main="with_simulated_PC")

plot(-log(SNPdataNoAssocAnalyze[,1]), main='without_PC', ylab='-Log(p)', xlab='SNPs', ylim=c(0,20))
plot(-log(SNPdataNoAssocAnalyze[,2]), main='with_computed_PC', ylab='-Log(p)', xlab='SNPs', ylim=c(0,20))
plot(-log(SNPdataNoAssocAnalyze[,3]), main='with_simulated_PC', ylab='-Log(p)', xlab='SNPs', ylim=c(0,20))

#SNPdataAssoc = scenario(populations=C1, SNP_association = c(0,0,0,5), stratification_rate = 0.5)
#SNPdataXXL = scenario(populations=C1, SNP_association = c(rep(0,50),4), stratification_rate = 0.05)