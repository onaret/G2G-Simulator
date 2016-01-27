rm(list = ls())
library(ggbiplot)

####Function
get_alternate_allele <- function(reference_alleles) {
  sapply(reference_alleles, function(reference_allele) {
    s1 = reference_allele*(1-fcoeff)/fcoeff
    s2 = (1-reference_allele)*(1-fcoeff)/fcoeff
    rbeta(n = 1, shape1 = s1,shape2 = s2)
  })
}

get_genotype_probability <- function(alternate_allele, R) {
  if(R!=0) {
    p_associated = c((1-alternate_allele)^2, 2*R*alternate_allele*(1-alternate_allele), R^2 * alternate_allele^2 )
    p_associated/(sum(p_associated) + (1-alternate_allele^2))
  }
  else c((1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2)
}

get_associated_probability <- function(alternate_allele, R) {
  p_associated = c((1-alternate_allele)^2, 2*R*alternate_allele*(1-alternate_allele), R^2 * alternate_allele^2 )
  scaled = sum(p_associated) + (1-alternate_allele^2)
  unname(p_associated/scaled)
}

generate_SNPs <- function(SNP_struct, populations) {
  a = c()
  ###Generating datas
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
  
  ###Formatting matrix
  row_names = c(unlist(sapply(1:length(populations), function(sample) {
    names = c()
    if(populations[[sample]]$case > 0) names = c(paste('Case', sample, 1:populations[[sample]]$case,sep='_'))
    if(populations[[sample]]$control > 0) names = c(names, paste('Control', sample, 1:populations[[sample]]$control,sep='_'))
    names
  })))
  
  generated_SNP = matrix(data, ncol = length(SNP_struct$SNP_association), nrow = length(row_names), dimnames = list(row_names, SNP_struct$SNP_names))
  as.matrix(generated_SNP[ order(rownames(generated_SNP)), ])
}

scenario <- function(populations, SNP_association, strat_percentage) {
  if(strat_percentage>1) strat_percentage = strat_percentage/100
  if(strat_percentage>100) strat_percentage = 1
  sort(SNP_association)
  stratified_SNP = min(trunc(strat_percentage * length(SNP_association)), table(SNP_association)["0"])

  #Generating stratified allele frequency
  if(stratified_SNP>0) {
    populations = lapply(1:length(populations), function(population){
      if(population %% 2 == 0) {
        RP = replicate(stratified_SNP, runif(1,0.7,0.9))
      }
      else {
        RP = replicate(stratified_SNP, runif(1,0.1,0.3))
      }
      populations[[population]]['alternate_allele'] = list(get_alternate_allele(RP))
      populations[[population]]
    })
  }


  #Generating unstratified allele frequency
  RP = replicate(length(SNP_association) - stratified_SNP, runif(1,0.4,0.5))
  populations = lapply(1:length(populations), function(population){
    populations[[population]]['alternate_allele'] = mapply(c, populations[[population]]['alternate_allele'], list(get_alternate_allele(RP)), SIMPLIFY=FALSE)
    populations[[population]]
  })
  SNP_names = c(sapply(1:stratified_SNP, function(snp_num){c(paste0('Snp_S_',snp_num,"_R",SNP_association[snp_num]))}), 
                sapply((stratified_SNP+1):length(SNP_association), function(snp_num){c(paste0('Snp_NS_',snp_num,"_R",SNP_association[snp_num]))}))
  SNP_struct = list(SNP_association=SNP_association, SNP_names=SNP_names)
  generated_SNP = generate_SNPs(SNP_struct = SNP_struct, populations = populations)
  list(generated_SNP = generated_SNP, population_structure = populations)
}

####Constants
C1 = list(P1 = c(case = 200, control = 400), P2 = c(case = 400, control = 200))
C2 = list(P1 = c(case = 400, control = 200), P2 = c(case = 200, control = 400))
C3 = list(P1 = c(case = 300, control = 0), P2 = c(case = 300, control = 600))
C4 = list(P1 = c(case = 300, control = 200), P2 = c(case = 200, control = 100), P3 = c(case = 100, control = 300))
C5 = list(P1 = c(case = 200, control = 0), P2 = c(case = 400, control = 200), P3 = c(case = 0, control = 400))

fcoeff  = 0.01 ##### Wright's coefficient for inbreeding

######Scenarios
SNPdata = scenario(populations=C1, SNP_association = c(0,0,0,4), strat_percentage = 50)

#Result_S1 = scenario(scenario="nonstrat", populations=C1, SNP_struct = c(0,0,0))
#Result_S2 = scenario(scenario="strat", populations=C1, SNP_struct = c(0,0,0))
#Result_S3 = scenario(scenario="nonstrat", populations=C1, SNP_struct = 4)
#Result_S4 = scenario(scenario="strat", populations=C1, SNP_struct = 4)

#Result_C1 = scenario(scenario="strat", populations=C1, SNP_struct = c(0,0,0,5))
#Result_C2 = scenario(scenario="nonstrat", populations=C1, SNP_struct = c(0,0,0,5))

#Result_C2 = scenario(scenario="nonstrat", populations=list(c(case = 200, control = 0), c(case = 400, control = 200), c(case = 0, control = 400), c(case=500, control=200)), SNP_struct = c(0,0,0,0,0))

##Alongside R, pass a Vector indicating if this SNP should be stratified

######Tests
##### SNP data with SNPs on the columns and samples on the rows. 
#SNPdata = cbind(Result_S1$generated_SNP, Result_S2$generated_SNP, Result_S3$generated_SNP, Result_S4$generated_SNP)
#SNPdata = cbind(Result_S2$generated_SNP, Result_S1$generated_SNP)
#SNPdata = cbind(Result_C1$generated_SNP, Result_C2$generated_SNP)

#sum(Result_S2[[1]][1:600,6])
#sum(Result_S2[[1]][600:1200,6])

#### Generate the pc vector for population group identification. 
PC_simulated = c(rep("A",200),rep("B",400),rep("A",400),rep("B",200))
y = c(rep(1,600),rep(0,600))   ###### Generate the response vector Y. This will have 600 ones and 600 zeros.

##########WHY
SNPdata = SNPdata[[1]]
summary(glm(y~SNPdata, family=binomial(logit)))        ###### Model without pc

sum(SNPdata[1:600,1])
sum(SNPdata[600:1200,1])

sum(SNPdata[1:600,3])
sum(SNPdata[600:1200,3])

sum(SNPdata[1:600,4])
sum(SNPdata[600:1200,4])
#############

ca = prcomp(SNPdata, scale. = TRUE) 
#print(ca)
#summary(ca)
PC_computed = ca$x[,1]<0

summary(glm(y~SNPdata, family=binomial(logit)))        ###### Model without pc
summary(glm(y~SNPdata+PC_computed))
summary(glm(y~SNPdata+PC_simulated))

ggbiplot(prcomp(SNPdata, scale. = TRUE))
ggbiplot(prcomp(SNPdata), obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = TRUE)
