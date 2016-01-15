####Function
createSNPsample <- function(size, pCases, pControls) {
  randomizedCaseProb = abs(c(rnorm(1, mean = pCases[1], sd = 0.1),   rnorm(1, mean = pCases[2], sd = 0.1)))
  randomizedContProb = abs(c(rnorm(1, mean = pControls[1], sd = 0.1),   rnorm(1, mean = pControls[2], sd = 0.01)))
  Cases = rbinom(size=2,n = size,prob = randomizedCaseProb)
  Cont = rbinom(size=2,n = size,prob = randomizedContProb)
  return(c(Cases, Cont))
}

createFullSNPSample <- function(popSize, causalSNPnb, neutralSNPnb) {
  #Generating population genotypes
  causalSNP = replicate(causalSNPnb, createSNPsample(popSize, c(0.4,0.5), c(0.1,0)))
  neutralSNP = replicate(neutralSNPnb, createSNPsample(popSize, c(0.1,0), c(0.1,0)))
  colnames(causalSNP) <-  paste0("SNPc", 1:causalSNPnb)
  colnames(neutralSNP) <-   paste0("SNPn", 1:neutralSNPnb)
  sample = cbind(causalSNP, neutralSNP)
  rownames(sample) <- c(paste0("Case", 1:popSize), paste0("Control", 1:popSize))
  return(sample)
}

computeBeta <- function(causalSNPnb, neutralSNPnb) {
  true_beta1 = replicate(causalSNPnb, rnorm(1,mean = 1,sd = 1))
  true_beta2 = rep(0.00001,neutralSNPnb)
  beta = as.matrix(c(true_beta1,true_beta2))
  colnames(beta) <- 'beta'
  return(beta)
}

computeResultsSet <- function(predictors, response) {
  resulset = matrix(ncol = 4, nrow = 0)
  for(i in 1:ncol(predictors)) {
    glm.out = glm(response~predictors[,i], family=binomial(logit))
    summary = summary(glm.out)
    ##Extract intercept p value
    int_pval=coef(summary)[,4][1]
    ##Extract coeff p-value p value
    coef_pval=signif(coef(summary)[,4][2], digits = 3)
    computed_beta=signif(coef(glm.out)[2], digits = 5)
    initial_beta=signif(beta[i], digits = 5)
    effect_size=exp(coef(glm.out)[2])
    resulset <- rbind(resulset, c(computed_beta,coef_pval, initial_beta, effect_size))
  }
  resulset = cbind(resulset, as.matrix(p.adjust(resulset[,2], method = 'bonferroni')))
  resulset = cbind(resulset, as.matrix(p.adjust(resulset[,2], method = 'BH')))
  colnames(resulset) <- c('computed_beta','p_val','initial_beta', 'effect_size', 'p_adj_BF', 'p_adj_BH')
  rownames(resulset) <- colnames(predictors)
  return(resulset)
}


####Script
popSize = 100
causalSNPnb = 1
neutralSNPnb = 499
SNPdata <- createFullSNPSample(popSize, causalSNPnb, neutralSNPnb)
beta <- computeBeta(causalSNPnb, neutralSNPnb)
z = SNPdata%*%beta
pr = 1/(1+exp(-z))

####Hypothesis testing
#This represent how could the sample be for an already fixed genotype and case/control status;
#According to the beta of the SNPs we randomly generated
# H0: beta = 0 the SNP is not associated with the outcome (case/control), it cannot differentiate case from control.
# H1: beta1 diff 0, SNP1 is associated with the outcome, it can differentiate case from control.
#Type I error: 0.05% chance to reject H0 when it is true
#With pval: On the 500 tested SNPs, there is an average risk of having 5% false positive among all asayed SNPs (0.05 * 500) 25 false positive. FDR
#With padj: On the 500 tested SNPs, with bon ferroni correction there is an average risk of 5% to have a false positive. FWER
#With padj: On the 500 tested SNPs, with BH correction there is an average number of 5% false positive among positive. FWER

#Effect size = for every 1 allele more, (score in SNP1 (causal) ) the odds for having the disease increase by exp(beta)


simulated_initial_case_control = rbinom(200,1,pr)

r1 = computeResultsSet(response = simulated_initial_case_control, predictors = SNPdata)

####Noise effect
# N : What if you add some noise to the data? How do these probabilities get affected? 
# Check if it is related to the logit function shape, where variation are sronger in the midle, if so, what is the effect on the noise ?
noise = rnorm(200, mean = 0,sd = 0.1)
noise2 = rnorm(200, mean = 0,sd = 0.2)
prnoised = 1/(1+exp(-z+noise))
prnoised2 = 1/(1+exp(-z+noise2))
par(mfrow=c(2,2))
curve(1/(1+exp(-x)), 0, 7)
plot(sort(prnoised), type='p', main='w/ noisy data 0.1', ylab='p')
plot(sort(prnoised2), type='p', main='w/ noisy data 0.2', ylab='p')
plot(sort(pr), type='p', main='w/o noisy data', ylab='p')

simulated_initial_case_control_noise1 = rbinom(200,1,prnoised)
simulated_initial_case_control_noise2 = rbinom(200,1,prnoised2)

r2 = computeResultsSet(response = simulated_initial_case_control_noise1, predictors = SNPdata)
r3 = computeResultsSet(response = simulated_initial_case_control_noise2, predictors = SNPdata)


plot(-log(r1[,2]), main='w/o noisy data MP', ylab='-Log(p)', xlab='SNPs', ylim=c(0,20))
plot(-log(r2[,2]), main='w noisy data MP', ylab='-Log(p)', xlab='SNPs', ylim=c(0,20))
plot(-log(r3[,2]), main='w noisier data MP', ylab='-Log(p)', xlab='SNPs', ylim=c(0,20))

qqplot(r1[,2], rnorm(200, mean = 0.95, sd = 0.2), plot.it = TRUE)
qqline(r1[,2])

sum(as.vector(r1[,2]), na.rm = TRUE)
sd(r1[,2], na.rm = TRUE)
summary(r1[,2])[4]


qqplot(r1[,2], rnorm(200, mean = summary(r1[,2])[4], sd = sd(r1[,2], na.rm = TRUE)), plot.it = TRUE)
