nsamp = 500
nvar = 5000
strat_snp = (10/100)*nvar
asso_snp = (10/100)*strat_snp
SNPdata = matrix(0,nsamp,nvar)
pop1_samp = nsamp/2
pop2_samp = nsamp/2

######## Assign samples to viral strains

virG = sample(c("A","B"),nsamp,replace=T)

####### Assign which ones will have y = 1 and which one will have y = 0. They can come from any viral strain randomly.

CC = virG
###In groupe A, half will be randomly 1
CC[sample(which(virG == "A"), size = (sum(virG == "A"))/2)] = 1
CC[sample(which(virG == "B"), size = (sum(virG == "B"))/2)] = 1
CC[which(CC == "A")] = 0
CC[which(CC == "B")] = 0


s = numeric(nsamp)

RP1 = runif(1,0.7,0.9) ##### ancestral allele frequency for population 1
RP2 = runif(1,0.1,0.3) ##### ancestral allele frequency for population 2

RP = c(RP1,RP2)

fcoeff  = 0.01 ##### Wright's coefficient for inbreeding

##### Generate allele frequencies for 2 different populations

AP = numeric(length(RP))

for(i in 1:length(RP)) {
  s1 = (RP[i]*(1-fcoeff))/fcoeff
  s2 = ((1-RP[i])*(1-fcoeff))/fcoeff
  AP[i] = rbeta(n = 1, shape1 = s1,shape2 = s2)
}

pCase1 = c((1-AP[1])^2,(2*AP[1])*(1-AP[1]),(AP[1]^2))      ###### cases from pop1
sCase1 = sample(which(CC==1),size = sum(CC==1)/2)
s[sCase1] = sample(0:2,size = length(sCase1),prob = pCase1,replace=T)

pCase2 = c((1-AP[2])^2,(2*AP[2])*(1-AP[2]),(AP[2]^2))      ###### cases from pop2
sCase2 = setdiff(which(CC==1),sCase1) 
s[sCase2] = sample(0:2,size = length(sCase2),prob = pCase2,replace=T)

pCont1 = c((1-AP[1])^2,(2*AP[1])*(1-AP[1]),(AP[1]^2))      ###### controls from pop1
sCont1 = sample(which(CC==0),size = sum(CC==0)/2)
s[sCont1] = sample(0:2,size = length(sCont1),prob = pCont1,replace=T)

pCont2 = c((1-AP[2])^2,(2*AP[2])*(1-AP[2]),(AP[2]^2))      ###### controls from pop1
sCont2 = setdiff(which(CC==0),sCont1)
s[sCont2] = sample(0:2,size = length(sCont2),prob = pCont2,replace=T)


############## Generate Y randomly in such a way that one population has more viral mutations than the other.

Y = numeric(nsamp)

for(i in 1:length(sCase1)){Y[sCase1[i]] = sample(0:1,1,prob = c(0.3,0.7))}

for(i in 1:length(sCase2)){Y[sCase2[i]] = sample(0:1,1,prob = c(0.7,0.3))}

for(i in 1:length(sCont1)){Y[sCont1[i]] = sample(0:1,1,prob = c(0.3,0.7))}

for(i in 1:length(sCont2)){Y[sCont2[i]] =sample(0:1,1,prob = c(0.7,0.3))}


pcs = numeric(nsamp)
pcs[sCase1] = "G1"
pcs[sCont1] = "G1"
pcs[pcs!="G1"] = "G2"

summary(glm(Y~s))
summary(glm(Y~s+pcs))
summary(glm(Y~s+pcs+virG))
