rm(list = ls())


generate_sample <- function() {
  
}

##Enum Category to pass as argument
##Custom categor with nbCase/nbContrl per population, data struct
####population = c(nbCase, nbControl)
####population 3 is optionanl

###### Generation for Category 1, 2 sub-populations

###### Generate catefory first and pass it as arguments to a generateSample Function.
###### Arguments, differentiation according to population structure, association with the disease.
###### If there is association wih the disease define r.


##### Generate allele frequencies for 2 different populations

generate_allele_frequencies <- function(RP) {
  AP = numeric(length(RP))  
  for(i in 1:length(RP)) {
    s1 = (RP[i]*(1-fcoeff))/fcoeff
    s2 = ((1-RP[i])*(1-fcoeff))/fcoeff
    AP[i] = rbeta(n = 1, shape1 = s1,shape2 = s2)
  }
  return(matrix(ncol = length(RP), data = c(RP,AP)))
}

generate_genotypes <- function(AF, population_struct, associated_with_disease) {
  if (associated_with_disease == FALSE) {
    for(i in seq(1, length(AF), by = 2) ) {
      pCase1 = c(AF[i]^2, (2*AF[i])*(1-AF[i+1]))   ##### for 0 the prob is (1-p2[1])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 
      Cases1 = rbinom(size=2,n = population_struct[i],prob = pCase1)
      ##### Cases from population 2
      pCase2 = c(AP[2]^2, (2*AP[2])*(1-AP[2]))   ##### for 0 the prob is (1-p2[2])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 
      Cases2 = rbinom(size=2,n = 400,prob = pCase2)
    }
  }
  else {
    
  }
}

####### Scenario 1: Random SNP with no differentiation and no association with the disease

RP1 = runif(1,0.4,0.5) ##### ancestral allele frequency for population 1
RP2 = RP1 ##### ancestral allele frequency for population 2

fcoeff  = 0.01 ##### Wright's coefficient for inbreeding

population_struct = c(200,400,400,200)

AF = generate_allele_frequencies(c(RP1,RP2))

generate_genotypes(AF, population_struct, associated_with_disease = TRUE)

######## Generate the genotypes for 200 Cases from pop1 and 400 Cases from pop2

##### Cases from population 1

pCase1 = c(AP[1]^2, (2*AP[1])*(1-AP[1]))   ##### for 0 the prob is (1-p2[1])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 

Cases1 = rbinom(size=2,n = 200,prob = pCase1)

##### Cases from population 2

pCase2 = c(AP[2]^2, (2*AP[2])*(1-AP[2]))   ##### for 0 the prob is (1-p2[2])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 

Cases2 = rbinom(size=2,n = 400,prob = pCase2)


######## Generate the genotypes for 400 Controls from pop1 and 200 Controls from pop2

##### Controls from population 1

pCont1 = c(AP[1]^2, (2*AP[1])*(1-AP[1]))   ##### for 0 the prob is (1-p2[1])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 

Cont1 = rbinom(size=2,n = 400,prob = pCont1)

##### Controls from population 2

pCont2 = c(AP[2]^2, (2*AP[2])*(1-AP[2]))   ##### for 0 the prob is (1-p2[2])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 

Cont2 = rbinom(size=2,n = 200,prob = pCont2)


###### SNP data

SNP1 = c(Cases1,Cases2,Cont1,Cont2)


#############################################################################################################################

####### Scenario 2: Random SNP with differentiation and no association with the disease

RP1 = runif(1,0.7,0.9) ##### ancestral allele frequency for population 1
RP2 = runif(1,0.1,0.3) ##### ancestral allele frequency for population 2

RP = c(RP1,RP2)

fcoeff  = 0.01 ##### Wright's coefficient for inbreeding

##### Generate allele frequencies for 2 different populations

AP = numeric(length(RP))

for(i in 1:length(RP)) {
     s1 = (RP[i]*(1-fcoeff))/fcoeff
     s2 = ((1-AP[i])*(1-fcoeff))/fcoeff
     AP[i] = rbeta(n = 1, shape1 = s1,shape2 = s2)
}

######## Generate the genotypes for 200 Cases from pop1 and 400 Cases from pop2

##### Cases from population 1

pCase1 = c(AP[1]^2, (2*AP[1])*(1-AP[1]))   ##### for 0 the prob is (1-p2[1])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 

Cases1 = rbinom(size=2,n = 200,prob = pCase1)

##### Cases from population 2

pCase2 = c(AP[2]^2, (2*AP[2])*(1-AP[2]))   ##### for 0 the prob is (1-p2[2])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 

Cases2 = rbinom(size=2,n = 400,prob = pCase2)


######## Generate the genotypes for 400 Controls from pop1 and 200 Controls from pop2

##### Controls from population 1

pCont1 = c(AP[1]^2, (2*AP[1])*(1-AP[1]))   ##### for 0 the prob is (1-p2[1])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 

Cont1 = rbinom(size=2,n = 400,prob = pCont1)

##### Controls from population 2

pCont2 = c(AP[2]^2, (2*AP[2])*(1-AP[2]))   ##### for 0 the prob is (1-p2[2])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 

Cont2 = rbinom(size=2,n = 200,prob = pCont2)


###### SNP data

SNP2 = c(Cases1,Cases2,Cont1,Cont2)

#############################################################################################################################

####### Scenario 3: Random SNP without differentiation and association with the disease

RP1 = runif(1,0.4,0.5) ##### ancestral allele frequency for population 1
RP2 = RP1 ##### ancestral allele frequency for population 2

RP = c(RP1,RP2)

fcoeff  = 0.01 ##### Wright's coefficient for inbreeding

##### Generate allele frequencies for 2 different populations

AP = numeric(length(RP))

for(i in 1:length(RP)) {
     s1 = (RP[i]*(1-fcoeff))/fcoeff
     s2 = ((1-AP[i])*(1-fcoeff))/fcoeff
     AP[i] = rbeta(n = 1, shape1 = s1,shape2 = s2)
}

######## Generate the genotypes for 200 Cases from pop1 and 400 Cases from pop2

#### For association let R = 1.5. R is the relative risk

##### Cases from population 1
R = 4

pCase1 = c((R^2 * AP[1]^2), (2*R*AP[1])*(1-AP[1]))   ##### for 0 the prob is (1-p2[1])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 
scaled = sum(pCase1) + (1-AP[1])^2
pCase1 = pCase1/scaled
Cases1 = rbinom(size=2,n = 200,prob = pCase1)

##### Cases from population 2

pCase2 = c((R^2 * AP[2]^2), (2*R*AP[2])*(1-AP[2]))   ##### for 0 the prob is (1-p2[2])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 
scaled = sum(pCase2) + (1-AP[2])^2
pCase2 = pCase2/scaled
Cases2 = rbinom(size=2,n = 400,prob = pCase2)


######## Generate the genotypes for 400 Controls from pop1 and 200 Controls from pop2

##### Controls from population 1

pCont1 = c(AP[1]^2, (2*AP[1])*(1-AP[1]))   ##### for 0 the prob is (1-p2[1])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 

Cont1 = rbinom(size=2,n = 400,prob = pCont1)

##### Controls from population 2

pCont2 = c(AP[2]^2, (2*AP[2])*(1-AP[2]))   ##### for 0 the prob is (1-p2[2])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 

Cont2 = rbinom(size=2,n = 200,prob = pCont2)


###### SNP data

SNP3 = c(Cases1,Cases2,Cont1,Cont2)

#############################################################################################################################

####### Scenario 4: Random SNP with differentiation and with association with the disease

RP1 = runif(1,0.7,0.9) ##### ancestral allele frequency for population 1
RP2 = runif(1,0.1,0.3) ##### ancestral allele frequency for population 2

RP = c(RP1,RP2)

fcoeff  = 0.01 ##### Wright's coefficient for inbreeding

##### Generate allele frequencies for 2 different populations

AP = numeric(length(RP))

for(i in 1:length(RP)) {
     s1 = (RP[i]*(1-fcoeff))/fcoeff
     s2 = ((1-AP[i])*(1-fcoeff))/fcoeff
     AP[i] = rbeta(n = 1, shape1 = s1,shape2 = s2)
}

######## Generate the genotypes for 200 Cases from pop1 and 400 Cases from pop2

#### For association let R = 1.5. R is the relative risk

##### Cases from population 1
R = 4

pCase1 = c((R^2 * AP[1]^2), (2*R*AP[1])*(1-AP[1]))   ##### for 0 the prob is (1-p2[1])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 
scaled = sum(pCase1) + (1-AP[1])^2
pCase1 = pCase1/scaled
Cases1 = rbinom(size=2,n = 200,prob = pCase1)

##### Cases from population 2

pCase2 = c((R^2 * AP[2]^2), (2*R*AP[2])*(1-AP[2]))   ##### for 0 the prob is (1-p2[2])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 
scaled = sum(pCase2) + (1-AP[2])^2
pCase2 = pCase2/scaled
Cases2 = rbinom(size=2,n = 400,prob = pCase2)


######## Generate the genotypes for 400 Controls from pop1 and 200 Controls from pop2

##### Controls from population 1

pCont1 = c(AP[1]^2, (2*AP[1])*(1-AP[1]))   ##### for 0 the prob is (1-p2[1])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 

Cont1 = rbinom(size=2,n = 400,prob = pCont1)

##### Controls from population 2

pCont2 = c(AP[2]^2, (2*AP[2])*(1-AP[2]))   ##### for 0 the prob is (1-p2[2])^2 .  The first value in the vector is the prob for 2 and the sencond value is the prob for 1. 

Cont2 = rbinom(size=2,n = 200,prob = pCont2)


###### SNP data

SNP4 = c(Cases1,Cases2,Cont1,Cont2)

###### Test the models for the SNP data 

SNPdata = cbind(SNP1,SNP2,SNP3,SNP4)      ##### SNP data with SNPs on the columns and samples on the rows. 
pc = c(rep("A",200),rep("B",400),rep("A",400),rep("B",200))    #### Generate the pc vector for population group identification. 
y = c(rep(1,600),rep(0,600))   ###### Generate the response vector Y. This will have 600 ones and 600 zeros.
summary(glm(y~SNPdata))        ###### Model without pc
summary(glm(y~SNPdata+pc))     ###### Model with pc



######################################################################################
#### Change the min max of uniform distribution for RP1 and RP2. For example in the fourth scenario RP1 comes from 0.7,0.9 and RP2 from 0.1,0.3. Instead make it RP1 from 0.1,0.3 and RP2 from 0.7 and 0.9. 
##### 
