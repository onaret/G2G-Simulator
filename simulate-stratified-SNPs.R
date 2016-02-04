rm(list = ls())
cat("\014")  
#library(rjson)
setwd("/home/zod/Documents/Workspace/EPFL")
####Function
compute_multiple <- function(populations= NULL, size = NULL, neutral, neutral_stratification_rate, causal_S, causal_NS, times) {
  sapply(1:times, function(time){
    print(paste(Sys.time(),"Executing scenario", time, sep=": "))
    res = scenario(populations, size, neutral, neutral_stratification_rate, causal_S, causal_NS)
    write(res, time)
  })
  print(paste(Sys.time(),"Done", sep=": "))}

scenario <- function(populations= NULL, size = NULL, neutral, neutral_stratification_rate, causal_S, causal_NS) {
  neutral_S = get_stratified_SNPs_qtt(neutral_stratification_rate, neutral)
  neutral_NS = neutral - neutral_S
  populations = generate_population_structure(populations, size)
  SNP_freq = generate_SNPs_frequencies(neutral_S, neutral_NS, causal_S, causal_NS, populations)
  SNPs = generate_SNPs(SNP_freq, populations)
  pvalues = analyse(SNPs, populations)
  SNP_struct = get_SNP_struct(neutral_S, neutral_NS, causal_S, causal_NS)
  summary_sim = summary_sim(pvalues, SNP_struct)
  trace_plot(pvalues)
  
  populations = rbind(populations, matrix(c(sum(populations[,"Case"]), sum(populations[,"Control"])),ncol = 2, nrow = 1, dimnames = list("Total") ) )
  list(#`SNPs` = SNPs,
    `populations` = populations,
    `SNP_freq` = SNP_freq,
    `SNP_struct` = SNP_struct,
    `pvalues` = pvalues,
    `summary_sim` = summary_sim)
}

get_stratified_SNPs_qtt <- function(stratification_rate, neutral_SNPs) {
  if(stratification_rate>1) stratification_rate = stratification_rate/100
  if(stratification_rate>100) stratification_rate = 1
  stratification_rate * neutral_SNPs}

generate_population_structure <- function(populations,size) {
  populations = {
    if(!is.list(populations))   {
      populations = sample(c(populations["min"]:populations["max"]), size=1, prob=rep(1/7,7), replace = TRUE)
      populations = min(size, populations)
      populations = runif(populations,0,1)
      populations = round(populations/sum(populations) * size)
      t(sapply(populations, function(strat_size, pop_id) {
        case_nb = round(strat_size*runif(1,0,1))
        #case_nb = round(strat_size*abs(rnorm(5,0.5,0.2)))
        c(case_nb,(strat_size-case_nb))}))}
    else if(is.null(populations) && !is.null(size)) {
      print("Population is not stratfied")
      case_nb = round(size*runif(1,0,1))
      c(case_nb,(strat_size-case_nb))}
    else if(is.null(populations) && is.null(size)) throw("You must set size or population arguments...")
    else populations}
  colnames(populations) <- c("Case", "Control")
  rownames(populations) <- paste0("P", 1:nrow(populations))
  populations}

generate_SNPs_frequencies <- function(neutral_S, neutral_NS, causal_S, causal_NS, populations) {
  print(paste(Sys.time(),"Generating allele frequency", sep=": "))
  AF = {
    if(neutral_S > 0) AF = sapply(1:nrow(populations), function(population) get_alternate_allele(replicate(neutral_S, runif(1,0.01,0.99))))
    if(neutral_NS > 0) AF = rbind(AF, replicate(nrow(populations), get_alternate_allele(replicate(neutral_NS, runif(1,0.4,0.5)))))
    if(length(causal_S) > 0) AF = rbind(AF, sapply(1:nrow(populations), function(population) get_alternate_allele(replicate(length(causal_S), runif(1,0.01,0.99)))))
    if(length(causal_NS) > 0) rbind(AF, replicate(nrow(populations), get_alternate_allele(replicate(length(causal_NS), runif(1,0.4,0.5)))))}
  
  AF = setNames(data.frame(AF), paste0("P", 1:nrow(populations), "_AF"))
  
  SNP_names = {
    if(neutral_S>0) SNP_names = sapply(1:neutral_S, function(snp_num)c(paste0('Snp_S_',snp_num,"_R0")))
    if(neutral_NS >0) SNP_names = c(SNP_names, sapply(1:neutral_NS, function(snp_num) paste0('Snp_NS_',snp_num,"_R0")))
    if(length(causal_S)>0) SNP_names = c(SNP_names, sapply(1:length(causal_S), function(snp_num)c(paste0('Snp_S_',snp_num,"_R", causal_S[snp_num]))))
    if(length(causal_NS)>0) c(SNP_names, sapply(1:length(causal_NS), function(snp_num)c(paste0('Snp_NS_',snp_num,"_R", causal_NS[snp_num]))))}
  
  data.frame(AF, `R`=c(rep(0,neutral_S + neutral_NS), causal_S, causal_NS), row.names = SNP_names)}

get_alternate_allele <- function(reference_alleles) {
  sapply(reference_alleles, function(reference_allele) {
    s1 = reference_allele*(1-fcoeff)/fcoeff
    s2 = (1-reference_allele)*(1-fcoeff)/fcoeff
    rbeta(n = 1, shape1 = s1,shape2 = s2)})}

generate_SNPs <- function(SNP_freq, populations) {
  print(paste(Sys.time(),"Generating SNPs dsitribution", sep=": "))
  data = unlist(sapply(1:nrow(SNP_freq), function(snp_num){
    sapply(1:nrow(populations), function(population) {
      pCase = get_genotype_probability(SNP_freq[snp_num,population], SNP_freq$R[snp_num])
      pControl = get_genotype_probability(SNP_freq[snp_num,population], 0)
      cases=sample(c(0,1,2), size = populations[population,"Case"], prob = pCase, replace = TRUE)
      controls=sample(c(0,1,2), size = populations[population,"Control"], prob = pControl, replace = TRUE)
      c(cases, controls)})}))
  
  sample_names = unlist(sapply(1:nrow(populations), function(population) {
    if(populations[population,"Case"] > 0) sample_names = c(paste('Case', 'P', population, 'sample', 1:populations[population,"Case"],sep='_'))
    if(populations[population,"Control"] > 0) c(sample_names, paste('Control', 'P', population, 'sample', 1:populations[population,"Control"],sep='_'))}))
  matrix(data, ncol = nrow(SNP_freq), nrow = length(sample_names), dimnames = list(sample_names, rownames(SNP_freq)))}

get_genotype_probability <- function(alternate_allele, R) {
  if(R!=0) {
    p_associated = c((1-alternate_allele)^2, 2*R*alternate_allele*(1-alternate_allele), R^2 * alternate_allele^2 )
    p_associated/(sum(p_associated) + (1-alternate_allele^2))}
  else c((1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2)}

analyse <- function(SNPs, populations) {
  simulated_PC = rep(1:nrow(populations), populations[1:nrow(populations),"Control"] + populations[1:nrow(populations), "Case"])
  simulated_PC = chartr("123456789", "ABCDEFGHI", simulated_PC)
  y = unlist(c(sapply(1:nrow(populations), function(pop_nb) c( rep(1, populations[pop_nb,"Case"]), rep(0, populations[pop_nb,"Control"])))))
  print(paste(Sys.time(),"Computing PCA", sep=": "))
  ca = prcomp(SNPs, scale. = FALSE) 
  
  print(paste(Sys.time(),"Computing logistic regression with computed PC", sep=": "))
  W_computed_PC = sapply(1:ncol(SNPs), function(SNP) coef(summary(glm(y~SNPs[,SNP]+ca$x[,1:5])))[,4][2])
  
  print(paste(Sys.time(),"Computing logistic regression with simulated  PC", sep=": "))
  W_simulated_PC = sapply(1:ncol(SNPs), function(SNP) coef(summary(glm(y~SNPs[,SNP]+simulated_PC)))[,4][2])
  
  print(paste(Sys.time(),"Computing logistic regression without PC", sep=": "))
  WO_PC = sapply(1:ncol(SNPs), function(SNP) coef(summary(glm(y~SNPs[,SNP])))[,4][2])
  
  apply(data.frame(WO_PC, W_simulated_PC, W_computed_PC, row.names = colnames(SNPs)),2 , function(pval) {
    signiff = c()
    signiff[threshold < pval] <- 1
    signiff[threshold*0.5 < pval & pval <= threshold] <- 2
    signiff[threshold*0.05 < pval & pval <= threshold*0.5] <- 3
    signiff[threshold*0.005 < pval & pval <= threshold*0.05] <- 4
    signiff[pval <= threshold*0.005] <- 5
    data.frame(pval, `signifficance` = factor(signiff, level=1:5, label=c("ns","*","**","***","+")))})}

get_SNP_struct <- function(neutral_S, neutral_NS, causal_S, causal_NS) {
  SNP_struct = matrix(ncol = 3, nrow = 4, dimnames = list(c("neutral_S","neutral_NS","causal_S","causal_NS"), c("Quantity", "From", "To")))
  SNP_struct[,"Quantity"] <- c(neutral_S, neutral_NS, length(causal_S), length(causal_NS))
  SNP_struct[,"From"] <- c(1, neutral_S + 1, neutral_S + neutral_NS + 1, neutral_S + neutral_NS + length(causal_S) + 1)
  SNP_struct[,"To"] <- c(neutral_S, neutral_S + neutral_NS, neutral_S + neutral_NS + length(causal_S), neutral_S + neutral_NS + length(causal_S) + length(causal_NS))
  SNP_struct}

summary_sim <- function(pvalues, SNP_struct) {
  sapply(pvalues, function(result) {
    extract_ns_row <- function(res_frag, cond) rownames(res_frag)[res_frag[,"signifficance"] != "ns"]
    FP_neutral_S = extract_ns_row(result[SNP_struct["neutral_S","From"]:SNP_struct["neutral_S","To"],])
    FP_neutral_NS = extract_ns_row(result[SNP_struct["neutral_NS","From"]:SNP_struct["neutral_NS","To"],])
    FN_causal_S = extract_ns_row(result[SNP_struct["causal_S","From"]:SNP_struct["causal_S","To"],])
    FN_causal_NS = extract_ns_row(result[SNP_struct["causal_NS","From"]:SNP_struct["causal_NS","To"],])
    total = c(`tot_FP_neutral_S` = length(FP_neutral_S),`tot_FP_neutral_NS` = length(FP_neutral_NS),`tot_FN_causal_S` = length(FN_causal_S), `tot_FN_causal_NS` = length(FN_causal_NS))
    res = list(`FP_neutral_S` = FP_neutral_S, `FP_neutral_NS` = FP_neutral_NS, `FN_causal_S` = FN_causal_S, `FN_causal_S` = FN_causal_NS, `total`=total)
    Filter(length, res)})}

trace_plot <- function(data, save = FALSE, file) {
  ifelse((save == FALSE), print(paste(Sys.time(),"Ploting...", sep=": ")), print(paste(Sys.time(),"Writting datas", sep=": "))) 
  
  if(save == TRUE) png(file, width = 960, height = 1440)
  par(mfrow=c(3,2))
  qq(data$WO_PC$pval,"without_PC")
  plot(-log(data$WO_PC$pval), main='without_PC', ylab='-Log(p)', xlab='SNPs', ylim=c(0,20))
  abline(h = -log10(threshold), untf = FALSE, col = "red")
  
  qq(data$W_simulated_PC$pval, "with_simulated_PC")
  plot(-log(data$W_simulated_PC$pval), main='with_simulated_PC', ylab='-Log(p)', xlab='SNPs', ylim=c(0,20))
  abline(h = -log10(threshold), untf = FALSE, col = "red")
  
  qq(data$W_computed_PC$pval,"with_computed_PC")
  plot(-log(data$W_computed_PC$pval), main='with_computed_PC', ylab='-Log(p)', xlab='SNPs', ylim=c(0,20))
  abline(h = -log10(threshold), untf = FALSE, col = "red")
  if(save == TRUE) dev.off()}

qq <- function(pvector, title="Quantile-quantile plot of p-values", spartan=F) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=0.25, main = title,xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)),ylim=c(0,max(e)))
  lines(e,e,col="red")
  plot}

write <- function(res, time) {
  end = Sys.time()
  trace_plot(res$pvalues, save = TRUE, file = paste0("gen-data/", end,"-N",time, "pvalues.png"))
  write.csv2(x = res$pvalues, file = paste0("gen-data/", end,"-N",time,"-pvalues.csv"))
  write.csv2(x = res$SNP_freq,file = paste0("gen-data/", end,"-N",time,"-SNP_freq.csv"))
  write.csv2(x = res$SNP_struct,file = paste0("gen-data/", end,"-N",time,"-SNP_struct.csv"))
  write.csv2(x = res$populations,file = paste0("gen-data/", end,"-N",time,"-populations.csv"))
  write.csv2(x = data.frame(`WO_PC`= res$summary_sim$WO_PC$total, `W_simulated_PC` = res$summary_sim$W_simulated_PC$total,`W_computed_PC` = res$summary_sim$W_computed_PC$total),file = paste0("gen-data/", end,"-N",time,"-summary_sim.csv"))
  #cat(toJSON(res$summary_sim), file = paste0("gen-data/", end,"-N",time,"-summary_sim.json"))
}

####Constants
C1 = list(`P1` = list(`case` = 201, `control` = 401), `P2` = list(`case` = 399, `control`  = 199))
C2 = list(`P1` = c(`case` = 400, `control` = 200), `P2` = c(`case` = 200, `control` = 400))
C3 = list(`P1` = c(`case` = 300, `control` = 0), `P2` = c(`case` = 300, `control` = 600))
C4 = list(`P1` = c(`case` = 300, `control` = 200), `P2` = c(`case` = 200, `control` = 100), `P3` = c(`case` = 100, `control` = 300))
C5 = list(`P1` = c(`case` = 200, `control` = 0), `P2` = c(`case` = 400, `control` = 200), `P3` = c(`case` = 0, `control` = 400))

fcoeff  = 0.01 ##### Wright's coefficient for inbreeding
threshold = 3.63*10^-8
######Scenarios
#SNPs_1M = scenario(populations=C1, SNP_association = c(rep(0,100000), seq(1,2, by = 0.1)), stratification_rate = 0.05)
#SNPs_test = scenario(populations = 5, size = 1200, SNP_association = c(rep(0,1200), seq(1,2, by = 0.1)), stratification_rate = 0.05)

#scenario(populations = 5, size = 1200, SNP_association = c(rep(0,1200), seq(1,2, by = 0.1)), stratification_rate = 0.05)
#compute_multiple(populations = 5, size = 1200, SNP_association = c(rep(0,1200), seq(1,2, by = 0.1)), stratification_rate = 0.05, times = 1)
#compute_multiple(populations = 5, size = 1200, SNP_association = c(rep(0,100000), seq(1,2, by = 0.01)), stratification_rate = 0.05, times = 20)

compute_multiple(populations = c(`min`=2, `max`=8),
                 size = 1200, 
                 neutral = 1200, 
                 neutral_stratification_rate = 0.05, 
                 causal_NS = seq(1,2, by = 0.05),
                 causal_S = seq(1,2, by = 0.05),
                 times = 1000)

#compute_multiple(
#  populations = 5, 
#  size = 1200, 
#  neutral_SNPs = 100000, 
#  neutral_stratification_rate = 0.05, 
#  causal_SNPs = seq(1,2, by = 0.01),
# causal_SNPs_stratified = seq(1,2, by = 0.01),
# times = 20)

####TODO
##Compute power difference with and without PC to find causal SNPs in function of R
##Compute lambda, (inflation)
##Optimize GLM
##Generate different stratified populations, compute power differences with and without PC,
#See how well the PC correct the stratification, independently of its complexity
##Gene burden test
##Skat + globaltest
#Simplify
#Compute lambda which is an inflation coefficient