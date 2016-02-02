#rm(list = ls())
cat("\014")  

####Function
compute_multiple <- function(populations= NULL, size = NULL, SNP_association, stratification_rate, times) {
  sapply(1:times, function(time){
    res = scenario(populations, size, SNP_association, stratification_rate)
    trace_plot(res$pvalues, save = TRUE, file = paste0("gen-data/", Sys.time(),"-N",time, "pvalues.png"))
    write.csv2(x = res$pvalues, file = paste0("gen-data/", Sys.time(),"-N",time,"-pvalues.csv"))
    write.csv2(x = res$SNP_struct,file = paste0("gen-data/", Sys.time(),"-N",time,"-study-infos.csv"))
    write.csv2(x = res$populations,file = paste0("gen-data/", Sys.time(),"-N",time,"-populations.csv"))})}

scenario <- function(populations= NULL, size = NULL, SNP_association, stratification_rate) {
  stratified_SNP_qtt = get_stratified_SNPs_qtt(stratification_rate, SNP_association)
  populations = generate_population_structure(populations, size)
  SNP_association = sort(SNP_association)
  SNP_struct = generate_SNPs_frequencies(stratified_SNP_qtt, SNP_association, populations)
  #SNPs = generate_SNPs(SNP_struct, populations)
  pvalues = analyse(generate_SNPs(SNP_struct, populations), populations)
  trace_plot(pvalues)
  populations = rbind(populations, matrix(c(sum(populations[,"Case"]), sum(populations[,"Control"])),ncol = 2, nrow = 1, dimnames = list("Total") ) )
  #print(paste(Sys.time(),"Done", sep=": "))
  #summary_sim = summary_sim(results,stratified_SNP_qtt,SNP_association)
  
  list(#`SNPs` = SNPs,
    `populations` = populations,
    `SNP_struct` = SNP_struct,
    `pvalues` = pvalues)}

get_stratified_SNPs_qtt <- function(stratification_rate, SNP_association) {
  if(stratification_rate>1) stratification_rate = stratification_rate/100
  if(stratification_rate>100) stratification_rate = 1
  min(trunc(stratification_rate * length(SNP_association)), table(SNP_association)["0"])}

generate_population_structure <- function(populations,size) {
  populations = {
    if(is.atomic(populations))   {
      populations = min(size, populations)
      populations=runif(populations,0,1)
      populations = round(populations/sum(populations) * size)
      t(sapply(populations, function(strat_size, pop_id) {
        case_nb = round(strat_size*runif(1,0,1))
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

generate_SNPs_frequencies <- function(stratified_SNP_qtt, SNP_association, populations) {
  print(paste(Sys.time(),"Generating allele frequency", sep=": "))
  AF = if(stratified_SNP_qtt>0) {
    sapply(1:nrow(populations), function(population){
      if(population %% 2 == 0) RP = replicate(stratified_SNP_qtt, runif(1,0.7,0.9))
      else RP = replicate(stratified_SNP_qtt, runif(1,0.1,0.3))
      get_alternate_allele(RP)})}
  
  AF = rbind(AF, if((length(SNP_association) - stratified_SNP_qtt) >0) {
    RP = replicate(length(SNP_association) - stratified_SNP_qtt, runif(1,0.4,0.5))
    replicate(nrow(populations), get_alternate_allele(RP))})
  AF = setNames(data.frame(AF), paste0("P", 1:nrow(populations), "_AF"))
  
  SNP_names = {
    if(stratified_SNP_qtt>0) SNP_names = sapply(1:stratified_SNP_qtt, function(snp_num)c(paste0('Snp_S_',snp_num,"_R",SNP_association[snp_num])))
    if((length(SNP_association) - stratified_SNP_qtt) >0) SNP_names = c(SNP_names, sapply((stratified_SNP_qtt+1):length(SNP_association), function(snp_num) paste0('Snp_NS_',snp_num,"_R",SNP_association[snp_num])))}
  
  data.frame(AF, SNP_association, row.names = SNP_names)}

get_alternate_allele <- function(reference_alleles) {
  sapply(reference_alleles, function(reference_allele) {
    s1 = reference_allele*(1-fcoeff)/fcoeff
    s2 = (1-reference_allele)*(1-fcoeff)/fcoeff
    rbeta(n = 1, shape1 = s1,shape2 = s2)})}

generate_SNPs <- function(SNP_struct, populations) {
  print(paste(Sys.time(),"Generating SNPs dsitribution", sep=": "))
  data = sapply(1:length(SNP_struct$SNP_association), function(snp_num){
    unlist(sapply(1:nrow(populations), function(sample) {
      pCase = get_genotype_probability(SNP_struct[snp_num,sample], SNP_struct$SNP_association[snp_num])
      pControl = get_genotype_probability(SNP_struct[snp_num,sample], 0)
      cases=sample(c(0,1,2), size = populations[sample,"Case"], prob = pCase, replace = TRUE)
      controls=sample(c(0,1,2), size = populations[sample,"Control"], prob = pControl, replace = TRUE)
      c(cases, controls)}))})
  
  sample_names = unlist(sapply(1:nrow(populations), function(sample) {
    if(populations[sample,"Case"] > 0) sample_names = c(paste('Case', sample, 1:populations[sample,"Case"],sep='_'))
    if(populations[sample,"Control"] > 0) c(sample_names, paste('Control', sample, 1:populations[sample,"Control"],sep='_'))}))
  matrix(data, ncol = nrow(SNP_struct), nrow = length(sample_names), dimnames = list(sample_names, rownames(SNP_struct)))}

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
  
  apply(data.frame(W_computed_PC, W_simulated_PC, WO_PC, row.names = colnames(SNPs)),2 , function(pval) {
    signiff = c()
    signiff[threshold < pval] <- 1
    signiff[threshold*0.5 < pval & pval <= threshold] <- 2
    signiff[threshold*0.05 < pval & pval <= threshold*0.5] <- 3
    signiff[threshold*0.005 < pval & pval <= threshold*0.05] <- 4
    signiff[pval <= threshold*0.005] <- 5
    data.frame(pval, `signifficance` = factor(signiff, level=1:5, label=c("ns","*","**","***","+")))})}


trace_plot <- function(data, save = FALSE, file) {
  print(paste(Sys.time(),"Ploting...", sep=": "))
  
  if(save == TRUE) png(file)
  par(mfrow=c(3,2))
  qq(data$WO_PC$pval,"without_PC")
  qq(data$W_simulated_PC$pval, "with_simulated_PC")
  qq(data$W_computed_PC$pval,"with_computed_PC")
  
  plot(-log(data$WO_PC$pval), main='without_PC', ylab='-Log(p)', xlab='SNPs', ylim=c(0,20))
  abline(h = -log10(threshold), untf = FALSE, col = "red")
  plot(-log(data$W_simulated_PC$pval), main='with_simulated_PC', ylab='-Log(p)', xlab='SNPs', ylim=c(0,20))
  abline(h = -log10(threshold), untf = FALSE, col = "red")
  plot(-log(data$W_computed_PC$pval), main='with_computed_PC', ylab='-Log(p)', xlab='SNPs', ylim=c(0,20))
  abline(h = -log10(threshold), untf = FALSE, col = "red")
  if(save == TRUE) dev.off()}

qq <- function(pvector, title="Quantile-quantile plot of p-values", spartan=F) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=0.25, main = title,xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)),ylim=c(0,max(e)))
  lines(e,e,col="red")
  plot}

summary_sim <- function(results,stratified_SNP_qtt,SNP_association) {
  before_associated = (length(SNP_association) - sum(SNP_association>0))  
  back = mapply(function(result, name) {
    fpstr = sum(result$signiff[0:stratified_SNP_qtt] != "ns")
    fpnonstr = sum(result$signiff[stratified_SNP_qtt:before_associated] != "ns")
    fn = sum(result$signiff[before_associated:length(SNP_association)] == "ns")
    matrix(c(fpstr,fpnonstr,fn),  dim=list(c("fpstr","fpnonstr","fn"), name))
  }, results, names(results) )
  back
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
compute_multiple(populations = 5, size = 1200, SNP_association = c(rep(0,100000), seq(1,2, by = 0.01)), stratification_rate = 0.05, times = 2)

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