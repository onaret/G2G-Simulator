###Data generation
generate_SNPs_for_GWAS <- function(populations, neutral=0, neutral_S_rate=0, causal_S=c(), causal_NS=c(), fst_strat) {
  neutral_S = neutral_S_rate * neutral
  SNP_params = get_SNP_frequencies_for_GWAS(neutral - neutral_S, neutral_S, causal_S, causal_NS, populations, fst_strat)
  if(trace == TRUE) print(paste(Sys.time(),"Generating SNPs dsitribution", sep=" : "))
  SNPs = apply(SNP_params, 1, function(SNP_param){
    p = get_genotype_probability(SNP_param[1:nrow(populations)], SNP_param["R"])
    unlist(lapply(1:nrow(populations), function(population) {
      cases=sample(c(0,1,2), size = populations[population,"Case"], prob = p[1:3, population], replace = TRUE)
      controls=sample(c(0,1,2), size = populations[population,"Control"], prob = p[4:6, population], replace = TRUE)
      c(cases, controls)}))})
  rownames(SNPs) <- unlist(lapply(1:nrow(populations), function(population) {
    c(if(populations[population,"Case"] > 0) {paste('Case', 'P', population, 'sample', 1:populations[population,"Case"],sep='_')},
      if(populations[population,"Control"] > 0) {paste('Control', 'P', population, 'sample', 1:populations[population,"Control"],sep='_')})}))
  list(`data` = SNPs, `params` = SNP_params)}

get_SNP_frequencies_for_GWAS <- function(neutral_NS, neutral_S, causal_S, causal_NS, populations, fst_strat) {
  if(trace == TRUE) print(paste(Sys.time(),"Generating allele frequency", sep=" : "))
  # nb_pop = length(levels(populations$Population))
  nb_pop = ncol(populations)
  f_neutral_NS = if(neutral_NS >0) data.frame(t(data.frame(replicate(neutral_NS, rep(get_AF(allele = runif(1,0.4,0.5), fst = fst_strat), nb_pop),simplify = F))), Stratified = F, Causal = F, R= NA)
  f_neutral_S = if(neutral_S>0) data.frame(t(data.frame(replicate(neutral_S, get_AF(allele = runif(1,0.4,0.5), fst = fst_strat, nb = nb_pop),simplify = F))), Stratified = T, Causal = F, R= NA)
  f_causal_NS = if(length(causal_NS)>0) data.frame(t(data.frame(replicate(length(causal_NS), rep(get_AF(allele = runif(1,0.4,0.5), fst = fst_strat), nb_pop),simplify = F))), Stratified = F, Causal = T, R= causal_NS)
  f_causal_S = if(length(causal_S)>0) data.frame(t(data.frame(replicate(length(causal_S), get_AF(allele = runif(1,0.4,0.5), fst = fst_strat, nb = nb_pop),simplify = F))), Stratified = T, Causal = T, R= causal_S)
  SNP_params = rbind(f_neutral_NS, f_neutral_S, f_causal_NS, f_causal_S)
  colnames(SNP_params) <- c(paste0("P", 1:nb_pop), "Stratified", "Causal", "R")
  SNP_names = c(if(neutral_NS >0) paste0('Snp_NS_',1:neutral_NS,"_R0"),
                if(neutral_S>0) paste0('Snp_S_',1:neutral_S,"_R0"),
                if(length(causal_NS)>0) paste0('Snp_NS_',1:length(causal_NS),"_R_", causal_NS),
                if(length(causal_S)>0) paste0('Snp_S_',1:length(causal_S),"_R_", causal_S))
  rownames(SNP_params) <- SNP_names
  SNP_params}

get_genotype_probability <- function(AFs, R = NA) {
  as.data.frame(lapply(AFs, function(alternate_allele) {
    c(`Case` =
        if(!is.na(R)) {
          p_associated = c((1-alternate_allele)^2, 2*as.numeric(R)*alternate_allele*(1-alternate_allele), as.numeric(R)^2 * alternate_allele^2 )
          p_associated/(sum(p_associated) + (1-alternate_allele^2))}
      else { 
        c( (1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2 )},
      `Control` = c( (1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2))}))}

###Data analyses
analyse_GWAS <- function(SNPs, populations, nb_pc) {
  simulated_PC = rep(1:nrow(populations), populations[1:nrow(populations),"Control"] + populations[1:nrow(populations), "Case"])
  simulated_PC = chartr("123456789", "ABCDEFGHI", simulated_PC)
  y = unlist(lapply(1:nrow(populations), function(nb_pop) c( rep(1, populations[nb_pop,"Case"]), rep(0, populations[nb_pop,"Control"]))))
  if(trace == TRUE) print(paste(Sys.time(),"Computing PCA", sep=" : "))
  SNP_PC = prcomp(SNPs, scale. = FALSE)
  SNP_PC = SNP_PC$x[,1:ifelse(ncol(SNP_PC$x)<nb_pc, ncol(SNP_PC$x), nb_pc)]
  if(trace == TRUE) print(paste(Sys.time(),"Computing logistic regression with computed PC", sep=" : "))
  W_computed_PC = apply(SNPs,2, function(SNP) coef(summary(glm(y~SNP+SNP_PC, family = binomial)))[,4][2])
  if(trace == TRUE) print(paste(Sys.time(),"Computing logistic regression with simulated  PC", sep=" : "))
  W_simulated_PC = apply(SNPs,2, function(SNP) coef(summary(glm(y~SNP+simulated_PC, family = binomial)))[,4][2])
  if(trace == TRUE) print(paste(Sys.time(),"Computing logistic regression without PC", sep=" : "))
  WO_correction = apply(SNPs,2, function(SNP) coef(summary(glm(y~SNP, family = binomial)))[,4][2])
  data.frame(WO_correction, W_simulated_PC, W_computed_PC, row.names = colnames(SNPs))}