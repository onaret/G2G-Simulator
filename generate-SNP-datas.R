####Function
multiple_SNP_scenarios <- function(populations, neutral, neutral_S_rat, causal_S=c(), causal_NS=c(), repetitions) {
  #This condition is ugly
  cl <- makeCluster(CPU, outfile = "output.txt", type = "FORK")
  if(is.list(populations) && is.list(populations[[1]])) {
    parLapply(cl, 1:(repetitions*length(populations)), function(repetition, populations, names){
      pop_index = repetition%%length(populations) + 1
      print(paste(Sys.time(),": For population", names[pop_index], "executing repetition", trunc(repetition/length(populations))+1))
      write(SNP_scenario(populations[pop_index], neutral = neutral, neutral_S_rate, causal_S, causal_NS), repetition, names[pop_index])
    },populations, names(populations) )}
  else {
    lapply(1:repetitions, function(repetition){
      if(trace == TRUE) print(paste(Sys.time(),": Executing repetition", repetition))
      write(SNP_scenario(populations, neutral, neutral_S_rate, causal_S, causal_NS), repetition)})}
  if(trace == TRUE) print(paste(Sys.time(),"Done", sep=" : "))
  stopCluster(cl)}

SNP_scenario <- function(populations, neutral=0, neutral_S_rate=0, causal_S=c(), causal_NS=c(), fst_strat) {
  populations = generate_population_structure_for_CC(populations)
  neutral_S = get_stratified_SNPs_qtt(neutral_S_rate,neutral)
  SNP_params = generate_SNPs_frequencies(neutral - neutral_S, neutral_S, causal_S, causal_NS, populations, fst_strat)
  SNPs = generate_SNPs(SNP_params, populations)
  pvalues = analyse_SNP(SNPs, populations)
  summary_sim = summary_sim(pvalues, SNP_params)
  trace_plot(pvalues)
  
  populations = rbind(populations, matrix(c(sum(populations[,"Case"]), sum(populations[,"Control"])),ncol = 2, nrow = 1, dimnames = list("Total") ) )
  list(#`SNPs` = SNPs,
    `study_design` = populations,
    `SNP_params` = SNP_params,
    `pvalues` = pvalues,
    `summary_sim` = summary_sim)}

get_stratified_SNPs_qtt <- function(stratification_rate, neutral_SNPs) {
  if(stratification_rate>1) stratification_rate = stratification_rate/100
  if(stratification_rate>100) stratification_rate = 1
  stratification_rate * neutral_SNPs}

generate_population_structure_for_CC <- function(populations) {
  populations = {
    if(!is.list(populations)) {
      possible_nb_pop = populations["min"]:populations["max"]
      size = populations["size"]
      populations = sample(c(possible_nb_pop), size=1, prob=rep(1/length(possible_nb_pop),length(possible_nb_pop)), replace = TRUE)
      populations = min(size, populations)
      populations = runif(populations,0,1)
      populations = round(populations/sum(populations) * size)
      populations[1] = populations[1] + size - sum(populations)
      t(sapply(populations, function(strat_size, pop_id) {
        case_nb = round(strat_size*runif(1,0,1))
        #case_nb = round(strat_size*abs(rnorm(5,0.5,0.2)))
        c(case_nb,(strat_size-case_nb))}))}
    else if(is.null(populations) && is.null(size)) stop("You must set size or population arguments...")
    else t(as.data.frame(populations))}
  colnames(populations) <- c("Case", "Control")
  rownames(populations) <- paste0("P", 1:nrow(populations))
  populations
  #do.call(rbind, lapply(1:nrow(populations), function(pop_num) {rbind(data.frame(CC = rep("Case", populations[pop_num,1]), Population = paste0("P",pop_num) ),data.frame(CC = rep("Control", populations[pop_num,2]), Population = paste0("P",pop_num) )) }))
  }

generate_SNPs_frequencies <- function(neutral_NS, neutral_S, causal_S, causal_NS, populations, fst_strat) {
  if(trace == TRUE) print(paste(Sys.time(),"Generating allele frequency", sep=" : "))
 # nb_pop = length(levels(populations$Population))
  nb_pop = ncol(populations)
  f_neutral_NS = if(neutral_NS >0) data.frame(Stratified = F, Causal = F, R= NA , t(data.frame(replicate(neutral_NS, rep(get_AF(allele = runif(1,0.4,0.5), fst = fst_strat), nb_pop),simplify = F))))
  f_neutral_S = if(neutral_S>0) data.frame(Stratified = T, Causal = F, R= NA , t(data.frame(replicate(neutral_S, get_AF(allele = runif(1,0.01,0.99), fst = fst_strat, nb = nb_pop),simplify = F))))
  f_causal_NS = if(length(causal_NS)>0) data.frame(Stratified = F, Causal = T, R= causal_NS , t(data.frame(replicate(length(causal_NS), rep(get_AF(allele = runif(1,0.4,0.5), fst = fst_strat), nb_pop),simplify = F))))
  f_causal_S = if(length(causal_S)>0) data.frame(Stratified = T, Causal = T, R= causal_S,  t(data.frame(replicate(length(causal_S), get_AF(allele = runif(1,0.01,0.99), fst = fst_strat, nb = nb_pop),simplify = F))))
  SNP_params = rbind(f_neutral_NS, f_neutral_S, f_causal_NS, f_causal_S)
  colnames(SNP_params) <- c("Stratified", "Causal", "R", paste0("P", 1:nb_pop))
  SNP_names = c(if(neutral_NS >0) paste0('Snp_NS_',1:neutral_NS,"_R0"),
    if(neutral_S>0) paste0('Snp_S_',1:neutral_S,"_R0"),
    if(length(causal_NS)>0) paste0('Snp_NS_',1:length(causal_NS),"_R_", causal_NS),
    if(length(causal_S)>0) paste0('Snp_S_',1:length(causal_S),"_R_", causal_S))
  rownames(SNP_params) <- SNP_names
  SNP_params}

generate_SNPs <- function(SNP_params, populations) {
  if(trace == TRUE) print(paste(Sys.time(),"Generating SNPs dsitribution", sep=" : "))
  data = data = apply(SNP_params, 1, function(SNP_param){
    p = get_SNP_probability(SNP_param[1:length(populations)], SNP_param["R"], CC=TRUE)
    unlist(lapply(1:nrow(populations), function(population) {
      cases=sample(c(0,1,2), size = populations[population,"Case"], prob = p[1:3, population], replace = TRUE)
      controls=sample(c(0,1,2), size = populations[population,"Control"], prob = p[4:6, population], replace = TRUE)
      c(cases, controls)}))})
  rownames(data) <- get_samples_name(populations)
  data}

get_samples_name <- function(populations) {
  unlist(lapply(1:nrow(populations), function(population) {
    c(if(populations[population,"Case"] > 0) {paste('Case', 'P', population, 'sample', 1:populations[population,"Case"],sep='_')},
      if(populations[population,"Control"] > 0) {paste('Control', 'P', population, 'sample', 1:populations[population,"Control"],sep='_')})}))}

get_SNP_probability <- function(AFs, R = NA, CC = FALSE) {
  as.data.frame(lapply(AFs, function(alternate_allele) {
    c(`Case` =
      if(!is.na(R)) {
        p_associated = c((1-alternate_allele)^2, 2*as.numeric(R)*alternate_allele*(1-alternate_allele), as.numeric(R)^2 * alternate_allele^2 )
        p_associated/(sum(p_associated) + (1-alternate_allele^2))}
      else { 
        c( (1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2 )},
      `Control` = c( (1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2))}))}

analyse_SNP <- function(SNPs, populations) {
  simulated_PC = rep(1:nrow(populations), populations[1:nrow(populations),"Control"] + populations[1:nrow(populations), "Case"])
  simulated_PC = chartr("123456789", "ABCDEFGHI", simulated_PC)
  y = unlist(lapply(1:nrow(populations), function(nb_pop) c( rep(1, populations[nb_pop,"Case"]), rep(0, populations[nb_pop,"Control"]))))
  if(trace == TRUE) print(paste(Sys.time(),"Computing PCA", sep=" : "))
  ca = prcomp(SNPs, scale. = FALSE) 
  if(trace == TRUE) print(paste(Sys.time(),"Computing logistic regression with computed PC", sep=" : "))
  W_computed_PC = apply(SNPs,2, function(SNP) coef(summary(glm(y~SNP+ca$x[,1:5])))[,4][2])
  if(trace == TRUE) print(paste(Sys.time(),"Computing logistic regression with simulated  PC", sep=" : "))
  W_simulated_PC = apply(SNPs,2, function(SNP) coef(summary(glm(y~SNP+simulated_PC)))[,4][2])
  if(trace == TRUE) print(paste(Sys.time(),"Computing logistic regression without PC", sep=" : "))
  WO_PC = apply(SNPs,2, function(SNP) coef(summary(glm(y~SNP)))[,4][2])
  parse_pvalues(data.frame(WO_PC, W_simulated_PC, W_computed_PC, row.names = colnames(SNPs)), threshold)}



####Constants
C0 = list(`P1` = c(`case` = 1250, `control` = 1250), `P2` = c(`case` = 1250, `control`  = 1250))
C1 = list(`P1` = c(`case` = 200, `control` = 400), `P2` = c(`case` = 400, `control`  = 200))
C2 = list(`P1` = c(`case` = 400, `control` = 200), `P2` = c(`case` = 200, `control` = 400))
C3 = list(`P1` = c(`case` = 300, `control` = 0), `P2` = c(`case` = 300, `control` = 600))
C4 = list(`P1` = c(`case` = 300, `control` = 200), `P2` = c(`case` = 200, `control` = 100), `P3` = c(`case` = 100, `control` = 300))
C5 = list(`P1` = c(`case` = 200, `control` = 0), `P2` = c(`case` = 400, `control` = 200), `P3` = c(`case` = 0, `control` = 400))
AllPop = list(`C1` = C1, `C2`= C2, `C3` = C3, `C4` = C4, `C5` = C5)
threshold = 5e-07

fcoeff  = 0.01 ##### Wright's coefficient for inbreeding
trace = TRUE
populations =  generate_population_structure_for_CC(C0)
SNP_scenario(C0,  neutral = 100000, 
             neutral_S_rate = 0.05, 
             causal_NS = seq(1,2, by = 0.05),
             causal_S = seq(1,2, by = 0.05), fcoeff)
