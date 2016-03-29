####Test scenario
test_scenario <- function(study_design, scenario, fst_pop_strat=NA, fst_pop_bias=NA, fst_strain_strat=NA, fst_strain_bias=NA, tag = NULL, get_viral = get_viral_output) {
  tag = if(!is.null(tag)) paste0("_", tag)
  scenario_name = deparse(substitute(scenario))
  tag = paste(scenario_name, tag_me(fst_pop_strat), tag_me(fst_pop_bias), tag_me(fst_strain_strat), tag_me(fst_strain_bias), "beta",scenario$beta, nrow(study_design), sep = "_")
  print(paste0(Sys.time(), " : Scenario ",tag))
  
  SNP = get_SNP(study_design, scenario$rep, unlist(scenario$S_Stratified), unlist(scenario$S_Partial_Stratification), fst_pop_strat, unlist(scenario$S_Biased), unlist(scenario$S_Partial_Bias), fst_pop_bias)
  aa = get_aa(study_design, scenario$rep, unlist(scenario$Y_Stratified), unlist(scenario$Y_Partial_Stratification), fst_strain_strat, unlist(scenario$Y_Biased), unlist(scenario$Y_Partial_Bias), fst_strain_bias,
              unlist(scenario$Associated_Strains), unlist(scenario$Associated_Populations), beta=scenario$beta, SNP)
  pvalues = analyse_scenario(SNP$SNP.data,aa$aa.data, study_design)
  summary_sim = summary_sim(pvalues, aa$aa.params)
  res = list(`study_design` = study_design, `params` = cbind(SNP$SNP.params, aa$aa.params), `pvalues` = pvalues, `summary_sim` = summary_sim, `scenario` = scenario)
  write(res, tag, paste0(getwd(),"/gen-data/",scenario_name,"/" ))
  res}

get_scenario <- function(rep, s_stratified = NA, s_partial_strat = NA, s_biased = NA, s_partial_bias = NA, y_stratified = NA, y_partial_strat = NA, y_biased = NA, y_partial_bias = NA, associated_strains = NA, associated_populations = NA, beta=NA) {
  associated = !is.na(associated_strains) | !is.na(associated_populations)
  tbl_df(data_frame(S_Stratified = list(s_stratified), S_Biased = list(s_biased), 
                    S_Partial_Stratification = list(s_partial_strat), S_Partial_Bias = list(s_partial_bias),
                    Y_Stratified = list(y_stratified), Y_Biased = list(y_biased),
                    Y_Partial_Stratification = list(y_partial_strat), Y_Partial_Bias = list(y_partial_bias),
                    Associated_Strains = list(associated_strains), Associated_Populations = list(associated_populations), beta, rep) )}

analyse_scenario <- function(SNPs, Y, study_design) {
  if(trace == TRUE) print(paste(Sys.time(),"Computing GLM", sep=" : "))
  WO_correction = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num])))[,4][2]))
  W_human_groups = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+study_design[,"Population"])))[,4][2]))
  W_strains_groups = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+study_design[,"Strain"])))[,4][2]))
  W_both_groups = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+study_design[,"Population"]+study_design[,"Strain"])))[,4][2]))
  if(trace == TRUE) print(paste(Sys.time(),"Parsing values", sep=" : "))   
  parse_pvalues(data.frame(WO_correction, W_human_groups, W_strains_groups, W_both_groups, row.names = colnames(SNPs)), threshold)}

##############
get_study_design <- function(sample_size, nb_strain, nb_pop) {
  get_structure <- function(sample_size, nb_substructure) {
    if(length(nb_substructure)>1) sample(1:length(nb_substructure), size = sample_size, prob = nb_substructure, replace = T)
    else sample(1:nb_substructure, size = sample_size, replace = T)}
  pop_structure = as.factor(paste0("P", get_structure(sample_size, nb_pop)))
  viral_pop_structure = as.factor(chartr("123456789", "ABCDEFGHI", get_structure(sample_size, nb_strain)))
  study_design = data.frame(`Population` = pop_structure, `Strain` = viral_pop_structure)
  study_design = arrange(study_design, Population, Strain)
  rownames(study_design) <- paste0(1:nrow(study_design),"_Pop_", study_design$Population,"_Strain_", study_design$Strain)
  study_design}

get_aa <- function(study_design, size, stratified = NA, partial_strat = NA, fst_strat=NA, biased = NA, partial_bias = NA, fst_bias=NA, associated_strains = NA, associated_populations =NA, beta=NA, associated_SNPs=NA) {
  id = get_id("aa",size)
  populations = levels(study_design$Population)
  strains = levels(study_design$Strain)
  associated = !is.na(associated_strains) | !is.na(associated_populations)
  associated_populations = if(!is.na(associated_populations)) if(associated_populations == "full") populations else if(associated_populations == "half") sample(populations, 1) else associated_populations
  associated_strains = if(!is.na(associated_strains)) if(associated_strains == "full") strains else if(associated_strains == "half") sample(strains, 1) else associated_strains
  scenarios = data_frame(id, Stratified = list(stratified), Biased = list(biased), Partial_Stratification = list(partial_strat), Partial_Bias = list(partial_bias),
                         Associated_Strains= list(associated_strains), Associated_Populations = list(associated_populations), associated, beta,
                         fst_strat, fst_bias, size)
  freq_aa = get_frequencies(scenarios, strains,populations, fst_strat, fst_bias)
  aa = get_viral_output(freq_aa, study_design, associated_SNPs$SNP.data)
  list(`aa.data` = aa,`aa.freq` = freq_aa, `aa.params` = scenarios)}

get_SNP <- function(study_design,size, stratified = NA, partial_strat = NA, fst_strat=NA, biased = NA, partial_bias = NA, fst_bias=NA) {
  id = get_id("SNP",size)
  scenarios = data_frame(Stratified = list(stratified), Biased = list(biased), Partial_Stratification = list(partial_strat), Partial_Bias = list(partial_bias), fst_strat, fst_bias, size)
  freq_SNP = get_frequencies(scenarios, levels(study_design$Population), levels(study_design$Strain), fst_strat, fst_bias)
  list(`SNP.data` = generate_SNPs(freq_SNP, study_design), `SNP.freq` = freq_SNP, `SNP.params` = scenarios)}

get_id <- function(prefix,size) {
  id = if(Sys.getenv(prefix) == "") 1 else as.numeric(Sys.getenv(prefix))
  args = list(id + size)
  names(args) = prefix
  do.call(Sys.setenv, args)
  id:(id+size)
  }

tag_me <- function(param) ifelse(is.na(param),"", paste(deparse(substitute(param)), round(param, digits=2), sep = "-"))

get_frequencies <- function(scenarios, main_subgroup, secondary_subgroup, fst_strat, fst_bias) {
  nb_main_subgroup = length(main_subgroup)
  nb_secondary_subgroup = length(secondary_subgroup)
  get_AF <- function(allele, fst, nb = 1) {
    if(is.na(fst)) stop("Specify fst")
    s1 = allele*(1-fst)/fst
    s2 = (1-allele)*(1-fst)/fst
    rbeta(n = nb, shape1 = s1,shape2 = s2)}
  
  strat = if(is.na(scenarios$Stratified) | scenarios$Stratified == "no") rep(0, nb_main_subgroup) 
  else if(scenarios$Stratified == "full" | scenarios$Stratified == "yes") sample(1:nb_main_subgroup, nb_main_subgroup)
  else rank(unlist(scenarios$Stratified)) ###Here we will use user defined order
  ###Here we make a filter that will later remove unstratified group, this act on seconday_subgroup
  p_strat = if(!is.na(scenarios$Partial_Stratification)) rep(strat, each = nb_secondary_subgroup) * rep(as.numeric(secondary_subgroup == scenarios$Partial_Stratification), each = nb_secondary_subgroup)
  
  bias = if(scenarios$Biased == "full" | scenarios$Biased == "yes") sample(1:nb_secondary_subgroup, nb_secondary_subgroup)
  else rank(unlist(scenarios$Biased)) ###Here we will use user defined order
  ###Here we make a filter that will later remove unbiased group, this act on seconday_subgroup
  p_bias = if(!is.na(scenarios$Partial_Bias)) rep(bias, nb_main_subgroup) * rep(as.numeric(main_subgroup == scenarios$Partial_Bias), each = nb_main_subgroup)
  
  freq = t(replicate(scenarios$size, {
    RF = runif(1, 0.1, 0.4)
    F_strat = if(length(unique(strat)) >1) sort(get_AF(RF, fst_strat, length(unique(strat))), decreasing = TRUE)[strat] else rep(RF, nb_main_subgroup)
    F_bias =  if(length(unique(bias)) >1) c(sapply(F_strat, function(RF) sort(get_AF(RF, fst_bias, length(unique(bias))), decreasing = TRUE)[bias])) else rep(F_strat, each = nb_secondary_subgroup)
    F_bias[!is.null(p_strat) & p_strat == 0] = rep(F_strat, each = nb_main_subgroup)[p_strat == 0]
    F_bias[!is.null(p_bias) & p_bias == 0] = rep(F_strat, each = nb_secondary_subgroup)[p_bias == 0]
    F_bias}))
  colnames(freq) <- paste0("F_", rep(main_subgroup, each = nb_secondary_subgroup),".", secondary_subgroup)
  rownames(freq) <- paste0("Nb_", 1:scenarios$size) 
  cbind(freq, scenarios)}

generate_SNPs <- function(SNP_params, study_design) {
  strains = levels(study_design$Strain)
  populations = levels(study_design$Population)
  if(trace == TRUE) print(paste(Sys.time(),"Generating SNPs dsitribution with viral properties", sep=" : "))
  data = matrix(nrow = nrow(study_design), ncol = nrow(SNP_params), dimnames = list(rownames(study_design), rownames(SNP_params)))
  SNP_params = select(SNP_params, +starts_with("F_"))
  sapply(1:nrow(SNP_params), function(snp_num) {
    p = as.data.frame(lapply(SNP_params[snp_num,], function(alternate_allele) {c((1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2)}))    
    unlist(lapply(populations, function(population) {
      unlist(lapply(strains, function(strain) {
        prob = unname(unlist(select(p, contains(paste0("F_", population,"." ,strain) ))))
        sample(c(0,1,2), size = nrow(filter(study_design, Population == population, Strain == strain)), prob = prob, replace = TRUE) })) })) })}

get_viral_output <- function(viral_params, study_design, associated_SNPs) {
  strains = levels(study_design$Strain)
  populations = levels(study_design$Population)
  if(trace == TRUE) print(paste(Sys.time(),"Generating Viral output", sep=" : "))
  sapply(1:nrow(viral_params), function(snp_num) {
    unlist(lapply(populations, function(population) {
      unlist(lapply(strains, function(strain) {
        filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
        AF = viral_params[snp_num,paste0("F_",strain ,".",population)]
        Y = sample(0:1, size = sum(filter), prob = c(1 - AF,AF), replace = TRUE)
        if(viral_params[snp_num,"associated"] && population %in% unlist(viral_params[snp_num,"Associated_Populations"]) && strain %in% unlist(viral_params[snp_num,"Associated_Strains"])) {
          z = as.matrix(associated_SNPs[,snp_num])[which(filter)]*viral_params[snp_num,"beta"]
          pr = 1/(1+exp(-z))
          Y = Y + unlist(lapply(pr -0.5, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))
          Y[Y>1] = 1
        }
        Y })) })) })}

scenario_full <- function(study_design, ...) {
  names = unlist(lapply(substitute(placeholderFunction(...))[-1], deparse))
  SNP_aa = list(...)
  names_SNP_aa = mapply(function(element, input_name) {
    colnames(element[[1]]) <- paste0(input_name, 1:ncol(element[[1]]))
    rownames(element[[2]]) <- paste0(input_name, 1:nrow(element[[2]]))
    element
    } ,  SNP_aa, names) 
  aa = do.call(cbind, lapply(SNP_aa, function(ele) ele$aa.data))
  SNP = do.call(cbind, lapply(SNP_aa, function(ele) ele$SNP.data))
  
  aa.params = do.call(rbind, lapply(SNP_aa, function(ele) ele$aa.params))
  SNP.params = do.call(rbind, lapply(SNP_aa, function(ele) ele$SNP.params))
  
  if(trace == TRUE) print(paste0(Sys.time()," : Computing PC"))
  SNP_PC = prcomp(SNP, .scale = FALSE)
  aa_PC = prcomp(aa, .scale = FALSE)
  cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
  if(trace == TRUE) print(paste0(Sys.time()," : Starting analysis on ",cl, " CPU(s)"))
  
  res = parApply(cl, as.array(aa), 2, function(Y) as.data.frame(analyse_full(SNP, Y, study_design, SNP_PC, aa_PC)))
  #res = apply(as.array(aa), 2, function(Y) analyse_viral(SNP, Y, study_design, SNP_PC, aa_PC))
  stopCluster(cl)
  res}

analyse_full <- function(SNP, Y, study_design, SNP_PC, aa_PC) {
  nb_PCs_SNP = ifelse(ncol(SNP_PC$x)<5, ncol(SNP_PC$x), 5)
  nb_PCs_aa = ifelse(ncol(aa_PC$x)<5, ncol(aa_PC$x), 5)
  
  if(trace == TRUE) print(paste(Sys.time(),"Computing GLM", sep=" : "))
  WO_correction = unlist(lapply(1:ncol(SNP), function(SNP_Y_num) coef(summary(glm(Y~SNP[,SNP_Y_num])))[,4][2]))
  W_human_groups = unlist(lapply(1:ncol(SNP), function(SNP_Y_num) coef(summary(glm(Y~SNP[,SNP_Y_num]+study_design[,"Population"])))[,4][2]))
  W_strains_groups = unlist(lapply(1:ncol(SNP), function(SNP_Y_num) coef(summary(glm(Y~SNP[,SNP_Y_num]+study_design[,"Strain"])))[,4][2]))
  W_both_groups = unlist(lapply(1:ncol(SNP), function(SNP_Y_num) coef(summary(glm(Y~SNP[,SNP_Y_num]+study_design[,"Population"]+study_design[,"Strain"])))[,4][2]))
  W_human_PCs = unlist(lapply(1:ncol(SNP), function(SNP_Y_num) coef(summary(glm(Y~SNP[,SNP_Y_num]+SNP_PC$x[,1:nb_PCs_SNP])))[,4][2]))
  W_strain_PC = unlist(lapply(1:ncol(SNP), function(SNP_Y_num) coef(summary(glm(Y~SNP[,SNP_Y_num]+aa_PC$x[,1:nb_PCs_aa])))[,4][2]))
  W_both_PC = unlist(lapply(1:ncol(SNP), function(SNP_Y_num) coef(summary(glm(Y~SNP[,SNP_Y_num]+SNP_PC$x[,1:nb_PCs_SNP]+aa_PC$x[,1:nb_PCs_aa])))[,4][2]))
  if(trace == TRUE) print(paste(Sys.time(),"Parsing values", sep=" : "))   
  threshold = 5e-5
  parse_pvalues(data.frame(WO_correction, W_human_groups, W_strains_groups, W_both_groups, W_human_PCs, W_strain_PC, W_both_PC, row.names = colnames(SNP)), threshold)}