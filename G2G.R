####Test scenario
test_G2G_setup <- function(study_design, scenario, fst_pop_strat=NA, fst_pop_bias=NA, fst_strain_strat=NA, fst_strain_bias=NA, tag = NULL, get_viral = get_viral_output) {
  tag = if(!is.null(tag)) paste0("_", tag)
  scenario_name = deparse(substitute(scenario))
  tag = paste(scenario_name, tag_me(fst_pop_strat), tag_me(fst_pop_bias), tag_me(fst_strain_strat), tag_me(fst_strain_bias), "beta",scenario$beta, nrow(study_design), sep = "_")
  print(paste0(Sys.time(), " : Scenario ",tag))
  
  SNP = get_SNP(study_design, scenario$rep, unlist(scenario$S_Stratified), unlist(scenario$S_Partial_Stratification), fst_pop_strat, unlist(scenario$S_Biased), unlist(scenario$S_Partial_Bias), fst_pop_bias)
  AA = get_aa(study_design, scenario$rep, unlist(scenario$Y_Stratified), unlist(scenario$Y_Partial_Stratification), fst_strain_strat, unlist(scenario$Y_Biased), unlist(scenario$Y_Partial_Bias), fst_strain_bias,
              unlist(scenario$Associated_Strains), unlist(scenario$Associated_Populations), beta=scenario$beta, SNP)
  
  threshold <<- 0.05/(ncol(SNP$SNP.data)+ncol(AA$AA.data))
  pvalues = analyse_G2G_setup(SNP$SNP.data,AA$AA.data, study_design)
  #summary_sim = summary_sim(pvalues, AA$AA.scenario)
  res = list(`study_design` = study_design, `params` = cbind(SNP$SNP.params, AA$AA.params), `pvalues` = pvalues, `scenario` = scenario)
  write(res, tag, paste0(getwd(),"/gen-data/",scenario_name,"/" ))
  res}

get_G2G_setup <- function(rep, s_stratified = NA, s_partial_strat = NA, s_biased = NA, s_partial_bias = NA, y_stratified = NA, y_partial_strat = NA, y_biased = NA, y_partial_bias = NA, associated_strains = NA, associated_populations = NA, beta=NA) {
  associated = !is.na(associated_strains) | !is.na(associated_populations)
  tbl_df(data_frame(S_Stratified = list(s_stratified), S_Biased = list(s_biased), 
                    S_Partial_Stratification = list(s_partial_strat), S_Partial_Bias = list(s_partial_bias),
                    Y_Stratified = list(y_stratified), Y_Biased = list(y_biased),
                    Y_Partial_Stratification = list(y_partial_strat), Y_Partial_Bias = list(y_partial_bias),
                    Associated_Strains = list(associated_strains), Associated_Populations = list(associated_populations), beta, rep) )}

analyse_G2G_setup <- function(SNPs, Y, study_design) {
  if(trace == TRUE) print(paste(Sys.time(),"Computing GLM", sep=" : "))
  WO_correction = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num])))[,4][2]))
  W_human_groups = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+study_design[,"Population"])))[,4][2]))
  W_strains_groups = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+study_design[,"Strain"])))[,4][2]))
  W_both_groups = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+study_design[,"Population"]+study_design[,"Strain"])))[,4][2]))
  if(trace == TRUE) print(paste(Sys.time(),"Parsing values", sep=" : "))   
  parse_pvalues(data.frame(WO_correction, W_human_groups, W_strains_groups, W_both_groups, row.names = colnames(SNPs)), threshold)}

##############Population
get_study_design <- function(sample_size, nb_pop, nb_strain) {
  get_structure <- function(sample_size, nb_substructure) {
    if(length(nb_substructure)>1) sample(1:length(nb_substructure), size = sample_size, prob = nb_substructure, replace = T)
    else sample(1:nb_substructure, size = sample_size, replace = T)}
  pop_structure = as.factor(paste0("P", get_structure(sample_size, nb_pop)))
  viral_pop_structure = as.factor(chartr("123456789", "ABCDEFGHI", get_structure(sample_size, nb_strain)))
  study_design = data.frame(`Population` = pop_structure, `Strain` = viral_pop_structure)
  study_design = arrange(study_design, Population, Strain)
  rownames(study_design) <- paste0(1:nrow(study_design),"_Pop_", study_design$Population,"_Strain_", study_design$Strain)
  study_design}

to_pop_structure <- function(study_design) {
  do.call(rbind, lapply(levels(study_design$Population), function(population) do.call(rbind, lapply(levels(study_design$Strain), function(strain) data_frame(strain,population, `nb` = sum(study_design[,"Strain"] == strain & study_design[,"Population"] == population))))))}

to_study_design <- function(structure) {
  do.call(rbind, lapply(1:nrow(structure), function(pop_num) {
    do.call(rbind, lapply(1:ncol(structure), function(strain_num) {
      data.frame(Strain = rep(chartr("123456789", "ABCDEFGHI", strain_num), structure[strain_num,pop_num]), Population = paste0("P",pop_num) ) }))})) }

get_aa <- function(study_design, size, stratified = NA, partial_strat = NA, fst_strat=NA, biased = NA, partial_bias = NA, fst_bias=NA, associated_strains = NA, associated_populations =NA, beta=NA, associated_SNPs=NA) {
  populations = levels(study_design$Population)
  strains = levels(study_design$Strain)
  
  associated = !is.na(associated_strains) | !is.na(associated_populations)
  if(associated && is.na(associated_SNPs)) stop("You must specify asociated SNPs")
  if(associated && ncol(associated_SNPs$SNP.data) != size) stop("You must define the same amount of associated SNPs as there is specified amino acids")
  associated_SNPs_id = if(associated) colnames(associated_SNPs$SNP.data) else NA
  associated_populations = if(!is.na(associated_populations)) if(associated_populations == "full") populations else if(associated_populations == "half") sample(populations, 1) else associated_populations
  associated_strains = if(!is.na(associated_strains)) if(associated_strains == "full") strains else if(associated_strains == "half") sample(strains, 1) else associated_strains
  
  id = get_id("aa",size)
  tag = c(paste0(round_me("Fstrat", fst_strat), part_me(partial_strat)), paste0(round_me("Fbias", fst_bias), part_me(partial_bias)), round_me("Beta",beta))
  tag = paste0(Filter(function(x) x!="", tag), collapse = "_")
  tag = if(tag=="") "homogenous" else tag
  
  AA.scenario = data_frame(Stratified = list(stratified), Biased = list(biased), Partial_Stratification = list(partial_strat), Partial_Bias = list(partial_bias),
                           Associated_Strains= list(associated_strains), Associated_Populations = list(associated_populations), associated, beta,
                           fst_strat, fst_bias, size, id = list(id), tag)
  
  AA.freq = get_frequencies(AA.scenario, strains,populations)
  rownames(AA.freq) <- id
  aa = get_viral_output(AA.freq, AA.scenario, study_design, associated_SNPs$SNP.data)
  colnames(aa) <- id
  
  list(`AA.data` = aa,`AA.associations` = data_frame(AA = as.factor(id), SNP= as.factor(associated_SNPs_id)), `AA.scenario` = AA.scenario)}

get_SNP <- function(study_design,size, stratified = NA, partial_strat = NA, fst_strat=NA, biased = NA, partial_bias = NA, fst_bias=NA) {
  populations = levels(study_design$Population)
  strains = levels(study_design$Strain)
  
  id =  get_id("SNP",size)
  tag = c(paste0(round_me("Fstrat", fst_strat), part_me(partial_strat)), paste0(round_me("Fbias", fst_bias), part_me(partial_bias)))
  tag = paste0(Filter(function(x) x!="", tag), collapse = "-")
  tag = if(tag=="") "homogenous" else tag
  
  SNP.scenario = data_frame(Stratified = list(stratified), Biased = list(biased), Partial_Stratification = list(partial_strat), Partial_Bias = list(partial_bias), fst_strat, fst_bias, size, id = list(id), tag)
  SNP.freq = get_frequencies(SNP.scenario, populations, strains)
  rownames(SNP.freq) <- id
  SNP.data = test_scenario(SNP.freq, study_design)
  colnames(SNP.data) <- id
  
  list(`SNP.data` = SNP.data, `SNP.scenario` = SNP.scenario)} #, `SNP.freq` = SNP.freq

get_id <- function(prefix,size) {
  id = if(Sys.getenv(prefix) == "") 0 else as.numeric(Sys.getenv(prefix))
  args = list(id + size)
  names(args) = prefix
  do.call(Sys.setenv, args)
  paste0(prefix, "_",(id+1):(id+size))}

round_me <- function(prefix, param) ifelse(is.na(param),"", paste0(prefix,"-",round(param, digits=2)))
part_me <- function(param) ifelse(is.na(param),"", paste0("_P-", paste0(param, collapse = ',')))
tag_me <- function(param) ifelse(is.na(param),"", paste(deparse(substitute(param)), round(param, digits=2), sep = "-"))

get_frequencies <- function(scenario, main_subgroup, secondary_subgroup) {
  nb_main_subgroup = length(main_subgroup)
  nb_secondary_subgroup = length(secondary_subgroup)
  
  strat = if(is.na(scenario$Stratified) | scenario$Stratified == "no") rep(0, nb_main_subgroup) 
  else if(scenario$Stratified == "full" | scenario$Stratified == "yes") sample(1:nb_main_subgroup, nb_main_subgroup)
  else rank(unlist(scenario$Stratified)) ###Here we will use user defined order
  ###Here we make a filter that will later remove unstratified group, this act on seconday_subgroup
  p_strat = if(!is.na(scenario$Partial_Stratification)) rep(strat, each = nb_secondary_subgroup) * rep(as.numeric(secondary_subgroup == scenario$Partial_Stratification), each = nb_secondary_subgroup)
  
  bias = if(scenario$Biased == "full" | scenario$Biased == "yes") sample(1:nb_secondary_subgroup, nb_secondary_subgroup)
  else rank(unlist(scenario$Biased)) ###Here we will use user defined order
  ###Here we make a filter that will later remove unbiased group, this act on seconday_subgroup
  p_bias = if(!is.na(scenario$Partial_Bias)) rep(bias, nb_main_subgroup) * rep(as.numeric(main_subgroup == scenario$Partial_Bias), each = nb_main_subgroup)
  
  freq = t(replicate(scenario$size, {
    RF = runif(1, 0.1, 0.4)
    F_strat = if(length(unique(strat)) >1) sort(get_AF(RF, scenario$fst_strat, length(unique(strat))), decreasing = TRUE)[strat] else rep(RF, nb_main_subgroup)
    F_bias =  if(length(unique(bias)) >1) c(sapply(F_strat, function(RF) sort(get_AF(RF, scenario$fst_bias, length(unique(bias))), decreasing = TRUE)[bias])) else rep(F_strat, each = nb_secondary_subgroup)
    F_bias[!is.null(p_strat) & p_strat == 0] = rep(F_strat, each = nb_main_subgroup)[p_strat == 0]
    F_bias[!is.null(p_bias) & p_bias == 0] = rep(F_strat, each = nb_secondary_subgroup)[p_bias == 0]
    F_bias}))
  colnames(freq) <- paste0("F_", rep(main_subgroup, each = nb_secondary_subgroup),".", secondary_subgroup)
  freq}

#Get alternate frequency according to Nicholas model
get_AF <- function(allele, fst, nb = 1) {
  if(is.na(fst)) stop("Specify fst")
  s1 = allele*(1-fst)/fst
  s2 = (1-allele)*(1-fst)/fst
  rbeta(n = nb, shape1 = s1,shape2 = s2)}

generate_SNPs_for_G2G <- function(SNP.freq, study_design) {
  nb_strains = length(levels(study_design$Strain))
  nb_populations = length(levels(study_design$Population))
  if(trace == TRUE) print(paste(Sys.time(),"Generating SNPs dsitribution with viral properties", sep=" : "))
  nb = summarize(group_by(study_design,Population, Strain), nb = n())$nb
  SNP = sapply(1:nrow(SNP.freq), function(snp_num) {
    p = as.data.frame(lapply(SNP.freq[snp_num,], function(alternate_allele) {c((1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2)}))    
    unlist(lapply(1:nb_populations, function(population_num) {
      iterator = ((population_num - 1 ) * nb_strains) + 1
      unlist(lapply(iterator:(iterator + nb_strains - 1), function(selector) {
        sample(c(0,1,2), size = nb[selector], prob = p[,selector], replace = TRUE) })) })) })
  matrix(data = SNP ,nrow = nrow(study_design), ncol = nrow(SNP.freq), dimnames = list(rownames(study_design), rownames(SNP.freq)))}

get_viral_output <- function(AA.freq, scenario, study_design, associated_SNPs) {
  strains = levels(study_design$Strain)
  populations = levels(study_design$Population)
  if(trace == TRUE) print(paste(Sys.time(),"Generating Viral output", sep=" : "))
  res = sapply(1:nrow(AA.freq), function(snp_num) {
    unlist(lapply(populations, function(population) {
      unlist(lapply(strains, function(strain) {
        filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
        AF = AA.freq[snp_num,paste0("F_",strain ,".",population)]
        Y = sample(0:1, size = sum(filter), prob = c(1 - AF,AF), replace = TRUE)
        if(scenario$associated && population %in% unlist(scenario$Associated_Populations) && strain %in% unlist(scenario$Associated_Strains)) {
          z = as.matrix(associated_SNPs[,snp_num])[which(filter)]*scenario$beta
          pr = 1/(1+exp(-z))
          Y = Y + unlist(lapply(pr -0.5, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))
          Y[Y>1] = 1}
        Y })) })) })
  rownames(res) <- paste0(rownames(study_design),"_",study_design[,"Population"],"_",study_design[,"Strain"]) 
  res}

parse_G2G_config <- function(...) {
  names = unlist(lapply(substitute(placeholderFunction(...))[-1], deparse))
  SNP_aa = list(...)
  SNP_aa = mapply(function(element, input_name) {
    element$input_name <- paste0(input_name,"_", 1:ncol(element[[1]]))
    element
  } ,  SNP_aa, names,SIMPLIFY = F) 
  AA = do.call(cbind, lapply(SNP_aa, function(ele) ele$AA.data))
  AA.scenarios = do.call(rbind, lapply(SNP_aa, function(ele) ele$AA.scenario))
  SNP = do.call(cbind, lapply(SNP_aa, function(ele) ele$SNP.data))
  SNP.scenarios = do.call(rbind, lapply(SNP_aa, function(ele) ele$SNP.scenario))
  associations = do.call(rbind, lapply(SNP_aa, function(ele) ele$AA.associations))
  list(`AA` = AA, `AA.scenarios` = AA.scenarios, `SNP` = SNP, `SNP.scenarios` = SNP.scenarios, `associations` = associations)}

analyse_G2G <- function(SNP, AA, study_design, nb_cpu, WO_correction = F, W_human_groups = F, W_strain_group = F, W_both_groups = F, W_human_PC = F, W_strain_PC = F, W_both_PC = F) { #SNP
  #analyse_G2G <- function(study_design) { #SNP
  if(trace == T) print(paste0(Sys.time()," : Computing PC"))
  SNP_PC = prcomp(SNP, .scale = F)
  AA_PC = prcomp(AA, .scale = F)
  SNP_PC = SNP_PC$x[,1:ifelse(ncol(SNP_PC$x)<5, ncol(SNP_PC$x), 5)]
  AA_PC = AA_PC$x[,1:ifelse(ncol(AA_PC$x)<5, ncol(AA_PC$x), 5)]
  #save(SNP_PC, AA_PC, file="savestate.RData")
  #load("savestate.RData")
  gc()
  if(trace == T) print(paste(Sys.time(),": Computing GLM on ",  nb_cpu, " CPU(s)"))
  cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
  ptm <- proc.time()
  
  analyse_AA <- function(Y, WO_correction, W_human_groups, W_strain_group, W_both_groups, W_human_PC, W_strain_PC, W_both_PC) {
    Filter(length, list(
      `WO_correction` = if(WO_correction == T) apply(SNP, 2, function(X) coef(summary(glm(Y~X)))[,4][2]),
      `W_human_groups` = if(W_human_groups == T) unlist(apply(SNP, 2, function(X) coef(summary(glm(Y~X+study_design[,"Population"])))[,4][2])),
      `W_strain_group` =  if(W_strain_group == T) apply(SNP, 2, function(X) coef(summary(glm(Y~X+study_design[,"Strain"])))[,4][2]),
      `W_both_groups` = if(W_both_groups == T) unlist(apply(SNP, 2, function(X) coef(summary(glm(Y~X+study_design[,"Population"]+study_design[,"Strain"])))[,4][2])),
      `W_human_PC` = if(W_human_PC == T) unlist(apply(SNP, 2, function(X) coef(summary(glm(Y~X+SNP_PC)))[,4][2])),
      `W_strain_PC` = if(W_strain_PC == T) unlist(apply(SNP, 2, function(X) coef(summary(glm(Y~X+AA_PC)))[,4][2])),
      `W_both_PC` = if(W_both_PC == T) unlist(apply(SNP, 2, function(X) coef(summary(glm(Y~X+SNP_PC+AA_PC)))[,4][2]))))}
  
  res = parApply(cl,AA, 2, analyse_AA, WO_correction, W_human_groups, W_strain_group, W_both_groups, W_human_PC, W_strain_PC, W_both_PC)
  time = proc.time() - ptm
  stopCluster(cl)
  if(trace == T) print(paste(Sys.time(),": Analysis took ", time, "Now saving data"))
  SNPcol = rep(colnames(SNP), length(res[[1]]))
  CorrectionCol = rep(names(res[[1]]), each = length(res[[1]][[1]]))
  res = do.call(cbind, mapply(function(aa, aa_id) {
    do.call(rbind, unname(lapply(aa, function(SNP) {
      setNames(data.frame(unname(SNP)), aa_id)})))
  }, res, names(res), SIMPLIFY = F))
  cbind(`SNP` = SNPcol, `Correction` =  CorrectionCol, res)}

plot_G2G <- function(res, associations, AA.scenarios, SNP.scenarios) {
  #save(res, associations, AA.scenarios, SNP.scenarios, file = "FullScenartio.RData")
  res_tidy = res %>% gather(AA, pvalue, -SNP, -Correction,factor_key = T)
  res_tidy = do.call(rbind, mapply(function(tag,id) {
    data.frame(SNP.tag = as.factor(tag), filter(res_tidy, SNP %in% unlist(id)))}, SNP.scenarios$tag, SNP.scenarios$id, SIMPLIFY = F))
  res_tidy = do.call(rbind, mapply(function(tag,id) {
    data.frame(AA.tag = as.factor(tag), filter(res_tidy, AA %in% unlist(id)))}, AA.scenarios$tag, AA.scenarios$id, SIMPLIFY = F))
  res_tidy$pvalue = -log10(res_tidy$pvalue)
  threshold = 0.05/(ncol(res)*nrow(res))
  
  plot_with_correction <- function(...,save = F) {
    correction=c(...)
    p <- ggplot(filter(res_tidy, Correction %in% correction), aes(AA.tag, pvalue, colour=SNP.tag, fill=Correction))
    p + geom_boxplot(outlier.color = "grey") + 
      geom_hline(yintercept = -log10(threshold), colour="red") + 
      labs(title = paste("pvalue for full model, with corrections:", paste(correction, collapse = " ")), y = "-log10(pval)", x="covariate") + 
      scale_y_continuous(limits = c(0, 15)) +
      geom_point(aes(y = pval_ass))
    if(save == T) ggsave(filename = paste0(getwd(),paste0(correction, collapse = "_"), Sys.time(), ".png"))}
  
  plot_with_correction_no_ass <- function(save = F) {
    p <- ggplot(res_tidy, aes(AA.tag, pvalue, fill=Correction))
    p + geom_boxplot(outlier.color = "grey") + 
      geom_hline(yintercept = -log10(threshold), colour="red") + 
      labs(title = paste("pvalue for full model, with corrections:", paste("", collapse = " ")), y = "-log10(pval)", x="covariate") + 
      scale_y_continuous(limits = c(0, 15))
    if(save == T) ggsave(filename = paste0(getwd(),paste0(collapse = "_"), Sys.time(), ".png"))}
  
  ###ASSOCIATION
  if(F) {
    associated = inner_join(associations, res_tidy, by = c("SNP", "AA"))
    mutated = mutate(associated, pval_ass = pvalue)
    associated$pval_ass <- mutated$pval_ass
    res_tidy = right_join(associated, res_tidy, by = c("SNP", "AA", "AA.tag", "SNP.tag", "Correction", "pvalue"))
    #levels(associations$AA) <- (levels(as.factor(res_tidy$AA)))
    #levels(associations$SNP) <- (levels(as.factor(res_tidy$SNP)))
    associations$associated <- T
    res_tidy = right_join(associations, res_tidy, by = c("SNP", "AA"))
  }
  
  ###???
  if(F){
    associated_pvalue = do.call(rbind, mapply(function(AA.id, SNP.id) {
      res = filter(res_tidy, SNP == SNP.id & AA == AA.id)
      data_frame(res$AA, res$SNP, res$pvalue, res$Correction)
    }, associations$AA, associations$SNP))}
  
  plot_with_correction_no_ass(save = T)
  #plot_with_correction("WO_correction", "W_human_groups", "W_human_PC", save = T)
  #plot_with_correction("WO_correction", "W_strain_group", "W_strain_PC", save = T)
  #plot_with_correction("WO_correction", "W_both_groups", "W_both_PC", save = T)
  #lapply(levels(res_tidy$Correction), plot_with_correction, save = T)
}

plot_G2G_exp <- function(res, associations, AA.scenarios, SNP.scenarios) {
  res_tidy = res %>% gather(AA, pvalue, -SNP, -Correction,factor_key = T)
  res_tidy = do.call(rbind, mapply(function(tag,id) {
    data.frame(SNP.tag = as.factor(tag), filter(res_tidy, SNP %in% unlist(id)))}, SNP.scenarios$tag, SNP.scenarios$id, SIMPLIFY = F))
  res_tidy = do.call(rbind, mapply(function(tag,id) {
    data.frame(AA.tag = as.factor(tag), filter(res_tidy, AA %in% unlist(id)))}, AA.scenarios$tag, AA.scenarios$id, SIMPLIFY = F))
  res_tidy$pvalue = -log10(res_tidy$pvalue)
  threshold = 0.05/(ncol(res)*nrow(res))
  p <- ggplot(res_tidy, aes(AA.tag, pvalue, fill=Correction))
  p + geom_boxplot(outlier.color = "grey") + 
    geom_hline(yintercept = -log10(threshold), colour="red") + 
    labs(title = paste("pvalue for full model, with corrections:", paste("", collapse = " ")), y = "-log10(pval)", x="covariate") + 
    scale_y_continuous(limits = c(0, 15))}