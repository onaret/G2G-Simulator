###G2G common methods
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

get_AA <- function(study_design, AA.scenarios, SNP.data, associated_SNPs = NULL) {
  populations = levels(study_design$Population)
  strains = levels(study_design$Strain)
  AA = lapply(1:nrow(AA.scenarios), function(scenario_num) {
    AA.scenario = AA.scenarios[scenario_num,, drop = F]
    associated_populations = unlist(AA.scenario$Associated_Populations)
    associated_strains = unlist(AA.scenario$Associated_Strains)
    AA.scenario$Associated_Populations = list(
      if(!is.na(associated_populations))
        if(associated_populations == "full") populations
        else if(associated_populations == "half") sample(populations, 1)
        else associated_populations)
    AA.scenario$Associated_Strains = list(
      if(!is.na(associated_strains))
        if(associated_strains == "full") strains
        else if(associated_strains == "half") sample(strains, 1)
        else associated_strains)
    AA.freq = get_frequencies(AA.scenario, strains,populations)
    ##If associated_SNPs is null we might be in a full G2G simulation
    associated_SNPs = if(is.null(associated_SNPs)) SNP.data[,unlist(AA.scenario$associated_SNPs), drop = F] 
      else associated_SNPs
    AA.data = generate_AAs(AA.freq, AA.scenario, study_design, associated_SNPs)
    colnames(AA.data) <- if("id" %in% colnames(AA.scenario)) unlist(AA.scenario$id)
    rownames(AA.freq) <- if("id" %in% colnames(AA.scenario)) unlist(AA.scenario$id)
    list(`data` = AA.data, `freq` = AA.freq)})
  data = do.call(cbind, lapply(AA, function(aa) aa$data))
  freq = do.call(rbind, lapply(AA, function(aa) aa$freq))
  list(`data` = data, `freq` = freq)}

get_SNP <- function(study_design, SNP.scenarios) {
  populations = levels(study_design$Population)
  strains = levels(study_design$Strain)
  SNP = lapply(1:nrow(SNP.scenarios), function(scenario_num) {
    SNP.scenario = SNP.scenarios[scenario_num,, drop = F]
    SNP.freq = get_frequencies(SNP.scenario, populations, strains)
    SNP.data = generate_SNPs(SNP.freq, study_design)
    colnames(SNP.data) <- if("id" %in% colnames(SNP.scenario)) unlist(SNP.scenario$id)
    rownames(SNP.freq) <- if("id" %in% colnames(SNP.scenario)) unlist(SNP.scenario$id)
    list(`data` = SNP.data, `freq` = SNP.freq)})
  data = do.call(cbind, lapply(SNP, function(snp) snp$data))
  freq = do.call(rbind, lapply(SNP, function(snp) snp$freq))
  list(`data` = data, `freq` = freq)}

get_frequencies <- function(scenario, main_subgroup, secondary_subgroup) {
  nb_main_subgroup = length(main_subgroup)
  nb_secondary_subgroup = length(secondary_subgroup)
  ##Things compared are on list, that is not natural!!!!
  strat = {if(is.na(scenario$Stratified) | scenario$Stratified == "no") rep(0, nb_main_subgroup) 
    else if(scenario$Stratified == "full" | scenario$Stratified == "yes") sample(1:nb_main_subgroup, nb_main_subgroup)
    else rank(unlist(scenario$Stratified))} ###Here we will use user defined order
  ###Here we make a filter that will later remove unstratified group, this act on seconday_subgroup
  p_strat = if(!is.na(scenario$Partial_Stratification)) rep(strat, each = nb_secondary_subgroup) * rep(as.numeric(secondary_subgroup == scenario$Partial_Stratification), nb_secondary_subgroup)
  
  bias = {if(is.na(scenario$Biased) | scenario$Biased == "no") rep(0, nb_main_subgroup)
    else if(scenario$Biased == "full" | scenario$Biased == "yes") sample(1:nb_secondary_subgroup, nb_secondary_subgroup)
    else rank(unlist(scenario$Biased))} ###Here we will use user defined order
  ###Here we make a filter that will later remove unbiased group, this act on seconday_subgroup
  p_bias = if(!is.na(scenario$Partial_Bias)) rep(bias, nb_main_subgroup) * rep(as.numeric(main_subgroup == scenario$Partial_Bias), each = nb_main_subgroup)
  
  #F_strat = if(length(unique(strat)) >1) sort(get_AF(RF, scenario$fst_strat, length(unique(strat))), decreasing = TRUE)[strat] else rep(get_AF(RF, 0.2,1), nb_main_subgroup)
  ###Bias first then strat
  
  lower = 0.1
  upper = 0.4
  freq = if(!is.null(p_strat)) {
    t(replicate(scenario$size, {
      RF = runif(1, lower, upper)
      F_bias = if(length(unique(bias)) >1)  
        sort(get_AF(RF, scenario$fst_bias, length(unique(bias))), decreasing = TRUE)[bias] 
      else 
        rep(RF, nb_secondary_subgroup) 
      ##IF pstrat case ! => I had to change things from rep(F_bias, each = nb_main_subgroup) to rep(F_bias, nb_main_subgroup), there mus be a reason!
      F_strat =  if(length(unique(strat)) >1) 
        c(sapply(F_bias, function(RF) sort(get_AF(RF, scenario$fst_strat, length(unique(strat))), decreasing = TRUE)[strat])) 
      else 
        rep(F_bias, nb_main_subgroup)
      F_strat[!is.null(p_strat) & p_strat == 0] = rep(F_bias, nb_main_subgroup)[p_strat == 0]
      F_strat}))
  }
  ###Strat first then bias
  else if(!is.null(p_bias) ||  (length(unique(strat)) >1)){
    t(replicate(scenario$size, {
      RF = if(!is.na(scenario$fst_strat) & length(unique(strat)) == 1){
        get_AF(runif(1, lower, upper), fst = scenario$fst_strat, nb = 1)}
      else {
        runif(1, lower, upper)}
      
      F_strat = if(length(unique(strat)) >1) {
        sort(get_AF(RF, scenario$fst_strat, length(unique(strat))), decreasing = TRUE)[strat] }
      else {
        rep(RF, nb_main_subgroup)}
      
      F_bias =  if(length(unique(bias)) >1) {
        c(sapply(F_strat, function(RF) sort(get_AF(RF, scenario$fst_bias, length(unique(bias))), decreasing = TRUE)[bias]))}
      else {
        rep(F_strat, each = nb_secondary_subgroup)}
      
      F_bias[!is.null(p_bias) & p_bias == 0] = rep(F_strat, each = nb_secondary_subgroup)[p_bias == 0]
      F_bias}))
  }
  else 	{
    do.call(rbind, replicate(scenario$size, rep(get_AF(runif(1, lower, upper), 0.02), nb_secondary_subgroup * nb_main_subgroup), simplify = F))}
  colnames(freq) <- paste0("F_", rep(main_subgroup, each = nb_secondary_subgroup),".", secondary_subgroup)
  freq}

#Get alternate frequency according to Nicholas model
get_AF <- function(allele, fst, nb = 1) {
  if(is.na(fst)) stop("Specify fst")
  s1 = allele*(1-fst)/fst
  s2 = (1-allele)*(1-fst)/fst
  get_ref <- function(ref_old=c()){
    ref_new = rbeta(n = 1, shape1 = s1,shape2 = s2)
    if(length(ref_old) == nb)
      ref_old
    else if(ref_new>0.05 && ref_new<0.95)
      get_ref(c(ref_new, ref_old))
    else
      get_ref(ref_old)}
  get_ref()
  #rbeta(n = nb, shape1 = s1,shape2 = s2)
}

generate_SNPs <- function(SNP.freq, study_design) {
  nb_strains = length(levels(study_design$Strain))
  nb_populations = length(levels(study_design$Population))
  #if(trace == TRUE) print(paste(Sys.time(),"Generating SNPs dsitribution with viral properties", sep=" : "))
  #nb = summarize(group_by(study_design,Population, Strain), nb = n())$nb
  nb = unname(unlist(to_pop_structure(study_design)[,"nb"]))
  SNP = unlist(lapply(1:nrow(SNP.freq), function(snp_num) {
    p = as.data.frame(lapply(SNP.freq[snp_num,], function(alternate_allele) {c((1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2)}))    
    unlist(lapply(1:nb_populations, function(population_num) {
      iterator = ((population_num - 1 ) * nb_strains) + 1
      unlist(lapply(iterator:(iterator + nb_strains - 1), function(selector) {
        if(nb[selector]>0) sample(c(0,1,2), size = nb[selector], prob = p[,selector], replace = TRUE) })) })) }))
  matrix(data = SNP ,nrow = nrow(study_design), ncol = nrow(SNP.freq), dimnames = list(rownames(study_design), rownames(SNP.freq)))}

generate_AAs <- function(AA.freq, scenario, study_design, associated_SNPs) {
  strains = levels(study_design$Strain)
  populations = levels(study_design$Population)
  #if(trace == TRUE) print(paste(Sys.time(),"Generating Viral output for", scenario$id_tag, sep=" : "))
  res = sapply(1:nrow(AA.freq), function(aa_num) {
    unlist(lapply(populations, function(population) {
      unlist(lapply(strains, function(strain) {
        filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
        if(sum(filter)>0) {
          Y_ass =	if(!is.null(associated_SNPs) && population %in% unlist(scenario$Associated_Populations) && strain %in% unlist(scenario$Associated_Strains)) {
            z = associated_SNPs[which(filter),,drop=F]%*%rep(scenario$beta, ncol(associated_SNPs))
            pr = 1/(1+exp(-z))
            unlist(lapply(pr -0.5, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))} 
          else numeric(sum(filter))
          AF = AA.freq[aa_num,paste0("F_",strain ,".",population)]
          Y = Y_ass + sample(0:1, size = sum(filter), prob = c(1 - AF,AF), replace = TRUE)
          Y[Y>1] = 1
          Y	}})) })) })
  rownames(res) <- paste0(rownames(study_design),"_",study_design[,"Population"],"_",study_design[,"Strain"]) 
  res}

###G2G single setup methods
analyse_G2G_setup <- function(SNP, AA, study_design, host_pc=NULL) {
  if(trace == TRUE) print(paste(Sys.time(),"Computing GLM", sep=" : "))
  pvalues = if(is.null(host_pc)) {
    list(`WO_correction` = unlist(lapply(1:ncol(AA), function(index) coef(summary(glm(AA[,index]~SNP[,index], family = binomial)))[,4][2])),
         `W_human_groups` = if(length(levels(study_design[,"Population"])) > 1) unlist(lapply(1:ncol(AA), function(index) coef(summary(glm(AA[,index]~SNP[,index]+study_design[,"Population"], family = binomial)))[,4][2])),
         `W_strain_groups` = if(length(levels(study_design[,"Strain"])) > 1) unlist(lapply(1:ncol(AA), function(index) coef(summary(glm(AA[,index]~SNP[,index]+study_design[,"Strain"], family = binomial)))[,4][2])),
         `W_both` = if(length(levels(study_design[,"Population"])) > 1 && length(levels(study_design[,"Strain"])) > 1) unlist(lapply(1:ncol(AA), function(index) coef(summary(glm(AA[,index]~SNP[,index]+study_design[,"Population"]+study_design[,"Strain"], family = binomial)))[,4][2])))
    
  } else {
    host_pc=host_pc$x[,1:5]
    list(`WO_correction` = unlist(lapply(1:ncol(AA), function(index) coef(summary(glm(AA[,index]~SNP[,index]), family = binomial))[,4][2])),
         `W_human_PCs` = if(length(levels(study_design[,"Population"])) > 1) unlist(lapply(1:ncol(AA), function(index) coef(summary(glm(AA[,index]~SNP[,index]+host_pc, family = binomial)))[,4][2])),
         `W_strain_groups` = if(length(levels(study_design[,"Strain"])) > 1) unlist(lapply(1:ncol(AA), function(index) coef(summary(glm(AA[,index]~SNP[,index]+study_design[,"Strain"], family = binomial)))[,4][2])),
         `W_both` = if(length(levels(study_design[,"Population"])) > 1 && length(levels(study_design[,"Strain"])) > 1) unlist(lapply(1:ncol(AA), function(index) coef(summary(glm(AA[,index]~SNP[,index]+host_pc+study_design[,"Strain"], family = binomial)))[,4][2])))
  }
  if(trace == TRUE) print(paste(Sys.time(),"Parsing values", sep=" : "))   
  pvalues = parse_pvalues(as.data.frame(Filter(length, pvalues)), threshold)
  list(`pvalues` = pvalues, `host_pc` = host_pc)}

###G2G full scenario methods
generate_id_tag <- function(fst_strat=NA, partial_strat=NA, fst_bias=NA, partial_bias=NA, beta=NA) {
  round_me <- function(prefix, param) ifelse(is.na(param),"", paste0(prefix,"-",round(param, digits=2)))
  part_me <- function(param) ifelse(is.na(param),"", paste0("_P-", paste0(param, collapse = ',')))
  tag = c(paste0(round_me("Fstrat", fst_strat), part_me(partial_strat)), paste0(round_me("Fbias", fst_bias), part_me(partial_bias)), round_me("Beta",beta))
  tag = paste0(Filter(function(x) x!="", tag), collapse = "_")
  tag = if(tag=="") "homogenous" else tag}

#@Not used
generate_biological_tag <- function(size, tag_size_range, id) {
  incrementer <-function() round(runif(1, min = tag_size_range["min"], max = tag_size_range["max"]))
  get_biotag_range <- function(sequence=c()) {
    if(size - sum(sequence) < tag_size_range["max"]) c(sequence, size - sum(sequence))
    else get_biotag_range(c(sequence, incrementer()))}
  bio_tag_range = get_biotag_range()
  get_bio_tag <- function (id, num=1) {
    res = setNames(list(id[1:bio_tag_range[num]]),  paste0(id[1],"_BT_",num,"_",bio_tag_range[num]))
    id = id[(bio_tag_range[num]+1):length(id)]
    if(length(bio_tag_range) >= num + 1) c(res, get_bio_tag(id, num+1)) else res}
  get_bio_tag(id)}

get_id <- function(prefix,size) {
  id = if(Sys.getenv(prefix) == "") 0 else as.numeric(Sys.getenv(prefix))
  args = list(id + size)
  names(args) = prefix
  do.call(Sys.setenv, args)
  (id+1):(id+size)}

set_env <- function(name, value) {
  args = list(value)
  names(args) = name
  do.call(Sys.setenv, args)}

analyse_G2G <- function(data, correction, analyse, nb_cpu = 1) {
  ##data contain both AA.data, SNP.data and study_design
  attach(data)
  attach(correction)
  attach(analyse)
  
  if(trace) print(paste0(Sys.time()," : Computing PC"))
  SNP_PC = if(W_human_PC || W_both_PC) {
    SNP_PC = prcomp(SNP.data)
    SNP_PC$x[,1:ifelse(ncol(SNP_PC$x)<5, ncol(SNP_PC$x), 5)]}
  
  AA_PC = if(W_strain_PC || W_both_PC){
    AA_PC = prcomp(AA.data)
    AA_PC$x[,1:ifelse(ncol(AA_PC$x)<5, ncol(AA_PC$x), 5)]}
  
  AA_NL_PC = if(W_non_linear_PC){
    AA_NL_PC = homals(AA.data, ndims = 5)
    matrix(unlist(AA_NL_PC$loadings), nrow = nrow(SNP), ncol = 5)} 
  
  #save(AA_PC, SNP_PC, AA_NL_PC, file ="PC-savestates.RData")
  
  count_ETA <- function() {
    if(trace) print(paste0(Sys.time(),": Working on amino acids batch [", counter,"-",counter+nb_cpu,"], on a total of ", nb_aa ," amino acids, ", (counter/nb_aa) * 100, "% has been done, ETA in ",
       trunc((nb_aa - counter)/(counter/((proc.time() - ptm)["elapsed"]))) %/% 86400, " Day(s), ",
       trunc(((nb_aa - counter)/(counter/((proc.time() - ptm)["elapsed"]))) %/% 3600) %% 24, " Hour(s), ",
       round(((nb_aa - counter)/(counter/((proc.time() - ptm)["elapsed"]))) %/% 60) %% 60, " Minute(s)"))
    counter <<- counter + nb_cpu}
  
  logistic_analyse <- function() {
    if(trace) print(paste(Sys.time(),": Computing logisic regression with ",  nb_cpu, " CPU(s)"))
    analyse_AA <- function(Y) {
      count_ETA()
      c(if(WO_correction) list(`Without correction` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X)))[,4][2])),
        if(W_human_group) list(`With human group` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+study_design[,"Population"], family = binomial)))[,4][2])),
        if(W_strain_group) list(`With strain group` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+study_design[,"Strain"], family = binomial)))[,4][2])),
        if(W_both_groups) list(`With both groups` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+study_design[,"Population"]+study_design[,"Strain"], family = binomial)))[,4][2])),
        if(W_human_PC) list(`With human PC` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+SNP_PC, family = binomial)))[,4][2])),
        if(W_strain_PC) list(`With strain PC` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+AA_PC, family = binomial)))[,4][2])),
        if(W_both_PC) list(`With both PC` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+SNP_PC+AA_PC, family = binomial)))[,4][2])),
        if(W_strain_groups_human_PC) list(`W_strain_groups_human_PC` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+SNP_PC+study_design[,"Strain"], family = binomial)))[,4][2])),
        if(W_non_linear_PC) list(`With non linear PC` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+AA_NL_PC+study_design[,"Population"], family = binomial)))[,4][2])))}
    
    counter <<- 0
    nb_aa <<- ncol(AA.data)
    ptm <<- proc.time()
    cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
    res = parApply(cl,AA.data, 2, analyse_AA)
    stopCluster(cl)
    print(paste0(Sys.time(),": ", nb_aa," has been analysed by logistic regression, 100% done!"))
    
    SNPcol = rep(1:ncol(SNP.data), length(res[[1]]))
    SNP_Tag = rep(unlist(mapply(function(tag, size) rep(tag,size), SNP.scenarios$bio_tag, SNP.scenarios$size)), length(res[[1]]))
    CorrectionCol = as.factor(rep(names(res[[1]]), each = ncol(SNP.data)))
    
    res = do.call(cbind, lapply(names(res), function(aa_id) {
      do.call(rbind, unname(lapply(res[[aa_id]], function(SNP) {
        setNames(data.frame(unname(SNP)), aa_id)})))}))
    res = cbind(`SNP` = SNPcol, `SNP_Tag` = SNP_Tag, `Correction` =  CorrectionCol, res)}
  
  SKAT_analyse <- function() {
    if(trace) print(paste(Sys.time(),": Computing SkAT with ",  nb_cpu, " CPU(s)"))
    skat_it = {
      if(skat_LW) function(HO) lapply(SNP_batch, function(batch) SKAT(SNP.data[,batch], HO, kernel = "linear.weighted")$p.value)
      else if(skato_LW) function(HO) lapply(SNP_batch, function(batch) SKAT(SNP.data[,batch], HO, kernel = "linear.weighted", method = "optimal.adj")$p.value)
      else if(skat_L) function(HO) lapply(SNP_batch, function(batch) SKAT(SNP.data[,batch], HO, kernel = "linear")$p.value)
      else if(skato_L) function(HO) lapply(SNP_batch, function(batch) SKAT(SNP.data[,batch], HO, kernel = "linear", method = "optimal.adj")$p.value)}
    
    analyse_AA <- function(Y) {c(
      if(WO_correction) list(`Without correction` = skat_it(SKAT_Null_Model(Y~1, out_type = "D"))),
      if(W_human_group)  list(`With human group` = skat_it(SKAT_Null_Model(Y ~ study_design[,"Population"], out_type = "D"))),
      if(W_strain_group) list(`With strain group` = skat_it(SKAT_Null_Model(Y ~ study_design[,"Strain"], out_type = "D"))),
      if(W_both_groups) list(`With both groups` = skat_it(SKAT_Null_Model(Y ~ study_design[,"Population"]+study_design[,"Strain"], out_type = "D"))),
      if(W_human_PC) list(`With human PC` = skat_it(SKAT_Null_Model(Y ~ SNP_PC, out_type = "D"))),
      if(W_strain_PC) list(`With strain PC` = skat_it(SKAT_Null_Model(Y ~ AA_PC, out_type = "D"))),
      if(W_both_PC) list(`With both PC` = skat_it(SKAT_Null_Model(Y ~ AA_PC + SNP_PC, out_type = "D"))),
      if(W_non_linear_PC) list(`With non linear PC` = skat_it(SKAT_Null_Model(Y ~ AA_NL_PC, out_type = "D"))))}
    
    flip_SNP <- function(SNP){
      nsample = nrow(SNP)
      apply(SNP,2, function(S) { 
        if(sum(S)/nsample > 1) 
          unlist(lapply(S, function(s) {
            if(s==0) 2 
            else if(s==2) 0 
            else 1}))
        else S })}
    
    #SNP = flip_SNP(SNP.data)
    bio_tag_key = unique(SNP.scenarios$bio_tag)
    SNP_batch = setNames(lapply(bio_tag_key, function(tag) unlist(filter(SNP.scenarios, bio_tag == tag)$id)), bio_tag_key)
    
    cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
    res = parApply(cl,AA.data, 2, analyse_AA)     #res = apply(AA, 2, analyse_AA)
    stopCluster(cl)
    
    SNP_Tag = rep(names(SNP_batch), length(res[[1]]))
    correction_col = as.factor(rep(names(res[[1]]),  each = length(SNP_batch)))
    res = as.data.frame(do.call(cbind,lapply(res, function(aa) {
      do.call(rbind, lapply(aa, function(correction) {
        do.call(rbind, lapply(correction, function(SNP_tag) {
          SNP_tag}))}))})))
    colnames(res) <- colnames(AA.data)
    rownames(res) <- NULL
    res = as.data.frame(cbind(`SNP_Tag` = SNP_Tag, `Correction` = correction_col, res))}
  
  GT_analyse <- function() {
    if(trace) print(paste(Sys.time(),": Computing global test with ",  nb_cpu, " CPU(s)"))
    analyse_AA <- function(Y) {c(
      if(WO_correction) list(`Without correction` = lapply(SNP_batch, function(batch) as.numeric(result(gt(Y, SNP.data[,batch]))["p-value"]))),
      if(W_human_group) list(`With human group` = lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ study_design[,"Population"], SNP.data[,batch]))["p-value"]))),
      if(W_strain_group) list(`With strain group` = lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ study_design[,"Strain"], SNP.data[,batch]))["p-value"]))),
      if(W_both_groups) list(`With both groups` = lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ study_design[,"Population"]+study_design[,"Strain"], SNP.data[,batch]))["p-value"]))),
      if(W_human_PC) list(`With human PC` = lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ SNP_PC, SNP.data[,batch]))["p-value"]))),
      if(W_strain_PC) list(`With strain PC` = lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ AA_PC, SNP.data[,batch]))["p-value"]))),
      if(W_both_PC) list(`With both PC` = lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ SNP_PC + AA_PC, SNP.data[,batch]))["p-value"]))),
      if(W_non_linear_PC) list(`With non linear PC` = lapply(SNP_batch, function(batch) as.numeric(result(gt(Y + AA_NL_PC, SNP.data[,batch]))["p-value"]))))}
    
    bio_tag_key = unique(SNP.scenarios$bio_tag)
    SNP_batch = setNames(lapply(bio_tag_key, function(tag) unlist(filter(SNP.scenarios, bio_tag == tag)$id)), bio_tag_key)
    
    cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
    res = parApply(cl,AA.data, 2, analyse_AA)     #res = apply(AA, 2, analyse_AA)
    stopCluster(cl)
    
    SNP_Tag = rep(names(SNP_batch), length(res[[1]]))
    correction_col = as.factor(rep(names(res[[1]]),  each = length(SNP_batch)))
    res = as.data.frame(do.call(cbind,lapply(res, function(aa) {
      do.call(rbind, lapply(aa, function(correction) {
        do.call(rbind, lapply(correction, function(SNP_tag) {
          SNP_tag}))}))})))
    colnames(res) <- colnames(AA.data)
    rownames(res) <- NULL
    res = as.data.frame(cbind(`SNP_Tag` = SNP_Tag, `Correction` = correction_col, res))}
  
  G2_analyse <- function() {
    bio_tag_key = unique(SNP.scenarios$bio_tag)
    SNP_batch_id = setNames(lapply(bio_tag_key, function(tag) unlist(filter(SNP.scenarios, bio_tag == tag)$id)), bio_tag_key)
    bio_tag_key = unique(AA.scenarios$bio_tag)
    AA_batch_id = setNames(lapply(bio_tag_key, function(tag) unlist(filter(AA.scenarios, bio_tag == tag)$id)), bio_tag_key)
    if(trace) print(paste(Sys.time(),": Computing G2 with ",  nb_cpu, " CPU(s)"))
    analyse_AA <- function(Y) {
      lapply(SNP_batch_id, function(X) G2(AA.data[,Y], SNP.data[,X],nperm = 2500)$std_pval)
    }
    cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
    res = parLapply(cl,AA_batch_id, analyse_AA)     #res = apply(AA, 2, analyse_AA)
    stopCluster(cl)
    res = do.call(cbind, lapply(res, function(AA) do.call(rbind, AA)))
    colnames(res) <- unique(AA.scenarios$bio_tag)
    cbind("SNP_Tag" = paste0("SNP_", unique(SNP.scenarios$bio_tag)), as.data.frame(res))}
  
  res = c(
    if(logistic) list(`logistic` = logistic_analyse()),
    if(gt) list(`gt` = GT_analyse()),
    if(G2_analysis) list(`G2_analysis` = G2_analyse()),
    if(skat_LW | skat_L | skato_LW | skato_L) `Skat` = SKAT_analyse())
  
  #save(AA_PC, SNP_PC, res, AA.data,SNP.data, AA.scenarios, SNP.scenarios, file ="res-savestates.RData")
  if(trace) print(paste(Sys.time()," : Analysis took ", (proc.time() - ptm)["elapsed"]/60, "minutes"))
  
  detach(data)
  detach(correction)
  detach(analyse)
  list(`data` = data, `results` = res, `PC` = SNP_PC)}