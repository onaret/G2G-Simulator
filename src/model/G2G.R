###G2G common methods
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
    AA.freq = get_frequencies(AA.scenario, strains, populations)
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
  else {
    do.call(rbind, replicate(scenario$size, 
                             rep(get_AF(runif(1, lower, upper), 0.02), 
                                 nb_secondary_subgroup * nb_main_subgroup), 
                             simplify = F))
  }
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
