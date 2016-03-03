scenario_viral <- function(sample_size, strains, nb_pop, SNPs, viral_aa=1, fcoeff_pop, fcoeff_vir, vir_bias, pop_bias, beta) {
  params = parse_params(sample_size, strains, nb_pop, viral_aa, fcoeff_pop,fcoeff_vir, vir_bias, pop_bias, beta)
  scenarios = unique(SNPs)
  SNP_params = parse_SNP_parameters(SNPs, params)
  scenarios["size"] = group_size(group_by(SNP_params,S_Stratified,S_Biased,Y_Stratified,Y_Biased,Viral_Association,R))
  different_scenarios = unique(select(SNP_params, -starts_with("S_P"), -starts_with("Y_P")))
  study_design = get_study_design(params)
  SNPs = generate_SNPs_with_viral(SNP_params, study_design, params)
  Y = get_viral_output(study_design, SNPs, SNP_params, params)
  pvalues = analyse_viral(SNPs,Y, study_design)
  summary_sim = summary_sim(pvalues, SNP_params)
  res = list(#`SNPs` = SNPs,
    `study_design` = study_design,
    `SNP_params` = SNP_params,
    `pvalues` = pvalues,
    `summary_sim` = summary_sim,
    `params` = params,
    `scenarios` = scenarios)
  write(res, tag = paste("size", sample_size,"fc_pop", round(fcoeff_pop, digits = 3), "beta", round(beta, digits = 3),sep = "-"))
  res}

parse_params <- function(sample_size, strains, nb_pop, viral_aa=1, fcoeff_pop,fcoeff_vir, vir_bias, pop_bias=NULL, beta) {
  pops = paste0("P", 1:nb_pop)
  nb_strains = length(strains)
  list(`fcoeff_pop` = fcoeff_pop, `fcoeff_vir` = fcoeff_vir, `vir_bias` = vir_bias, `pop_bias` = pop_bias, 
      `viral_aa` = viral_aa, `nb_strains` = nb_strains, `strains` = strains, `nb_pop` = nb_pop, `pops` = pops, `sample_size` = sample_size, `beta` = beta)}

get_study_design <- function(params) {
  ####Case is not usefull in vir study
  CC_structure = get_CC_structure(params$sample_size, params$nb_pop)
  pop_structure = get_pop_structure(params$sample_size, params$nb_pop)
  viral_pop_structure = get_viral_structure(params$sample_size, params$strains)
  study_design = data.frame(`Population` = pop_structure, `Case` = CC_structure, `Strain` = viral_pop_structure)
  study_design = arrange(study_design, Population, Strain, Case)
  rownames(study_design) <- get_samples_name_from_SD(study_design)
  study_design}

get_pop_structure <- function(sample_size, nb_pop) paste0("P", rep(x = 1:nb_pop, times = 1, each = sample_size/nb_pop))
get_CC_structure <- function(sample_size, nb_pop) c(replicate(nb_pop, c(rep(TRUE, sample_size/(nb_pop*2)), rep(FALSE, sample_size/(nb_pop*2)))))
get_viral_structure <- function(sample_size, strains) sample(strains, size = sample_size, replace = T)

get_samples_name_from_SD <- function(study_design) {
  paste0(1:nrow(study_design),"_", study_design[,1],"_", ifelse(study_design[,2], "Case", "Control"),"_",study_design[,3])}

parse_SNP_parameters <- function(SNPs, params) {
  ####TODO : remove `Viral Association` when CC study
  print(paste(Sys.time(),"Generating SND ad Y frequencies", sep=" : "))
  lvl = c("yes","no","half")
  prototype = data.frame(`S_Stratified` = factor(levels = lvl), `S_Biased` = factor(levels = lvl), `Y_Stratified` = factor(levels =lvl), `Y_Biased` = factor(levels = lvl), `Viral_Association` = character(), `R` = numeric())
  SNPs = rbind.fill(SNPs,prototype)
  SNPs[,"Y_Stratified"][SNPs[,"Y_Biased"] == "half"] = "yes"
  SNPs[,"Y_Stratified"][is.na(SNPs[,"Y_Stratified"])] = "no"
  
  SNPs[,"Y_Biased"][SNPs[,"Y_Stratified"] == "half"] = "yes"
  SNPs[,"Y_Biased"][is.na(SNPs[,"Y_Biased"])] = "no"
  
  SNPs[,"S_Stratified"][SNPs[,"S_Biased"] == "half"] = "yes"
  SNPs[,"S_Stratified"][is.na(SNPs[,"S_Stratified"])] = "no"
  
  SNPs[,"S_Biased"][SNPs[,"S_Stratified"] == "half"] = "yes"
  SNPs[,"S_Biased"][is.na(SNPs[,"S_Biased"])] = "no"
  
  SNPs["Causal"]=!is.na(SNPs[,"R"])
  y_freq = t(data.frame(apply(select(SNPs, -R, -Viral_Association, -Causal),1, function(snp) get_AFs(params, snp['Y_Stratified'], params$fcoeff_vir, snp['Y_Biased'], params$vir_bias))))
  s_freq = t(data.frame(apply(select(SNPs, -R, -Viral_Association, -Causal),1, function(snp) get_AFs(params, snp['S_Stratified'], params$fcoeff_pop, snp['S_Biased'], params$pop_bias))))
  colnames(s_freq) <- paste0("S_", "P", as.vector(t(replicate(params$nb_strains, 1:params$nb_pop))),".", params$strains)
  colnames(y_freq) <- paste0("Y_", "P",1:params$nb_pop,".", as.vector(t(replicate(params$nb_pop, params$strains))))
  SNPs = cbind(s_freq, y_freq, SNPs)
  rownames(SNPs) <- paste0('Snp_S_', 1:nrow(SNPs), '_N') 
  SNPs}

get_AFs <- function(params, stratified, fcoeff=0, biased, bias =0) {
  get_AF <- function(allele, fcoeff) {
    s1 = allele*(1-fcoeff)/fcoeff
    s2 = (1-allele)*(1-fcoeff)/fcoeff
    rbeta(n = 1, shape1 = s1,shape2 = s2)}
  AFs = {
    if(biased == "half") {
      t = sort(replicate(3,get_AF(runif(1, 0.1, 0.5), fcoeff)), decreasing = TRUE)
      c(t,tail(t, n=1)) }
    else if(stratified == "half") {
      t = sort(replicate(3,get_AF(runif(1, 0.1, 0.5), fcoeff)), decreasing = TRUE)
      c(head(t, n=1), tail(t, n=1), tail(head(t, n=2), n=1), tail(t, n=1)) }
    else if(stratified == "no" && biased == "no") rep(get_AF(runif(1, 0.1, 0.5), fcoeff),params$nb_pop * params$nb_strains)
    else if(stratified == "no" && biased == "yes") rep(sort(replicate(params$nb_strains, get_AF(runif(1, 0.1, 0.5), bias)), decreasing = TRUE), params$nb_pop)
    else if(stratified == "yes" && biased == "no") sort(as.vector(replicate(params$nb_pop, rep(get_AF(runif(1, 0.1, 0.5), fcoeff), params$nb_strains))), decreasing = TRUE) 
    else if(stratified == "yes" && biased == "yes") {
      stratified_afs = sort(as.vector(replicate(params$nb_pop, rep(get_AF(runif(1, 0.1, 0.5), fcoeff), params$nb_strains))), decreasing = TRUE)
      st = unlist(lapply(stratified_afs, function(stratified_af) get_AF(stratified_af,  bias)))
      c(sort(st[1:2], decreasing = TRUE),  sort(st[3:4], decreasing = TRUE))} } }

generate_SNPs_with_viral <- function(SNP_params, study_design, params) {
  print(paste(Sys.time(),"Generating SNPs dsitribution with viral properties", sep=" : "))
  data = matrix(nrow = nrow(study_design), ncol = nrow(SNP_params), dimnames = list(rownames(study_design), rownames(SNP_params)))
  SNP_params = select(SNP_params, +starts_with("S_P"))
  sapply(1:nrow(SNP_params), function(snp_num) {
    p = get_SNP_probability(SNP_params[snp_num,], CC=FALSE)
    unlist(lapply(params$pops, function(population) {
      unlist(lapply(params$strains, function(strain) {
        prob = unname(unlist(select(p, contains(paste0("S_", population,"." ,strain) ))))
        sample(c(0,1,2), size = nrow(filter(study_design, Population == population, Strain == strain)), prob = prob, replace = TRUE) })) })) })}

get_viral_output <- function(study_design, SNPs, SNP_params, params) {
  print(paste(Sys.time(),"Generating Viral output", sep=" : "))
  Y = sapply(1:nrow(SNP_params), function(snp_num) {
    Y = unlist(lapply(params$pops, function(population) {
      unlist(lapply(params$strains, function(strain) {
        AF = SNP_params[snp_num,paste0("Y_", population,".",strain)]
        sample(0:1, size = nrow(filter(study_design, Population == population, Strain == strain)), prob = c(1 - AF,AF), replace = TRUE) })) }))
    if(SNP_params[snp_num,"Causal"]) {
      filter = study_design[,"Strain"] %in% unlist(strsplit(x = as.character(SNP_params[snp_num,"Viral_Association"]), split = ""))
      z = numeric(sum(filter))
      z = as.matrix(SNPs[,snp_num])[which(filter)]*params$beta
      pr = 1/(1+exp(-z))
      Y[filter] = unlist(lapply(pr, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))}
    Y})
  Y}

analyse_viral <- function(SNPs, Y, study_design) {
  print(paste(Sys.time(),"Computing PC from SNPs for population stratification", sep=" : "))
  #SNPs_PC = prcomp(SNPs, .scale = FALSE)
  #nb_PCs_strains = ifelse(ncol(SNPs_PC$x)<5, ncol(SNPs_PC$x), 5)
  print(paste(Sys.time(),"Computing PC from Y for viral stratification", sep=" : "))
  #strain_PC = prcomp(Y, .scale = FALSE)
  #nb_PCs_strains = ifelse(ncol(strain_PC$x)<5, ncol(strain_PC$x), 5)
  print(paste(Sys.time(),"Computing GLM", sep=" : "))
  WO_PC = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num])))[,4][2]))
  W_S_Pop_PC = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+study_design[,"Population"])))[,4][2]))
  W_S_Vir_PC = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+study_design[,"Strain"])))[,4][2]))
  #W_C_Pop_PC = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+SNPs_PC$x[,1:5])))[,4][2]))
  W_S_Pop_S_Vir_PC = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+study_design[,"Population"]+study_design[,"Strain"])))[,4][2]))
  #W_C_Pop_C_Vir_PC = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+SNPs_PC$x[,1:5]+strain_PC$x[,1:nb_PCs_strains])))[,4][2]))
  print(paste(Sys.time(),"Parsing values", sep=" : "))
  parse_pvalues(data.frame(WO_PC, W_S_Pop_PC, W_S_Vir_PC, W_S_Pop_S_Vir_PC, row.names = colnames(SNPs)), threshold)}
