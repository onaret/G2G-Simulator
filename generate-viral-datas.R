get_SNP_input <- function(S_Stratified = NA, S_Biased = NA,Y_Stratified = NA, Y_Biased = NA, Associated_Strains = NA, Associated_Populations = NA, size = size, Ordered_Bias = FALSE) {
  lvl = c("full","no","half", "P1", "P2", "A", "B")
  if(is.na(Associated_Strains) & !is.na(Associated_Populations)) Associated_Strains = "full"
  if(is.na(Associated_Populations) & !is.na(Associated_Strains)) Associated_Populations = "full"
  #  if(is.na(Associated_Strains) & !is.na(Associated_Populations)) throw("Define associated_Populations")
  # if(is.na(Associated_Populations) & !is.na(Associated_Strains)) throw("Define Associated_Strains")
  tbl_df(data.frame(`S_Stratified` = factor(rep(S_Stratified, size), levels = lvl), `S_Biased` = factor(S_Biased, levels = lvl), `Y_Stratified` = factor(Y_Stratified, levels = lvl), 
                    `Y_Biased` = factor(Y_Biased, levels = lvl), `Associated_Strains` = Associated_Strains, `Associated_Populations` = Associated_Populations, `Ordered_Bias` = Ordered_Bias,
                    `size` = size, stringsAsFactors = FALSE) )}

scenario_viral <- function(sample_size, nb_strains, nb_pop, scenarios, fcoeff_pop=NA, fcoeff_vir=NA, vir_bias=NA, pop_bias=NA, beta=NA, tag = NULL, get_viral = get_viral_output) {
  params = parse_params(sample_size, nb_strains, nb_pop, fcoeff_pop,fcoeff_vir, vir_bias, pop_bias, beta, scenarios)
  tag = if(!is.null(tag)) paste0("_", tag)
  scenario_name = paste0(deparse(substitute(scenarios)),tag)
  tag_me <- function(param) ifelse(is.na(param),"", paste(deparse(substitute(param)), round(param, digits=2), sep = "-"))
  tag = paste(scenario_name, tag_me(fcoeff_pop), tag_me(pop_bias), tag_me(fcoeff_vir), tag_me(vir_bias), tag_me(beta), tag_me(sample_size), sep = "-")
  
  print(paste0(Sys.time(), " : Scenario ", scenario_name))
  SNP_params = parse_SNP_parameters(scenarios, params)
  study_design = get_study_design(params, nb_strains, nb_pop)
  SNPs = generate_SNPs_with_viral(SNP_params, study_design, params)
  Y = get_viral(study_design, SNPs, SNP_params, params)
  pvalues = analyse_viral(SNPs,Y, study_design)
  summary_sim = summary_sim(pvalues, SNP_params)
  res = list(#`SNPs` = SNPs,
    `study_design` = study_design,
    `SNP_params` = SNP_params,
    `pvalues` = pvalues,
    `summary_sim` = summary_sim,
    `params` = params,
    `scenarios` = scenarios)
  write(res, tag, paste0(getwd(),"/gen-data/",scenario_name,"/" ))
  res}

parse_params <- function(sample_size, nb_strains, nb_pop, fcoeff_pop,fcoeff_vir, vir_bias, pop_bias=NULL, beta, scenarios) {
  #Can t because of later parsing of SNP, to improve
  #if(is.na(fcoeff) & scenarios[,"S_Stratified"] != no) stop("Specify fcoeff")
  #if(is.na(fcoeff_vir) & scenarios[,"Y_Stratified"] != no) stop("Specify fcoeff_vir")
  #if(is.na(pop_bias) & scenarios[,"S_Biased"] != no) stop("Specify pop_bias")
  #if(is.na(vir_bias) & scenarios[,"Y_Biased"] != no) stop("Specify vir_bias")
  nb_strains = nb_strains[1]
  pops = paste0("P", 1:nb_pop)
  strains = chartr("123456789", "ABCDEFGHI", 1:nb_strains)
  list(`fcoeff_pop` = fcoeff_pop, `fcoeff_vir` = fcoeff_vir, `vir_bias` = vir_bias, `pop_bias` = pop_bias, `nb_strains` = nb_strains, 
       `strains` = strains, `nb_pop` = nb_pop[1], `pops` = pops, `sample_size` = sample_size, `beta` = beta)}

get_study_design <- function(params, nb_strains, nb_pop) {
  ####Case is not usefull in vir study
  CC_structure = get_CC_structure(params$sample_size, nb_pop)
  pop_structure = get_pop_structure(params$sample_size, nb_pop)
  viral_pop_structure = get_viral_structure(params$sample_size, params$strains, nb_strains)
  study_design = data.frame(`Population` = pop_structure, `Case` = CC_structure, `Strain` = viral_pop_structure)
  study_design = arrange(study_design, Population, Strain, Case)
  rownames(study_design) <- get_samples_name_from_SD(study_design)
  study_design}

get_pop_structure <- function(sample_size, nb_pop) paste0("P", rep(x = 1:nb_pop, times = 1, each = sample_size/nb_pop))
get_CC_structure <- function(sample_size, nb_pop) sample(c(T, F), sample_size, T)
#get_viral_structure <- function(sample_size, strains, nb_strains) ifelse(length(nb_strains)>1, sample(strains, sample_size, T, nb_strains[2:length(nb_strains)]), sample(strains, sample_size, T))
get_viral_structure <- function(sample_size, strains, nb_strains) {
  if(length(nb_strains)>1)  sample(strains, size = sample_size, prob = nb_strains[2:length(nb_strains)],replace = T)  
  else sample(strains, size = sample_size, replace = T) } 
get_samples_name_from_SD <- function(study_design)  paste0(1:nrow(study_design),"_", study_design[,1],"_", ifelse(study_design[,2], "Case", "Control"),"_",study_design[,3])

parse_SNP_parameters <- function(scenarios, params) {
  ####TODO : remove `Viral Association` when CC study
  if(trace == TRUE) print(paste(Sys.time(),"Generating SNP and Y frequencies", sep=" : "))
  #lvl = c("full","no","half")
  #prototype = data.frame(`S_Stratified` = factor(levels = lvl), `S_Biased` = factor(levels = lvl), `Y_Stratified` = factor(levels =lvl), `Y_Biased` = factor(levels = lvl), `Associated_Strains` = character(), `R` = numeric())
  #SNPs = rbind_all(SNPs,prototype)
  scenarios[,"Y_Stratified"][scenarios[,"Y_Biased"] == "half"] = "full"
  scenarios[,"Y_Stratified"][is.na(scenarios[,"Y_Stratified"])] = "no"
  
  scenarios[,"Y_Biased"][scenarios[,"Y_Stratified"] == "half"] = "full"
  scenarios[,"Y_Biased"][is.na(scenarios[,"Y_Biased"])] = "no"
  
  scenarios[,"S_Stratified"][scenarios[,"S_Biased"] == "half"] = "full"
  scenarios[,"S_Stratified"][is.na(scenarios[,"S_Stratified"])] = "no"
  
  scenarios[,"S_Biased"][scenarios[,"S_Stratified"] == "half"] = "full"
  scenarios[,"S_Biased"][is.na(scenarios[,"S_Biased"])] = "no"
  scenarios["Causal"]=!is.na(scenarios[,"Associated_Strains"])
  
  scenarios[,"Associated_Populations"][scenarios[,"Causal"] & (is.na(scenarios[,"Associated_Populations"]) | scenarios[,"Associated_Populations"] == "full")] = list(params$pops)
  scenarios[,"Associated_Populations"][scenarios[,"Associated_Populations"] == "half"] = params$pops
  
  scenarios[,"Associated_Strains"][scenarios[,"Causal"] & (is.na(scenarios[,"Associated_Strains"]) | scenarios[,"Associated_Strains"] == "full")] = list(params$strains)
  scenarios[,"Associated_Strains"][scenarios[,"Associated_Strains"] == "half"] = sample(params$strains, nrow(scenarios[,"Associated_Strains"] == "half"), 1)
  
  #sort(rep(params$pops, sum(scenarios[,"Associated_Populations"] == "half")/params$nb_pop))
  y_freq = t(data.frame(apply(select(scenarios, -Associated_Strains,-Associated_Populations, -Causal),1, function(snp) get_AFs(params, snp['Y_Stratified'], params$fcoeff_vir, snp['Y_Biased'], params$vir_bias, snp['Ordered_Bias'], runif(1, 0.1, 0.4) ))))
  s_freq = t(data.frame(apply(select(scenarios, -Associated_Strains,-Associated_Populations, -Causal),1, function(snp) get_AFs(params, snp['S_Stratified'], params$fcoeff_pop, snp['S_Biased'], params$pop_bias, snp['Ordered_Bias'], runif(1, 0.1, 0.4) ))))
  colnames(s_freq) <- paste0("S_", "P", as.vector(t(replicate(params$nb_strains, 1:params$nb_pop))),".", params$strains)
  colnames(y_freq) <- paste0("Y_", "P",1:params$nb_pop,".", as.vector(t(replicate(params$nb_pop, params$strains))))
  scenarios = cbind(s_freq, y_freq, scenarios)
  rownames(scenarios) <- paste0('Snp_S_', 1:nrow(scenarios), '_N') 
  scenarios}

get_AFs <- function(params, stratified, fcoeff=0, biased, bias =0, order_bias, RF) {
  get_AF <- function(allele, fcoeff) {
    if(is.na(fcoeff)) stop("Specify fcoeff")
    s1 = allele*(1-fcoeff)/fcoeff
    s2 = (1-allele)*(1-fcoeff)/fcoeff
    rbeta(n = 1, shape1 = s1,shape2 = s2)}
  #  RF = runif(1, 0.1, 0.5) 
  if(stratified == "no" && biased == "no") rep(RF,params$nb_pop * params$nb_strains)
  else if(order_bias == FALSE) {
    if(stratified == "no") {
      if(biased == "full") rep(replicate(params$nb_strains, get_AF(RF, bias)), params$nb_pop)
      else {
        #as.numeric(strsplit(biased, "P")[[1]][2])
        if(biased == "P1"| biased == "A") rep(sort(replicate(params$nb_strains, get_AF(RF, bias)), decreasing = TRUE), params$nb_pop)
        else if(biased == "P2"| biased == "B") rep(sort(replicate(params$nb_strains, get_AF(RF, bias)), decreasing = FALSE), params$nb_pop) }  }
    else if(biased == "no") {
      if (stratified == "full")  as.vector(replicate(params$nb_pop, rep(get_AF(RF, fcoeff), params$nb_strains)))
      else {
        if(stratified == "P1" | stratified == "A" )  sort(as.vector(replicate(params$nb_pop, rep(get_AF(RF, fcoeff), params$nb_strains))), decreasing = TRUE)
        else if(stratified == "P2" | stratified == "B")  sort(as.vector(replicate(params$nb_pop, rep(get_AF(RF, fcoeff), params$nb_strains))), decreasing = FALSE)  }  }
    else if (stratified == "full" & biased == "full") {
      stratified_afs = as.vector(replicate(params$nb_pop, rep(get_AF(RF, fcoeff), params$nb_strains)))
      unlist(lapply(stratified_afs, function(stratified_af) get_AF(stratified_af,  bias)))}
    else {
      stratified_afs = sort(as.vector(replicate(params$nb_pop, rep(get_AF(RF, fcoeff), params$nb_strains))), decreasing = TRUE)
      st = unlist(lapply(stratified_afs, function(stratified_af) get_AF(stratified_af,  bias)))
      c(sort(st[1:2], decreasing = TRUE),  sort(st[3:4], decreasing = TRUE))} }
  else {
    if(biased == "half") {
      t = sort(replicate(3,get_AF(RF, bias)), decreasing = TRUE)
      c(t,tail(t, n=1)) }
    else if(stratified == "half") {
      t = sort(replicate(3,get_AF(RF, fcoeff)), decreasing = TRUE)
      c(head(t, n=1), tail(t, n=1), tail(head(t, n=2), n=1), tail(t, n=1)) }
    else if(stratified == "no" && biased == "full") rep(sort(replicate(params$nb_strains, get_AF(RF, bias)), decreasing = TRUE), params$nb_pop)
    else if(stratified == "full" && biased == "no") sort(as.vector(replicate(params$nb_pop, rep(get_AF(RF, fcoeff), params$nb_strains))), decreasing = TRUE) 
    else if(stratified == "full" && biased == "full") {
      stratified_afs = sort(as.vector(replicate(params$nb_pop, rep(get_AF(RF, fcoeff), params$nb_strains))), decreasing = TRUE)
      st = unlist(lapply(stratified_afs, function(stratified_af) get_AF(stratified_af,  bias)))
      c(sort(st[1:2], decreasing = TRUE),  sort(st[3:4], decreasing = TRUE))}
  }
}

generate_SNPs_with_viral <- function(SNP_params, study_design, params) {
  if(trace == TRUE) print(paste(Sys.time(),"Generating SNPs dsitribution with viral properties", sep=" : "))
  data = matrix(nrow = nrow(study_design), ncol = nrow(SNP_params), dimnames = list(rownames(study_design), rownames(SNP_params)))
  SNP_params = select(SNP_params, +starts_with("S_P"))
  sapply(1:nrow(SNP_params), function(snp_num) {
    p = get_SNP_probability(SNP_params[snp_num,], CC=FALSE)
    unlist(lapply(params$pops, function(population) {
      unlist(lapply(params$strains, function(strain) {
        prob = unname(unlist(select(p, contains(paste0("S_", population,"." ,strain) ))))
        sample(c(0,1,2), size = nrow(filter(study_design, Population == population, Strain == strain)), prob = prob, replace = TRUE) })) })) })}

get_viral_output_sp <- function(study_design, SNPs, SNP_params, params) {
  if(trace == TRUE) print(paste(Sys.time(),"Generating Viral output", sep=" : "))
  sapply(1:nrow(SNP_params), function(snp_num) {
    Y = unlist(lapply(params$pops, function(population) {
      unlist(lapply(params$strains, function(strain) {
        filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
        AF = SNP_params[snp_num,paste0("Y_", population,".",strain)]
        Y = sample(0:1, size = sum(filter), prob = c(1 - AF,AF), replace = TRUE)
        if(SNP_params[snp_num,"Causal"] && population %in% unlist(SNP_params[snp_num,"Associated_Populations"]) && strain %in% unlist(SNP_params[snp_num,"Associated_Strains"])) {
          z = as.matrix(SNPs[,snp_num])[which(filter)]*params$beta
          pr = 1/(1+exp(-z))
          Y = Y + unlist(lapply(pr -0.5, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))
          Y[Y>1] = 1
        }
        Y })) })) 
    Y[study_design[,"Population"] == "P2"] = Y[study_design[,"Population"] == "P1"]
    Y
  })}

get_viral_output_sp2 <- function(study_design, SNPs, SNP_params, params) {
  if(trace == TRUE) print(paste(Sys.time(),"Generating Viral output", sep=" : "))
  sapply(1:nrow(SNP_params), function(snp_num) {
    Y = unlist(lapply(params$pops, function(population) {
      unlist(lapply(params$strains, function(strain) {
        filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
        AF = SNP_params[snp_num,paste0("Y_", population,".",strain)]
        Y = sample(0:1, size = sum(filter), prob = c(1 - AF,AF), replace = TRUE)
        if(SNP_params[snp_num,"Causal"] && population %in% unlist(SNP_params[snp_num,"Associated_Populations"]) && strain %in% unlist(SNP_params[snp_num,"Associated_Strains"])) {
          z = as.matrix(SNPs[,snp_num])[which(filter)]*params$beta
          pr = 1/(1+exp(-z))
          Y = Y + unlist(lapply(pr -0.5, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))
          Y[Y>1] = 1
        }
        Y })) })) 
    Y[study_design[,"Population"] == "P2"] = Y[study_design[,"Population"] == "P1"]
    
    Y2 = unlist(lapply(params$pops, function(population) {
      unlist(lapply(params$strains, function(strain) {
        filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
        AF = SNP_params[snp_num,paste0("Y_", population,".",strain)]
        Y2 = numeric(sum(filter))
        if(SNP_params[snp_num,"Causal"] && population %in% unlist(SNP_params[snp_num,"Associated_Populations"]) && strain %in% unlist(SNP_params[snp_num,"Associated_Strains"])) {
          z = as.matrix(SNPs[,snp_num])[which(filter)]*params$beta
          pr = 1/(1+exp(-z))
          Y2 = unlist(lapply(pr -0.5, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))}
        Y2})) }))
    
    Y[study_design[,"Population"] == "P2"] = Y[study_design[,"Population"] == "P2"] + Y2[study_design[,"Population"] == "P2"]
    Y[Y>1] = 1
    Y
  })}

get_viral_output_sp3 <- function(study_design, SNPs, SNP_params, params) {
  if(trace == TRUE) print(paste(Sys.time(),"Generating Viral output", sep=" : "))
  sapply(1:nrow(SNP_params), function(snp_num) {
    Y = unlist(lapply(params$pops, function(population) {
      unlist(lapply(params$strains, function(strain) {
        filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
        AF = SNP_params[snp_num,paste0("Y_", population,".",strain)]
        Y = sample(0:1, size = sum(filter), prob = c(1 - AF,AF), replace = TRUE)
        if(SNP_params[snp_num,"Causal"] && population %in% unlist(SNP_params[snp_num,"Associated_Populations"]) && strain %in% unlist(SNP_params[snp_num,"Associated_Strains"])) {
          z = as.matrix(SNPs[,snp_num])[which(filter)]*params$beta
          pr = 1/(1+exp(-z))
          Y = Y + unlist(lapply(pr -0.5, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))
          Y[Y>1] = 1
        }
        Y })) }))
    
    Y2 = unlist(lapply(params$pops, function(population) {
      unlist(lapply(params$strains, function(strain) {
        filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
        AF = SNP_params[snp_num,paste0("Y_", population,".",strain)]
        Y2 = numeric(sum(filter))
        if(SNP_params[snp_num,"Causal"] && population %in% unlist(SNP_params[snp_num,"Associated_Populations"]) && strain %in% unlist(SNP_params[snp_num,"Associated_Strains"])) {
          z = as.matrix(SNPs[,snp_num])[which(filter)]*params$beta
          pr = 1/(1+exp(-z))
          Y2 = unlist(lapply(pr -0.5, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))}
        Y2})) }))
    
    Y[study_design[,"Population"] == "P2"] = Y[study_design[,"Population"] == "P1"] + Y2[study_design[,"Population"] == "P2"]
    Y[Y>1] = 1
    Y[sample(which(Y[study_design[,"Population"] == "P1"] == 1), size = sum(Y2[study_design[,"Population"] == "P2"]))] = 0
    Y
  })}


get_viral_output <- function(study_design, SNPs, SNP_params, params) {
  if(trace == TRUE) print(paste(Sys.time(),"Generating Viral output", sep=" : "))
  sapply(1:nrow(SNP_params), function(snp_num) {
    unlist(lapply(params$pops, function(population) {
      unlist(lapply(params$strains, function(strain) {
        filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
        AF = SNP_params[snp_num,paste0("Y_", population,".",strain)]
        Y = sample(0:1, size = sum(filter), prob = c(1 - AF,AF), replace = TRUE)
        if(SNP_params[snp_num,"Causal"] && population %in% unlist(SNP_params[snp_num,"Associated_Populations"]) && strain %in% unlist(SNP_params[snp_num,"Associated_Strains"])) {
          z = as.matrix(SNPs[,snp_num])[which(filter)]*params$beta
          pr = 1/(1+exp(-z))
          Y = Y + unlist(lapply(pr -0.5, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))
          Y[Y>1] = 1
        }
        Y })) })) })}

analyse_viral <- function(SNPs, Y, study_design) {
  #if(trace == TRUE) print(paste(Sys.time(),"Computing PC from SNPs for population stratification", sep=" : "))
  #SNPs_PC = prcomp(SNPs, .scale = FALSE)
  #nb_PCs_strains = ifelse(ncol(SNPs_PC$x)<5, ncol(SNPs_PC$x), 5)
  #if(trace == TRUE) print(paste(Sys.time(),"Computing PC from Y for viral stratification", sep=" : "))
  #strain_PC = prcomp(Y, .scale = FALSE)
  #nb_PCs_strains = ifelse(ncol(strain_PC$x)<5, ncol(strain_PC$x), 5)
  if(trace == TRUE) print(paste(Sys.time(),"Computing GLM", sep=" : "))
  WO_PC = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num])))[,4][2]))
  W_S_Pop_PC = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+study_design[,"Population"])))[,4][2]))
  W_S_Vir_PC = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+study_design[,"Strain"])))[,4][2]))
  #W_C_Pop_PC = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+SNPs_PC$x[,1:5])))[,4][2]))
  W_S_Pop_S_Vir_PC = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+study_design[,"Population"]+study_design[,"Strain"])))[,4][2]))
  #W_C_Pop_C_Vir_PC = unlist(lapply(1:ncol(Y), function(SNP_Y_num) coef(summary(glm(Y[,SNP_Y_num]~SNPs[,SNP_Y_num]+SNPs_PC$x[,1:5]+strain_PC$x[,1:nb_PCs_strains])))[,4][2]))
  if(trace == TRUE) print(paste(Sys.time(),"Parsing values", sep=" : "))   
  parse_pvalues(data.frame(WO_PC, W_S_Pop_PC, W_S_Vir_PC, W_S_Pop_S_Vir_PC, row.names = colnames(SNPs)), threshold)}