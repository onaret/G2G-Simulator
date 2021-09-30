### G2G common methods
to_pop_structure <- function(study_design) {
  do.call(rbind, lapply(levels(study_design$Host), function(host) {
    do.call(rbind, lapply(levels(study_design$Pathogen), function(pathogen) { 
      tibble(pathogen, host, `nb` = sum(study_design[,"Pathogen"] == pathogen & study_design[,"Host"] == host))
    })) 
  }))
}

get_study_design <- function(structure) {
  structure = as.data.frame(structure)
  study_design = do.call(rbind, lapply(1:ncol(structure), function(pop_num) {
    do.call(rbind, lapply(1:nrow(structure), function(pathogen_num) {
      if(structure[pathogen_num, pop_num] > 0) {
        data.frame(Host = paste0("P",pop_num), 
                   Pathogen = rep(chartr("123456789", "ABCDEFGHI", pathogen_num), structure[pathogen_num, pop_num]) )
      }
    }))
  }))
  study_design$Host = as.factor(study_design$Host)
  study_design$Pathogen = as.factor(study_design$Pathogen)
  study_design
}


### G2G single setup methods
test_G2G_setup <- function(study_design, scenario, fst_host_strat=NA, fst_host_bias=NA, fst_pathogen_strat=NA, fst_pathogen_bias=NA, 
                           tag = "unnamed", sup_SNP_for_PC = NA) {
  #scenario_name = paste0(deparse(substitute(scenario)), if(!is.null(tag)) paste0("_", tag))
  #tag_me <- function(param) ifelse(is.na(param),"", paste(deparse(substitute(param)), round(param, digits=3), sep = "-"))
  #tag = paste(tag, tag_me(fst_host_strat), tag_me(fst_host_bias), tag_me(fst_pathogen_strat), tag_me(fst_pathogen_bias), "beta",scenario$beta, nrow(study_design), sep = "_")
  #ADD WARNING: When f coeff is defined but no stratification is defined
  tag_cm = paste(tag, "fps", fst_host_strat, "fpb", fst_host_bias, "fss", fst_pathogen_strat, "fsb", fst_pathogen_bias, "beta", scenario$beta, "s_size", nrow(study_design), sep = "_")
  
  print(paste0(Sys.time(), " : Scenario ",tag_cm))
  SNP = get_SNP(study_design, 
                tibble(`size` = scenario$rep, `Stratified` = scenario$s_stratified, `Partial_Stratification` = scenario$s_partial_stratification, `fst_strat` = fst_host_strat, `Biased`= scenario$s_biased, `Partial_Bias`= scenario$s_partial_bias,`fst_bias` = fst_host_bias))
  AA = apply(SNP$data, 2, function(snp.data) {
    get_AA(study_design,
           tibble(`size` = 1, `Stratified` = scenario$a_stratified, `Partial_Stratification` = scenario$a_partial_stratification, `fst_strat` = fst_pathogen_strat, `Biased`= scenario$a_biased, `Partial_Bias`= scenario$a_partial_bias, `fst_bias` = fst_pathogen_bias,
                  `Associated_pathogens`=scenario$associated_pathogens, `Associated_hosts`= scenario$associated_hosts, beta=scenario$beta), 
           associated_SNPs = as.matrix(snp.data))})
  AA = list(`data` = do.call(cbind, lapply(AA, function(aa) aa$data)), `freq` = do.call(cbind, lapply(AA, function(aa) aa$freq)))
  
  threshold <<- 0.05/((ncol(SNP$data)*ncol(AA$data)))
  host_pc = if(!is.na(sup_SNP_for_PC)) prcomp(cbind(SNP$data, do.call(rbind, sup_SNP_for_PC)), scale. = FALSE) else NULL
  res = analyse_G2G_setup(SNP$data, AA$data, study_design, host_pc)
  df = as.data.frame(t(do.call(rbind, lapply(res$pvalues, function(cor) cor$pval))))
  res = list(`study_design` = study_design, `pvalues` = res$pvalues, `scenario` = scenario, `pvalues_short` = df, `host_pc` = host_pc, `AA` = AA, `SNP` = SNP) 
  write_G2G_setup(res, tag_cm, paste0(getwd(),"/gen-data/",tag,"/" ))
  
  invisible(res)}

get_G2G_setup <- function(rep, s_stratified = NA, s_partial_strat = NA, s_biased = NA, s_partial_bias = NA, a_stratified = NA, a_partial_strat = NA, a_biased = NA, a_partial_bias = NA, associated_pathogens = NA, associated_hosts = NA, beta=NA) {
  associated = !is.na(associated_pathogens) | !is.na(associated_hosts)
  tbl_df(tibble(s_stratified = list(s_stratified), s_biased = list(s_biased), 
                s_partial_stratification = list(s_partial_strat), s_partial_bias = list(s_partial_bias),
                a_stratified = list(a_stratified), a_biased = list(a_biased),
                a_partial_stratification = list(a_partial_strat), a_partial_bias = list(a_partial_bias),
                associated_pathogens = list(associated_pathogens), associated_hosts = list(associated_hosts), beta, rep))}


### G2G full scenario methods
#TODO: add a "s" to method parse_G2G_config, it allow to make different G2G_config and replicate them.
#TODO: Remove use of environnement variable using recursive approach
#TODO: merge get_G2G_data with analyze_G2G
#TODO: beta in association()

#Here receive last AA_id, last SNP_id, turn replicate to for loops to set last id after last AA or SNP call.
#@`...`: AA() | SNP() | association()
get_G2G_conf <- function(..., bio_tag=NA, replicate = 1) {
  calls = match.call(expand.dots = FALSE)$`...`
  
  #Should be cleaned, first time using functionnal patern
  res = lapply(1:replicate, function(bio_tag_suffix) {
    lapply(calls, function(call){
      if(call[1] == "SNP()" || call[1] == "AA()") {
        call$bio_tag = if(!is.na(bio_tag)) paste(bio_tag, bio_tag_suffix, sep = "_") else paste(call$bio_tag, bio_tag_suffix,sep = "_")
        eval(call)
      }
      else if(call[1] == "association()"){
        call_mod = lapply(call[2:length(call)], function(ass_call){
          if(ass_call[1] == "SNP()" || ass_call[1] == "AA()") {
            ass_call$bio_tag = if(!is.na(bio_tag)) paste(bio_tag, bio_tag_suffix, sep = "_") else paste(ass_call$bio_tag, bio_tag_suffix,sep = "_")
            ass_call
          }
          else
            ass_call
        })
        do.call(association, call_mod)}
      else {
        stop(paste0(call[1], " is not a valid function call"))
      }
    })
  })
  ##Could merge scenario having same id_tag and different biotag, need to recalculate 'Size', and merge 'id'
  AA.scenarios = bind_rows(lapply(unlist(res,recursive =F), function(res) res$AA.scenarios))
  SNP.scenarios = bind_rows(lapply(unlist(res,recursive =F), function(res) res$SNP.scenarios))
  
  res = list(`AA.scenarios` = AA.scenarios, `SNP.scenarios` = SNP.scenarios)  }

SNP <- function(size, stratified = NA, partial_strat = NA, fst_strat=NA, biased = NA, partial_bias = NA, fst_bias=NA, bio_tag=NA) {
  list(`SNP.scenarios` = do.call(rbind, lapply(fst_strat, function(fst_strat) {
    do.call(rbind, lapply(fst_bias, function(fst_bias) {
      id =  get_id("SNP", size)
      id_tag = generate_id_tag(fst_strat, partial_strat, fst_bias, partial_bias)
      #      bio_tag = if(is.character(bio_tag)) setNames(list(id), bio_tag)  else generate_biological_tag(size, bio_tag, id)
      tibble(Stratified = list(stratified), Biased = list(biased), Partial_Stratification = list(partial_strat), Partial_Bias = list(partial_bias), `fst_strat` = fst_strat, `fst_bias`= fst_bias, size, id = list(id), id_tag, bio_tag)}))})))}

AA <- function(size, stratified = NA, partial_strat = NA, fst_strat=NA, biased = NA, partial_bias = NA, fst_bias=NA, 
               associated_pathogens = NA, associated_hosts = NA, beta=NA, bio_tag = NA) {
  
  if(!is.na(beta)) {
    associated_pathogens = "full"
    associated_hosts = "full"
  }
  
  list(`AA.scenarios` = do.call(rbind, lapply(fst_strat, function(fst_strat) {
    do.call(rbind, lapply(fst_bias, function(fst_bias) {
      do.call(rbind, lapply(beta, function(beta) {
        id =  get_id("AA", size)
        id_tag = generate_id_tag(fst_strat, partial_strat, fst_bias, partial_bias, beta)
        #        bio_tag = if(is.character(bio_tag)) setNames(list(id), bio_tag) else generate_biological_tag(size, bio_tag,id)
        tibble(Stratified = list(stratified), Biased = list(biased), 
               Partial_Stratification = list(partial_strat), Partial_Bias = list(partial_bias),
               Associated_pathogens = list(associated_pathogens), Associated_hosts = list(associated_hosts), 
               beta, `fst_strat` = fst_strat, `fst_bias` = fst_bias, size, id = list(id), id_tag, bio_tag)}))}))})))}

#@`...`: AA() | SNP()
association <- function(..., replicate = 1) {
  saved_eval = substitute(list(...))
  res = replicate(replicate, {
    res = eval(saved_eval)
    AA.scenarios = bind_rows(lapply(res, function(ele) ele$AA.scenarios))
    SNP.scenarios = bind_rows(lapply(res, function(ele) ele$SNP.scenarios))
    #Here make a thirs object, association table instead of putting it into AA.scenario
    AA.scenarios$associated_SNPs = list(unlist(SNP.scenarios$id))
    AA.scenarios$associated_SNP_tag = list(unique(unlist(SNP.scenarios$bio_tag)))
    list(`AA.scenarios` = AA.scenarios, `SNP.scenarios` = SNP.scenarios)}, simplify = F)
  list(`AA.scenarios` = bind_rows(lapply(res, function(ele) ele$AA.scenarios)), `SNP.scenarios` = bind_rows(lapply(res, function(ele) ele$SNP.scenarios)))}

#@`...`: G2G_conf()
get_G2G_data <- function(study_design, ...) {
  set_env("AA", 0)
  set_env("SNP", 0)
  #Turn it to calls list to set last AA_id and last SNP_id after last G2G_conf call
  s = list(...)
  AA.scenarios = bind_rows(lapply(s, function(ele) ele$AA.scenarios))
  SNP.scenarios = bind_rows(lapply(s, function(ele) ele$SNP.scenarios))
  SNP = get_SNP(study_design, SNP.scenarios)
  AA = get_AA(study_design, AA.scenarios, SNP$data)
  list(`AA.scenarios` = AA.scenarios,`AA.data` = AA$data,`AA.freq` = AA$freq, `SNP.scenarios` = SNP.scenarios, `SNP.data` = SNP$data, `SNP.freq` = SNP$freq, `study_design` = study_design)}

get_correction <-function(WO_correction = F, W_host_group = F, W_pathogen_group = F, W_both_groups = F, W_host_PC = F, W_pathogen_PC = F, W_both_PC = F, W_pathogen_groups_host_PC = F, W_non_linear_PC =F) {
  list(`WO_correction` = WO_correction, `W_host_group` = W_host_group, `W_pathogen_group` = W_pathogen_group, `W_both_groups` = W_both_groups, `W_host_PC` = W_host_PC, `W_pathogen_PC` = W_pathogen_PC, `W_both_PC` = W_both_PC, `W_pathogen_groups_host_PC` = W_pathogen_groups_host_PC, `W_non_linear_PC` = W_non_linear_PC)}

get_analyse <- function(logistic = F, skat_L  = F, skato_L = F, skat_LW  = F, skato_LW = F, gt = F, G2_analysis = F) {
  list(`logistic` = logistic, `skat_L` = skat_L, `skato_L` = skato_L, `skat_LW` = skat_LW, `skato_LW` = skato_LW, `gt` = gt, `G2_analysis` = G2_analysis)}

analyse_G2G <- function(G2G_data, correction, analyse = NA, nb_cpu = 1) {
  ##G2G_data contain both AA.data, SNP.data and study_design
  if(is.na(analyse)){
    analyse = get_analyse(logistic = T)
  }
  
  attach(G2G_data)
  attach(correction)
  attach(analyse)
  
  if(trace) print(paste0(Sys.time()," : Computing PC"))
  SNP_PC = if(W_host_PC || W_both_PC) {
    SNP_PC = prcomp(SNP.data)
    SNP_PC$x[,1:ifelse(ncol(SNP_PC$x)<5, ncol(SNP_PC$x), 5)]}
  
  AA_PC = if(W_pathogen_PC || W_both_PC){
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
      restest = c(if(WO_correction) list(`Without correction` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X)))[,4][2])),
                  if(W_host_group) list(`With host group` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+study_design[,"Host"], family = binomial)))[,4][2])),
                  if(W_pathogen_group) list(`With pathogen group` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+study_design[,"Pathogen"], family = binomial)))[,4][2])),
                  if(W_both_groups) list(`With both groups` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+study_design[,"Host"]+study_design[,"Pathogen"], family = binomial)))[,4][2])),
                  if(W_host_PC) list(`With host PC` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+SNP_PC, family = binomial)))[,4][2])),
                  if(W_pathogen_PC) list(`With pathogen PC` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+AA_PC, family = binomial)))[,4][2])),
                  if(W_both_PC) list(`With both PC` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+SNP_PC+AA_PC, family = binomial)))[,4][2])),
                  if(W_pathogen_groups_host_PC) list(`W_pathogen_groups_host_PC` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+SNP_PC+study_design[,"Pathogen"], family = binomial)))[,4][2])),
                  if(W_non_linear_PC) list(`With non linear PC` = apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+AA_NL_PC+study_design[,"Host"], family = binomial)))[,4][2])))
      restest
    }
    
    counter <<- 0
    nb_aa <<- ncol(AA.data)
    ptm <<- proc.time()
    
    cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
    res = parApply(cl, AA.data, 2, analyse_AA)
    stopCluster(cl)
    
    #res = apply(AA.data, 2, analyse_AA)
    
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
      if(W_host_group)  list(`With host group` = skat_it(SKAT_Null_Model(Y ~ study_design[,"Host"], out_type = "D"))),
      if(W_pathogen_group) list(`With pathogen group` = skat_it(SKAT_Null_Model(Y ~ study_design[,"Pathogen"], out_type = "D"))),
      if(W_both_groups) list(`With both groups` = skat_it(SKAT_Null_Model(Y ~ study_design[,"Host"]+study_design[,"Pathogen"], out_type = "D"))),
      if(W_host_PC) list(`With host PC` = skat_it(SKAT_Null_Model(Y ~ SNP_PC, out_type = "D"))),
      if(W_pathogen_PC) list(`With pathogen PC` = skat_it(SKAT_Null_Model(Y ~ AA_PC, out_type = "D"))),
      if(W_both_PC) list(`With both PC` = skat_it(SKAT_Null_Model(Y ~ AA_PC + SNP_PC, out_type = "D"))),
      if(W_non_linear_PC) list(`With non linear PC` = skat_it(SKAT_Null_Model(Y ~ AA_NL_PC, out_type = "D"))))}
    
    flip_SNP <- function(SNP){
      nsample = nrow(SNP)
      apply(SNP,2, function(S){
        if(sum(S)/nsample > 1) 
          unlist(lapply(S, function(s){
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
      if(W_host_group) list(`With host group` = lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ study_design[,"Host"], SNP.data[,batch]))["p-value"]))),
      if(W_pathogen_group) list(`With pathogen group` = lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ study_design[,"Pathogen"], SNP.data[,batch]))["p-value"]))),
      if(W_both_groups) list(`With both groups` = lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ study_design[,"Host"]+study_design[,"Pathogen"], SNP.data[,batch]))["p-value"]))),
      if(W_host_PC) list(`With host PC` = lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ SNP_PC, SNP.data[,batch]))["p-value"]))),
      if(W_pathogen_PC) list(`With pathogen PC` = lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ AA_PC, SNP.data[,batch]))["p-value"]))),
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
  
  detach(G2G_data)
  detach(correction)
  detach(analyse)
  
  list(`data` = data, `results` = res, `PC` = SNP_PC)
}
