###G2G common methods
to_pop_structure <- function(study_design) {
  do.call(rbind, lapply(levels(study_design$Population), function(population) {
    do.call(rbind, lapply(levels(study_design$Strain), function(strain) { 
      data_frame(strain,population, `nb` = sum(study_design[,"Strain"] == strain & study_design[,"Population"] == population))
    })) 
  }))
}

get_study_design <- function(structure) {
  structure = as.data.frame(structure)
  study_design = do.call(rbind, lapply(1:ncol(structure), function(pop_num) {
    do.call(rbind, lapply(1:nrow(structure), function(strain_num) {
      if(structure[strain_num,pop_num] > 0) {
        data.frame(Population = paste0("P",pop_num), Strain = rep(chartr("123456789", "ABCDEFGHI", strain_num), structure[strain_num,pop_num]) )
      }
    }))
  }))
  study_design$Population = as.factor(study_design$Population)
  study_design$Strain = as.factor(study_design$Strain)
  study_design
}


###G2G single setup methods
test_G2G_setup <- function(study_design, scenario, fst_pop_strat=NA, fst_pop_bias=NA, fst_strain_strat=NA, fst_strain_bias=NA, 
                           tag = "unnamed", sup_SNP_for_PC = NA) {
  #scenario_name = paste0(deparse(substitute(scenario)), if(!is.null(tag)) paste0("_", tag))
  #tag_me <- function(param) ifelse(is.na(param),"", paste(deparse(substitute(param)), round(param, digits=3), sep = "-"))
  #tag = paste(tag, tag_me(fst_pop_strat), tag_me(fst_pop_bias), tag_me(fst_strain_strat), tag_me(fst_strain_bias), "beta",scenario$beta, nrow(study_design), sep = "_")
  #ADD WARNING: When f coeff is defined but no stratification is defined
  tag_cm = paste(tag, "fps", fst_pop_strat, "fpb", fst_pop_bias, "fss", fst_strain_strat, "fsb", fst_strain_bias, "beta", scenario$beta, "s_size", nrow(study_design), sep = "_")
  
  print(paste0(Sys.time(), " : Scenario ",tag_cm))
  SNP = get_SNP(study_design, 
                data_frame(`size` = scenario$rep, `Stratified` = scenario$s_stratified, `Partial_Stratification` = scenario$s_partial_stratification, `fst_strat` = fst_pop_strat, `Biased`= scenario$s_biased, `Partial_Bias`= scenario$s_partial_bias,`fst_bias` = fst_pop_bias))
  AA = apply(SNP$data, 2, function(snp.data) {
    get_AA(study_design,
           data_frame(`size` = 1, `Stratified` = scenario$a_stratified, `Partial_Stratification` = scenario$a_partial_stratification, `fst_strat` = fst_strain_strat, `Biased`= scenario$a_biased, `Partial_Bias`= scenario$a_partial_bias, `fst_bias` = fst_strain_bias,
                      `Associated_Strains`=scenario$associated_strains, `Associated_Populations`= scenario$associated_populations, beta=scenario$beta), 
           associated_SNPs = as.matrix(snp.data))})
  AA = list(`data` = do.call(cbind, lapply(AA, function(aa) aa$data)), `freq` = do.call(cbind, lapply(AA, function(aa) aa$freq)))
  
  threshold <<- 0.05/((ncol(SNP$data)*ncol(AA$data)))
  host_pc = if(!is.na(sup_SNP_for_PC)) prcomp(cbind(SNP$data, do.call(rbind, sup_SNP_for_PC)), scale. = FALSE) else NULL
  res = analyse_G2G_setup(SNP$data, AA$data, study_design, host_pc)
  df = as.data.frame(t(do.call(rbind, lapply(res$pvalues, function(cor) cor$pval))))
  res = list(`study_design` = study_design, `pvalues` = res$pvalues, `scenario` = scenario, `pvalues_short` = df, `host_pc` = host_pc, `AA` = AA, `SNP` = SNP) 
  write_G2G_setup(res, tag_cm, paste0(getwd(),"/gen-data/",tag,"/" ))
  
  invisible(res)}

get_G2G_setup <- function(rep, s_stratified = NA, s_partial_strat = NA, s_biased = NA, s_partial_bias = NA, a_stratified = NA, a_partial_strat = NA, a_biased = NA, a_partial_bias = NA, associated_strains = NA, associated_populations = NA, beta=NA) {
  associated = !is.na(associated_strains) | !is.na(associated_populations)
  tbl_df(data_frame(s_stratified = list(s_stratified), s_biased = list(s_biased), 
                    s_partial_stratification = list(s_partial_strat), s_partial_bias = list(s_partial_bias),
                    a_stratified = list(a_stratified), a_biased = list(a_biased),
                    a_partial_stratification = list(a_partial_strat), a_partial_bias = list(a_partial_bias),
                    associated_strains = list(associated_strains), associated_populations = list(associated_populations), beta, rep))}


###G2G full scenario methods
#TODO: add a "s" to method parse_G2G_config, it allow to make different G2G_config and replicate them.
#TODO: Remove use of environnement variable using recursive approach
#TODO merge parse_G2G_conf with analyze_G2G
#TODO beta in association()

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

AA <- function(size, stratified = NA, partial_strat = NA, fst_strat=NA, biased = NA, partial_bias = NA, fst_bias=NA, 
               associated_strains = NA, associated_populations = NA, beta=NA, bio_tag = NA) {
  
  if(!is.na(beta)) {
    associated_strains = "full"
    associated_populations = "full"
  }
  
  list(`AA.scenarios` = do.call(rbind, lapply(fst_strat, function(fst_strat) {
    do.call(rbind, lapply(fst_bias, function(fst_bias) {
      do.call(rbind, lapply(beta, function(beta) {
        id =  get_id("AA", size)
        id_tag = generate_id_tag(fst_strat, partial_strat, fst_bias, partial_bias, beta)
        #        bio_tag = if(is.character(bio_tag)) setNames(list(id), bio_tag) else generate_biological_tag(size, bio_tag,id)
        data_frame(Stratified = list(stratified), Biased = list(biased), 
                   Partial_Stratification = list(partial_strat), Partial_Bias = list(partial_bias),
                   Associated_Strains = list(associated_strains), Associated_Populations = list(associated_populations), 
                   beta, `fst_strat` = fst_strat, `fst_bias` = fst_bias, size, id = list(id), id_tag, bio_tag)}))}))})))}

SNP <- function(size, stratified = NA, partial_strat = NA, fst_strat=NA, biased = NA, partial_bias = NA, fst_bias=NA, bio_tag=NA) {
  list(`SNP.scenarios` = do.call(rbind, lapply(fst_strat, function(fst_strat) {
    do.call(rbind, lapply(fst_bias, function(fst_bias) {
      id =  get_id("SNP", size)
      id_tag = generate_id_tag(fst_strat, partial_strat, fst_bias, partial_bias)
      #      bio_tag = if(is.character(bio_tag)) setNames(list(id), bio_tag)  else generate_biological_tag(size, bio_tag, id)
      data_frame(Stratified = list(stratified), Biased = list(biased), Partial_Stratification = list(partial_strat), Partial_Bias = list(partial_bias), `fst_strat` = fst_strat, `fst_bias`= fst_bias, size, id = list(id), id_tag, bio_tag)}))})))}

correction <-function(WO_correction = F, W_human_group = F, W_strain_group = F, W_both_groups = F, W_human_PC = F, W_strain_PC = F, W_both_PC = F, W_strain_groups_human_PC = F, W_non_linear_PC =F) {
  list(`WO_correction` = WO_correction, `W_human_group` = W_human_group, `W_strain_group` = W_strain_group, `W_both_groups` = W_both_groups, `W_human_PC` = W_human_PC, `W_strain_PC` = W_strain_PC, `W_both_PC` = W_both_PC, `W_strain_groups_human_PC` = W_strain_groups_human_PC, `W_non_linear_PC` = W_non_linear_PC)}

analyse <- function(logistic = F, skat_L  = F, skato_L = F, skat_LW  = F, skato_LW = F, gt = F, G2_analysis = F) {
  list(`logistic` = logistic, `skat_L` = skat_L, `skato_L` = skato_L, `skat_LW` = skat_LW, `skato_LW` = skato_LW, `gt` = gt, `G2_analysis` = G2_analysis)}
