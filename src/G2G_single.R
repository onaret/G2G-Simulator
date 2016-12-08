source("G2G.R")

test_G2G_setup <- function(study_design, scenario, fst_pop_strat=NA, fst_pop_bias=NA, fst_strain_strat=NA, fst_strain_bias=NA, tag = "unnamed", get_viral = generate_AAs_standard, sup_SNP_for_PC = NULL) {
	#scenario_name = paste0(deparse(substitute(scenario)), if(!is.null(tag)) paste0("_", tag))
	#tag_me <- function(param) ifelse(is.na(param),"", paste(deparse(substitute(param)), round(param, digits=3), sep = "-"))
	#tag = paste(tag, tag_me(fst_pop_strat), tag_me(fst_pop_bias), tag_me(fst_strain_strat), tag_me(fst_strain_bias), "beta",scenario$beta, nrow(study_design), sep = "_")
  #ADD WARNING: When f coeff is defined but no stratification is defined
  tag_cm = paste(tag, "fps", fst_pop_strat, "fpb", fst_pop_bias, "fss", fst_strain_strat, "fsb", fst_strain_bias, "beta", scenario$beta, "s_size", nrow(study_design), sep = "_")
  	
  print(paste0(Sys.time(), " : Scenario ",tag_cm))
	SNP = get_SNP(study_design, 
	      data_frame(`size` = scenario$rep, `Stratified` = scenario$s_stratified, `Partial_Stratification` = scenario$s_partial_stratification, `fst_strat` = fst_pop_strat, `Biased`= scenario$s_biased, `Partial_Bias`= scenario$s_partial_bias,`fst_bias` = fst_pop_bias))
	AA = apply(SNP, 2, function(snp.data) {
		get_AA(study_design,
		       data_frame(`size` = 1, `Stratified` = scenario$a_stratified, `Partial_Stratification` = scenario$a_partial_stratification, `fst_strat` = fst_strain_strat, `Biased`= scenario$a_biased, `Partial_Bias`= scenario$a_partial_bias, `fst_bias` = fst_strain_bias,
										`Associated_Strains`=scenario$associated_strains, `Associated_Populations`= scenario$associated_populations, beta=scenario$beta), 
		       associated_SNPs = as.matrix(snp.data), generate_AAs = get_viral)})
	
	threshold <<- 0.05/((ncol(SNP)*ncol(AA)))
	
	host_pc = if(!is.null(sup_SNP_for_PC)) prcomp(cbind(SNP, sup_SNP_for_PC), scale. = FALSE) else NULL
	res = analyse_G2G_setup(SNP, AA, study_design, host_pc)
	df = as.data.frame(t(do.call(rbind, lapply(res$pvalues, function(cor) cor$pval))))
	res = list(`study_design` = study_design, `pvalues` = res$pvalues, `scenario` = scenario, `pvalues_short` = df) 
	write_G2G_setup(res, tag_cm, paste0(getwd(),"/../gen-data/",tag,"/" ))
	
	invisible(res)}

get_G2G_setup <- function(rep, s_stratified = NA, s_partial_strat = NA, s_biased = NA, s_partial_bias = NA, a_stratified = NA, a_partial_strat = NA, a_biased = NA, a_partial_bias = NA, associated_strains = NA, associated_populations = NA, beta=NA) {
	associated = !is.na(associated_strains) | !is.na(associated_populations)
	tbl_df(data_frame(s_stratified = list(s_stratified), s_biased = list(s_biased), 
										s_partial_stratification = list(s_partial_strat), s_partial_bias = list(s_partial_bias),
										a_stratified = list(a_stratified), a_biased = list(a_biased),
										a_partial_stratification = list(a_partial_strat), a_partial_bias = list(a_partial_bias),
										associated_strains = list(associated_strains), associated_populations = list(associated_populations), beta, rep) )}

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

write_G2G_setup <- function(res, tag = NULL, output_dir = paste0(getwd(), "/")) {
  dir.create(paste(output_dir), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  #plot_G2G_setup(res, out = TRUE, save = TRUE, file = paste0(output_dir, tag, "-boxplots_out"))
  plot_G2G_setup(res, out = FALSE, save = TRUE, title = tag, file = paste0(output_dir, tag, "-boxplots"))
  #write.csv(x = summary(res$pvalues_shor),file = paste0(output_dir, tag, "-summary", end, ".csv"))
  write.table(x = res$pvalues_short, file = paste0(output_dir, tag, "-pvalues", ".csv"))}

###@G2SR : G2G setup results
plot_G2G_setup <- function(G2SR, save = FALSE, out = TRUE, lim = 15, title ="", file = Sys.time()) {
	if(trace == TRUE) ifelse((save == FALSE), print(paste(Sys.time(),"Ploting...", sep=" : ")), print(paste(Sys.time(),"Writting plots", sep=" : "))) 
  G2GSR = -log10(select(as.data.frame(G2SR$pvalues), ends_with("pval")))
	G2GSRD = select(as.data.frame(G2SR$pvalues), ends_with("pval_diff"))
	t=mapply(function(res, name) {
		colnames(res) = gsub("WO_correction.pval", "Without correction", colnames(res))
		colnames(res) = gsub("W_human_groups.pval", "With human groups", colnames(res))
		colnames(res) = gsub("W_human_PCs.pval", "With human PCs", colnames(res))
		colnames(res) = gsub("W_strain_groups.pval", "With strain groups", colnames(res))
		colnames(res) = gsub("W_both.pval", "With both human and strain covariates", colnames(res))
		res = res %>% gather(factor_key = T)
		p <- ggplot(res, aes(key, value))
		p <- if(out == TRUE) {
			# compute lower and upper whiskers
			ylim1 = boxplot.stats(res$value)$stats[c(1, 5)]
			# scale y limits based on ylim1
			p + coord_cartesian(ylim = ylim1*1.05)
		} else p
		p + geom_boxplot(outlier.color = "grey", aes(fill=key)) + #geom_hline(yintercept = -log10(threshold), colour="red") + 
			labs(title = paste(name, title," with different covariates"), y = "-log10(pval)", x="covariate") + 
			theme(axis.text.x=element_blank())
		#if(save == TRUE) ggsave(filename = paste0(file,name,".png"), width = 14, height = 18)
		p <- ggplot(res, aes(key, value))
		p + geom_boxplot(outlier.color = "grey", aes(fill=key)) +
			labs(title = paste(name, title, " with different covariates"), y = "-log10(pval)", x="covariate") + 
			coord_cartesian(ylim = c(0, lim)) +
			theme(axis.text.x=element_blank())
		if(save == TRUE) ggsave(filename = paste0(file, name, "lim-",lim,".png"), width = 14, height = 18)
	}, list(G2GSR, G2GSRD), c("_res", "_difference_from_correction"))}
