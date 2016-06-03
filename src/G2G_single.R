source("G2G.R")

test_G2G_setup <- function(study_design, scenario, fst_pop_strat=NA, fst_pop_bias=NA, fst_strain_strat=NA, fst_strain_bias=NA, tag = NULL, get_viral = generate_AAs_standard) {
	scenario_name = paste0(deparse(substitute(scenario)), if(!is.null(tag)) paste0("_", tag))
	tag_me <- function(param) ifelse(is.na(param),"", paste(deparse(substitute(param)), round(param, digits=2), sep = "-"))
	tag = paste(scenario_name, tag_me(fst_pop_strat), tag_me(fst_pop_bias), tag_me(fst_strain_strat), tag_me(fst_strain_bias), "beta",scenario$beta, nrow(study_design), sep = "_")
	print(paste0(Sys.time(), " : Scenario ",tag))
	SNP = get_SNP(study_design, data_frame(`size` = scenario$rep, `Stratified` = scenario$S_Stratified, `Partial_Stratification` = scenario$S_Partial_Stratification, `fst_strat` = fst_pop_strat, `Biased`= scenario$S_Biased, `Partial_Bias`= scenario$S_Partial_Bias,`fst_bias` = fst_pop_bias))
	AA = apply(SNP, 2, function(snp.data) {
		get_AA(study_design, data_frame(`size` = 1, `Stratified` = scenario$A_Stratified, `Partial_Stratification` = scenario$Y_Partial_Stratification, `fst_strat` = fst_strain_strat, `Biased`= scenario$A_Biased, `Partial_Bias`= scenario$Y_Partial_Bias, `fst_bias` = fst_strain_bias,
																		`Associated_Strains`=scenario$Associated_Strains, `Associated_Populations`= scenario$Associated_Populations, beta=scenario$beta), associated_SNPs = as.matrix(snp.data), generate_AAs = get_viral)})
	
	threshold <<- 0.05/((ncol(SNP)*ncol(AA)))
	pvalues = analyse_G2G_setup(SNP, AA, study_design)
	
	###TODO fix student
	df = as.data.frame(t(do.call(rbind, lapply(pvalues, function(cor) cor$pval))))
	student = if(ncol(df) == 4) list(
		`Without correction VS With human groups` = t.test(x = df[,"WO_correction"], y = df[,"W_human_groups"])$p.value,
		`Without correction VS With strains groups` = t.test(x = df[,"WO_correction"], y = df[,"W_strain_groups"])$p.value,
		`Without correction VS With both groups` = t.test(x = df[,"WO_correction"], y = df[,"W_both_groups"])$p.value,
		`With both groups VS With strains groups` = t.test(x = df[,"W_both_groups"], y = df[,"W_strain_groups"])$p.value,
		`With both groups VS With human groups` = t.test(x = df[,"W_both_groups"], y = df[,"W_human_groups"])$p.value)
	student = data.frame("pvalues" = unlist(student))
	
	res = list(`study_design` = study_design, `pvalues` = pvalues, `scenario` = scenario, `student`= student, `pvalues_short` = df) #, `params` = cbind(SNP$SNP.params, AA$AA.params)
	write_G2G_setup(res, tag, paste0(getwd(),"/../gen-data/",scenario_name,"/" ))
	invisible(res)}

get_G2G_setup <- function(rep, s_stratified = NA, s_partial_strat = NA, s_biased = NA, s_partial_bias = NA, a_stratified = NA, y_partial_strat = NA, a_biased = NA, y_partial_bias = NA, associated_strains = NA, associated_populations = NA, beta=NA) {
	associated = !is.na(associated_strains) | !is.na(associated_populations)
	tbl_df(data_frame(S_Stratified = list(s_stratified), S_Biased = list(s_biased), 
										S_Partial_Stratification = list(s_partial_strat), S_Partial_Bias = list(s_partial_bias),
										A_Stratified = list(a_stratified), A_Biased = list(a_biased),
										Y_Partial_Stratification = list(y_partial_strat), Y_Partial_Bias = list(y_partial_bias),
										Associated_Strains = list(associated_strains), Associated_Populations = list(associated_populations), beta, rep) )}

analyse_G2G_setup <- function(SNP, AA, study_design) {
	if(trace == TRUE) print(paste(Sys.time(),"Computing GLM", sep=" : "))
	pvalues = list(`WO_correction` = unlist(lapply(1:ncol(AA), function(index) coef(summary(glm(AA[,index]~SNP[,index])))[,4][2])),
								 `W_human_groups` = if(length(levels(study_design[,"Population"])) > 1) unlist(lapply(1:ncol(AA), function(index) coef(summary(glm(AA[,index]~SNP[,index]+study_design[,"Population"])))[,4][2])),
								 `W_strain_groups` = if(length(levels(study_design[,"Strain"])) > 1) unlist(lapply(1:ncol(AA), function(index) coef(summary(glm(AA[,index]~SNP[,index]+study_design[,"Strain"])))[,4][2])),
								 `W_both_groups` = if(length(levels(study_design[,"Population"])) > 1 && length(levels(study_design[,"Strain"])) > 1) unlist(lapply(1:ncol(AA), function(index) coef(summary(glm(AA[,index]~SNP[,index]+study_design[,"Population"]+study_design[,"Strain"])))[,4][2])))
	if(trace == TRUE) print(paste(Sys.time(),"Parsing values", sep=" : "))   
	parse_pvalues(as.data.frame(Filter(length, pvalues)), threshold)}

###@G2SR : G2G setup results
plot_G2G_setup <- function(G2SR, save = FALSE, out = TRUE, lim = 15, file = Sys.time()) {
	if(trace == TRUE) ifelse((save == FALSE), print(paste(Sys.time(),"Ploting...", sep=" : ")), print(paste(Sys.time(),"Writting plots", sep=" : "))) 
	G2GSR = -log10(select(as.data.frame(G2SR$pvalues), ends_with("pval")))
	G2GSRD = select(as.data.frame(G2SR$pvalues), ends_with("pval_diff"))
	t=mapply(function(res, name) {
		colnames(res) = gsub("WO_correction.pval", "Without correction", colnames(res))
		colnames(res) = gsub("W_human_groups.pval", "With human groups", colnames(res))
		colnames(res) = gsub("W_strain_groups.pval", "With strain groups", colnames(res))
		colnames(res) = gsub("W_both_groups.pval", "With human and strain groups", colnames(res))
		res = res %>% gather(factor_key = T)
		p <- ggplot(res, aes(key, value))
		p <- if(out == TRUE) {
			# compute lower and upper whiskers
			ylim1 = boxplot.stats(res$value)$stats[c(1, 5)]
			# scale y limits based on ylim1
			p + coord_cartesian(ylim = ylim1*1.05)
		} else p
		p + geom_boxplot(outlier.color = "grey", aes(fill=key)) + #geom_hline(yintercept = -log10(threshold), colour="red") + 
			labs(title = paste(name, " with different covariates"), y = "-log10(pval)", x="covariate") + 
			theme(axis.text.x=element_blank())
		if(save == TRUE) ggsave(filename = paste0(file,name,".png"), width = 14, height = 18)
		p <- ggplot(res, aes(key, value))
		p + geom_boxplot(outlier.color = "grey", aes(fill=key)) +
			labs(title = paste(name, " with different covariates"), y = "-log10(pval)", x="covariate") + 
			coord_cartesian(ylim = c(0, lim)) +
			theme(axis.text.x=element_blank())
		if(save == TRUE) ggsave(filename = paste0(file,name,"lim-",lim,".png"), width = 14, height = 18)
	}, list(G2GSR, G2GSRD), c("pvalues", "pvalues relative difference from corection"))}

write_G2G_setup <- function(res, tag = NULL, output_dir = paste0(getwd(), "/")) {
	#end = paste0("-(", Sys.time(),")")
	end = ""
	dir.create(paste(output_dir), showWarnings = TRUE, recursive = FALSE, mode = "0777")
	plot_G2G_setup(res, out = TRUE, save = TRUE, file = paste0(output_dir, tag, "-boxplots_out", end))
	plot_G2G_setup(res, out = FALSE, save = TRUE, file = paste0(output_dir, tag, "-boxplots", end))
	write.table(x = res$student, file = paste0(output_dir, tag, "-students-(", end, ".csv"))
	write.table(x = res$pvalues_short, file = paste0(output_dir, tag, "-pvalues", end, ".csv"))}