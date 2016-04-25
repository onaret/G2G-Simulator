library(reshape2)

##@G
to_pop_structure <- function(study_design) {
	do.call(rbind, lapply(levels(study_design$Population), function(population) do.call(rbind, lapply(levels(study_design$Strain), function(strain) data_frame(strain,population, `nb` = sum(study_design[,"Strain"] == strain & study_design[,"Population"] == population))))))}

##@G
to_study_design <- function(structure) {
	do.call(rbind, lapply(1:nrow(structure), function(pop_num) {
		do.call(rbind, lapply(1:ncol(structure), function(strain_num) {
			data.frame(Strain = rep(chartr("123456789", "ABCDEFGHI", strain_num), structure[strain_num,pop_num]), Population = paste0("P",pop_num) ) }))})) }

##@G
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

test_G2G_setup <- function(study_design, scenario, fst_pop_strat=NA, fst_pop_bias=NA, fst_strain_strat=NA, fst_strain_bias=NA, tag = NULL, get_viral = generate_AAs) {
	scenario_name = paste0(deparse(substitute(scenario)), if(!is.null(tag)) paste0("_", tag))
	tag_me <- function(param) ifelse(is.na(param),"", paste(deparse(substitute(param)), round(param, digits=2), sep = "-"))
	tag = paste(scenario_name, tag_me(fst_pop_strat), tag_me(fst_pop_bias), tag_me(fst_strain_strat), tag_me(fst_strain_bias), "beta",scenario$beta, nrow(study_design), sep = "_")
	print(paste0(Sys.time(), " : Scenario ",tag))
	SNP = get_SNP(study_design, data_frame(`size` = scenario$rep, `Stratified` = scenario$S_Stratified, `Partial_Stratification` = scenario$S_Partial_Stratification, `fst_strat` = fst_pop_strat, `Biased`= scenario$S_Biased, `Partial_Bias`= scenario$S_Partial_Bias,`fst_bias` = fst_pop_bias))
	#AA = get_AA(study_design, data_frame(`size` = scenario$rep, `Stratified` = scenario$Y_Stratified, `Partial_Stratification` = scenario$Y_Partial_Stratification, `fst_strat` = fst_strain_strat, `Biased`= scenario$Y_Biased, `Partial_Bias`= scenario$Y_Partial_Bias, `fst_bias` = fst_strain_bias,
	#						`Associated_Strains`=scenario$Associated_Strains, `Associated_Populations`= scenario$Associated_Populations, beta=scenario$beta), SNP)
	
	AA = apply(SNP, 2, function(snp.data) {
		get_AA(study_design, data_frame(`size` = 1, `Stratified` = scenario$Y_Stratified, `Partial_Stratification` = scenario$Y_Partial_Stratification, `fst_strat` = fst_strain_strat, `Biased`= scenario$Y_Biased, `Partial_Bias`= scenario$Y_Partial_Bias, `fst_bias` = fst_strain_bias,
		`Associated_Strains`=scenario$Associated_Strains, `Associated_Populations`= scenario$Associated_Populations, beta=scenario$beta), associated_SNPs = as.matrix(snp.data))
	})
	
	threshold <<- 0.05/(ncol(SNP+ncol(AA)))
	pvalues = analyse_G2G_setup(SNP, AA, study_design)
	#summary_sim = summary_sim(pvalues, AA$AA.scenario)
	res = list(`study_design` = study_design, `pvalues` = pvalues, `scenario` = scenario) #, `params` = cbind(SNP$SNP.params, AA$AA.params)
	write_G2G_setup(res, tag, paste0(getwd(),"/../gen-data/",scenario_name,"/" ))
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

###@G2SR : G2G setup results
plot_G2G_setup <- function(G2SR, save = FALSE, out = TRUE, file = Sys.time()) {
	if(trace == TRUE) ifelse((save == FALSE), print(paste(Sys.time(),"Ploting...", sep=" : ")), print(paste(Sys.time(),"Writting plots", sep=" : "))) 
	G2GSR = -log10(select(as.data.frame(G2SR$pvalues), ends_with("pval")))
	G2GSRD = select(as.data.frame(G2SR$pvalues), ends_with("pval_diff"))
	t=mapply(function(res, name) {
		colnames(res) <- c("Without correction", "With human groups", "With viral groups", "With human and viral groups")
		res = melt(res, measure.vars = 1:4)
		p <- ggplot(res, aes(variable, value))
		p <- if(out == TRUE) {
			# compute lower and upper whiskers
			ylim1 = boxplot.stats(res$value)$stats[c(1, 5)]
			# scale y limits based on ylim1
			p + coord_cartesian(ylim = ylim1*1.05)
		} else p
		p + geom_boxplot(outlier.color = "grey", aes(fill=variable)) + geom_hline(yintercept = -log10(threshold), colour="red") + 
			labs(title = paste(name, "pvalues with different covariates"), y = "-log10(pval)", x="covariate") + 
			theme(axis.text.x=element_blank(), axis.text.y = element_text(size=24), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36))# + 
			#if(!is.null(lim)) scale_y_continuous(limits = lim) +
			#guides(fill=FALSE)
		if(save == TRUE) ggsave(filename = paste0(file,name,".png"), width = 14, height = 18)
	}, list(G2GSR, G2GSRD), c("G2GR", "G2GRD"))}

write_G2G_setup <- function(res, tag = NULL, output_dir = paste0(getwd(), "/")) {
	end = Sys.time()
	dir.create(paste(output_dir), showWarnings = TRUE, recursive = FALSE, mode = "0777")
	plot_GWAS(res$pvalues, save = TRUE, file = paste0(output_dir, "-plots-(", end, ")"))
	plot_G2G_setup(res, out = TRUE, save = TRUE, file = paste0(output_dir, tag, "-boxplots_out-(", end, ")"))
	plot_G2G_setup(res, out = FALSE, save = TRUE, file = paste0(output_dir, tag, "-boxplots-(", end, ")"))
	#res$SNP_params[,"Associated_Populations"] = vapply(res$SNP_params[,"Associated_Populations"], paste, collapse = ", ", character(1L))
	#res$SNP_params[,"Associated_Strains"] = vapply(res$SNP_params[,"Associated_Strains"], paste, collapse = ", ", character(1L))
	#write.table(x = cbind(res$SNP_params,res$pvalues), file = paste0(output_dir, tag, "-SNP_params-pvals-(", end, ").csv"))
	#write.table(x = res$study_design, file = paste0(output_dir, tag, "-study_design-(", end, ").csv"))
	#write.table(x = res$params, file = paste0(output_dir, tag, "-params-(", end, ").csv"))
	#write.table(x = res$scenario, file = paste0(output_dir, tag, "-scenario-(", end, ").csv"))
	#write.table(x =res$summary_sim, file = paste0(output_dir, tag, "-summary-(", end, ").csv"))
	}