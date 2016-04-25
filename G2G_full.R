parse_G2G_config <- function(study_design, ...) {
	set_env("AA", 0)
	set_env("SNP", 0)
	AA.scenarios = rbind_all(lapply(list(...), function(ele) ele$AA))
	SNP.scenarios = rbind_all(lapply(list(...), function(ele) ele$SNP))
	SNP.data = get_SNP(study_design, SNP.scenarios)
	AA.data = get_AA(study_design, AA.scenarios, SNP.data)
	list(`AA.scenarios` = AA.scenarios,`AA.data` = AA.data, `SNP.data` = SNP.data, `SNP.scenarios` = SNP.scenarios)}

G2G_conf<- function(...,bio_tag, replicate = 1) {
	calls = match.call(expand.dots = FALSE)$`...`
	res = lapply(paste0(bio_tag,"_", 1:replicate), function(bio_tag) {
		lapply(calls, function(call){
			if(call[1] == "SNP()" || call[1] == "AA()") {
				call$bio_tag = bio_tag
				eval(call)}
			else if(call[1] == "association()"){
				call_mod = lapply(call[2:length(call)], function(ass_call){
					ass_call$bio_tag = bio_tag
					ass_call})
				do.call(association, call_mod)}
			else {stop(paste0(call[1], " is not a valid function call"))}})})
	##Could merge scenario having same id_tag and different biotag, need to recalculate 'Size', and merge 'id'
	AA = rbind_all(lapply(unlist(res,recursive =F), function(res) res$AA))
	SNP = rbind_all(lapply(unlist(res,recursive =F), function(res) res$SNP))
	list(`AA` = AA, `SNP` = SNP)}

association <- function(...) {
	res = list(...)
	AA = rbind_all(lapply(res, function(ele) ele$AA))
	SNP = rbind_all(lapply(res, function(ele) ele$SNP))
	AA$associated_SNPs = SNP$id
	AA$associated_SNP_tag = names(unlist(SNP$bio_tag, recursive = F))
	list(`AA` = AA, `SNP` = SNP)	}

AA <- function(size, stratified = NA, partial_strat = NA, fst_strat=NA, biased = NA, partial_bias = NA, fst_bias=NA, associated_strains = NA, associated_populations =NA, beta=NA, bio_tag) {
	list(`AA` = do.call(rbind, lapply(fst_strat, function(fst_strat) {
		do.call(rbind, lapply(fst_bias, function(fst_bias) {
			do.call(rbind, lapply(beta, function(beta) {
				id =  get_id("AA",size)
				id_tag = generate_id_tag(fst_strat, partial_strat, fst_bias, partial_bias, beta)
				bio_tag = if(is.character(bio_tag)) setNames(list(id), bio_tag) else generate_biological_tag(size, bio_tag,id)
				data_frame(Stratified = list(stratified), Biased = list(biased), Partial_Stratification = list(partial_strat), Partial_Bias = list(partial_bias),
									 Associated_Strains= list(associated_strains), Associated_Populations = list(associated_populations), beta, `fst_strat` = fst_strat, `fst_bias` = fst_bias, size, id = list(id), id_tag, bio_tag = list(bio_tag))}))}))})))}

SNP <- function(size, stratified = NA, partial_strat = NA, fst_strat=NA, biased = NA, partial_bias = NA, fst_bias=NA, bio_tag) {
	list(`SNP` = do.call(rbind, lapply(fst_strat, function(fst_strat) {
		do.call(rbind, lapply(fst_bias, function(fst_bias) {
			id =  get_id("SNP",size)
			id_tag = generate_id_tag(fst_strat, partial_strat, fst_bias, partial_bias)
			bio_tag = if(is.character(bio_tag)) setNames(list(id), bio_tag)  else generate_biological_tag(size, bio_tag, id)
			data_frame(Stratified = list(stratified), Biased = list(biased), Partial_Stratification = list(partial_strat), Partial_Bias = list(partial_bias), `fst_strat` = fst_strat, `fst_bias`= fst_bias, size, id = list(id), id_tag, bio_tag = list(bio_tag))}))})))}

generate_id_tag <- function(fst_strat=NA, partial_strat=NA, fst_bias=NA, partial_bias=NA, beta=NA) {
	round_me <- function(prefix, param) ifelse(is.na(param),"", paste0(prefix,"-",round(param, digits=2)))
	part_me <- function(param) ifelse(is.na(param),"", paste0("_P-", paste0(param, collapse = ',')))
	tag = c(paste0(round_me("Fstrat", fst_strat), part_me(partial_strat)), paste0(round_me("Fbias", fst_bias), part_me(partial_bias)), round_me("Beta",beta))
	tag = paste0(Filter(function(x) x!="", tag), collapse = "_")
	tag = if(tag=="") "homogenous" else tag}

generate_biological_tag <- function(size, tag_size_range, id) {
	incrementer <-function() round(runif(1, min = tag_size_range["min"], max = tag_size_range["max"]))
	get_biotag_range <- function(sequence=c()) {
		if(size - sum(sequence) < tag_size_range["max"]) c(sequence, size - sum(sequence))
		else get_biotag_range(c(sequence, incrementer()))}
	bio_tag_range = get_biotag_range()
	get_bio_tag <- function (id, num=1) {
		res = setNames(list(id[1:bio_tag_range[num]]),	paste0(id[1],"_BT_",num,"_",bio_tag_range[num]))
		id = id[(bio_tag_range[num]+1):length(id)]
		if(length(bio_tag_range) >= num + 1) c(res, get_bio_tag(id, num+1)) else res}
	get_bio_tag(id)}

get_id <- function(prefix,size) {
	id = if(Sys.getenv(prefix) == "") 0 else as.numeric(Sys.getenv(prefix))
	args = list(id + size)
	names(args) = prefix
	do.call(Sys.setenv, args)
	paste0(prefix, "_",(id+1):(id+size))}

set_env <- function(name, value) {
	args = list(value)
	names(args) = name
	do.call(Sys.setenv, args)}

###associated_SNPs used for G2G setup
get_AA <- function(study_design, AA.scenarios, SNP.data, associated_SNPs = NULL) {
	populations = levels(study_design$Population)
	strains = levels(study_design$Strain)
	do.call(cbind, lapply(1:nrow(AA.scenarios), function(scenario_num) {
		AA.scenario = AA.scenarios[scenario_num,]
		associated_populations = unlist(AA.scenario$Associated_Populations)
		associated_strains = unlist(AA.scenario$Associated_Strains)
		AA.scenario$Associated_Populations = list(if(!is.na(associated_populations)) if(associated_populations == "full") populations else if(associated_populations == "half") sample(populations, 1) else associated_populations)
		AA.scenario$Associated_Strains = list(if(!is.na(associated_strains)) if(associated_strains == "full") strains else if(associated_strains == "half") sample(strains, 1) else associated_strains)
		AA.freq = get_frequencies(AA.scenario, strains,populations)
		associated_SNPs = if(is.null(associated_SNPs)) SNP.data[,unlist(AA.scenario$associated_SNPs), drop = F] else associated_SNPs
		AA.data = generate_AAs(AA.freq, AA.scenario, study_design, associated_SNPs)
		colnames(AA.data) <- unlist(AA.scenario$id)
		AA.data}))}

get_SNP <- function(study_design, SNP.scenarios) {
	populations = levels(study_design$Population)
	strains = levels(study_design$Strain)
	do.call(cbind, lapply(1:nrow(SNP.scenarios), function(scenario_num) {
		SNP.scenario = SNP.scenarios[scenario_num,]
		SNP.freq = get_frequencies(SNP.scenario, populations, strains)
		SNP.data = generate_SNPs_for_G2G(SNP.freq, study_design)
		colnames(SNP.data) <- unlist(SNP.scenario$id)
		SNP.data}))}

get_frequencies <- function(scenario, main_subgroup, secondary_subgroup) {
	nb_main_subgroup = length(main_subgroup)
	nb_secondary_subgroup = length(secondary_subgroup)
	##Things comapred are on list, that is not natural!!!!
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
##@G
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
	SNP = unlist(lapply(1:nrow(SNP.freq), function(snp_num) {
		p = as.data.frame(lapply(SNP.freq[snp_num,], function(alternate_allele) {c((1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2)}))    
		unlist(lapply(1:nb_populations, function(population_num) {
			iterator = ((population_num - 1 ) * nb_strains) + 1
			unlist(lapply(iterator:(iterator + nb_strains - 1), function(selector) {
				sample(c(0,1,2), size = nb[selector], prob = p[,selector], replace = TRUE) })) })) }))
	matrix(data = SNP ,nrow = nrow(study_design), ncol = nrow(SNP.freq), dimnames = list(rownames(study_design), rownames(SNP.freq)))}

generate_AAs <- function(AA.freq, scenario, study_design, associated_SNPs) {
	strains = levels(study_design$Strain)
	populations = levels(study_design$Population)
	if(trace == TRUE) print(paste(Sys.time(),"Generating Viral output for", scenario$id_tag, sep=" : "))
	res = sapply(1:nrow(AA.freq), function(aa_num) {
		unlist(lapply(populations, function(population) {
			unlist(lapply(strains, function(strain) {
				filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
				AF = AA.freq[aa_num,paste0("F_",strain ,".",population)]
				Y = sample(0:1, size = sum(filter), prob = c(1 - AF,AF), replace = TRUE)
				if(!is.null(associated_SNPs) && population %in% unlist(scenario$Associated_Populations) && strain %in% unlist(scenario$Associated_Strains)) {
					z = associated_SNPs[which(filter),,drop=F]%*%rep(scenario$beta, ncol(associated_SNPs))
					pr = 1/(1+exp(-z))
					Y = Y + unlist(lapply(pr -0.5, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))
					Y[Y>1] = 1}
				Y })) })) })
	rownames(res) <- paste0(rownames(study_design),"_",study_design[,"Population"],"_",study_design[,"Strain"]) 
	res}

analyse_G2G <- function(SNP, AA, study_design, nb_cpu, WO_correction = F, W_human_groups = F, W_strain_group = F, W_both_groups = F, W_human_PC = F, W_strain_PC = F, W_both_PC = F, W_non_linear_PC =F, analyse="logistic") {
	if(trace == T) print(paste0(Sys.time()," : Computing PC"))
	
	SNP_PC = if(W_human_PC || W_both_PC == T) {
		SNP_PC = prcomp(SNP, .scale = F)
		SNP_PC$x[,1:ifelse(ncol(SNP_PC$x)<5, ncol(SNP_PC$x), 5)]}
	
	AA_PC = if(W_strain_PC || W_both_PC == T){
		AA_PC = prcomp(AA, .scale = F)
		AA_PC$x[,1:ifelse(ncol(AA_PC$x)<5, ncol(AA_PC$x), 5)]}
	
	AA_NL_PC = if(W_non_linear_PC == T){
		AA_NL_PC = homals(AA, ndims = 5)
		AA_NL_PC_num = homals(AA_df, level = "numerical")
		matrix(unlist(AA_NL_PC$loadings), nrow = nrow(SNP), ncol = 5)} 
	
	save(AA_PC, SNP_PC, file ="savestates.RData")
	
	logistic_analyse <- function() {
		cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
		
		analyse_AA <- function(Y) {
			Filter(length, list(
				`Without correction` = if(WO_correction == T) apply(SNP, 2, function(X) coef(summary(glm(Y~X)))[,4][2]),
				`With human group` = if(W_human_groups == T) apply(SNP, 2, function(X) coef(summary(glm(Y~X+study_design[,"Population"])))[,4][2]),
				`With strain group` =  if(W_strain_group == T) apply(SNP, 2, function(X) coef(summary(glm(Y~X+study_design[,"Strain"])))[,4][2]),
				`With both groups` = if(W_both_groups == T) apply(SNP, 2, function(X) coef(summary(glm(Y~X+study_design[,"Population"]+study_design[,"Strain"])))[,4][2]),
				`With human PC` = if(W_human_PC == T) apply(SNP, 2, function(X) coef(summary(glm(Y~X+SNP_PC)))[,4][2]),
				`With strain PC` = if(W_strain_PC == T) apply(SNP, 2, function(X) coef(summary(glm(Y~X+AA_PC)))[,4][2]),
				`With both PC` = if(W_both_PC == T) apply(SNP, 2, function(X) coef(summary(glm(Y~X+SNP_PC+AA_PC)))[,4][2]),
				`With non linear PC` = if(W_non_linear_PC == T) apply(SNP, 2, function(X) coef(summary(glm(Y~X+AA_NL_PC)))[,4][2])))}
		
		res = parApply(cl,AA, 2, analyse_AA)
		SNPcol = rep(colnames(SNP), length(res[[1]]))
		CorrectionCol = rep(names(res[[1]]), each = length(res[[1]][[1]]))
		res = do.call(cbind, lapply(names(res), function(aa_id) {
			do.call(rbind, unname(lapply(res[[aa_id]], function(SNP) {
				setNames(data.frame(unname(SNP)), aa_id)})))}))
		cbind(`SNP` = SNPcol, `Correction` =  CorrectionCol, res)}
	
	SKAT_analyse <- function(skatO=F) {
		analyse_AA <- function(Y) {
			Filter(length, list(
				`Without correction` = if(WO_correction == T)	skat_it(SKAT_Null_Model(Y~NULL, out_type = "D")),
				`With human group` = if(W_human_group == T) skat_it(SKAT_Null_Model(Y ~ study_design[,"Population"], out_type = "D")),
				`With strain group` =  if(W_strain_group == T) skat_it(SKAT_Null_Model(Y ~ study_design[,"Strain"], out_type = "D")),
				`With both groups` = if(W_both_groups == T) skat_it(SKAT_Null_Model(Y ~ study_design[,"Population"]+study_design[,"Strain"], out_type = "D")),
				`With human PC` = if(W_strain_group == T) skat_it(SKAT_Null_Model(Y ~ SNP_PC, out_type = "D")),
				`With strain PC` = if(W_strain_PC == T) skat_it(SKAT_Null_Model(Y ~ AA_PC, out_type = "D")),
				`With both PC` = if(W_both_PC == T) skat_it(SKAT_Null_Model(Y ~ AA_PC + SNP_PC, out_type = "D")),
				`With non linear PC` = if(W_non_linear_PC == T) skat_it(SKAT_Null_Model(Y ~ AA_NL_PC, out_type = "D"))))}
		
		flip_SNP <- function(SNP){
			nsample = nrow(SNP)
			apply(SNP,2, function(S) { 
				if(sum(S)/nsample > 1) 
					unlist(lapply(S, function(s) {
						if(s==0) 2 
						else if(s==2) 0 
						else 1}))
				else S })}
		
		SNP = flip_SNP(SNP)
		bio_tag = unlist(SNP.scenarios$bio_tag,recursive = F)
		unique_key_list = setNames(lapply(unique(names(bio_tag)), function (key) which(names(bio_tag) == key)), unique(names(bio_tag)))
		bio_tag = lapply(unique_key_list, function(unique_key_list_index) unname(unlist(lapply(unique_key_list_index, function(index) bio_tag[index]))))
		
		skat_it = {if(!skatO) function(HO) lapply(bio_tag, function(tag) SKAT(SNP[,tag], HO)$p.value)
			else function(HO) lapply(bio_tag, function(tag) SKAT(SNP[,tag], HO, method = "optimal.adj")$p.value)}
		
		cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
		res = parApply(cl,AA, 2, analyse_AA)
		
		SNPcol = rep(names(bio_tag), length(res[[1]]))
		CorrectionCol = rep(names(res[[1]]),  each = length(names(bio_tag)))
		res = do.call(cbind,lapply(res, function(aa) {
			do.call(rbind, lapply(aa, function(correction) {
				do.call(rbind, lapply(correction, function(SNP_tag) {
					SNP_tag}))}))}))
		colnames(res) <- colnames(AA)
		rownames(res) <- NULL
		as.data.frame(cbind(`SNP` = SNPcol, `Correction` = CorrectionCol, res))}
	
	if(trace == T) print(paste(Sys.time(),": Computing ", analyse, "with ",  nb_cpu, " CPU(s)"))
	ptm <- proc.time()
	res = {if(analyse == "logistic") logistic_analyse()
		else if(analyse == "skat") SKAT_analyse()
		else if(analyse == "skato") SKAT_analyse(T)
		else stop(paste(analyse, " is incorrect value for analyse parameter"))}
	save(AA_PC, SNP_PC, res, file ="savestates.RData")
	#save(AA_PC, SNP_PC, res, AA,SNP, AA.scenarios, SNP.scenarios, AA, SNP, file ="savestates.RData")
	if(trace == T) print(paste(Sys.time(),": Analysis took ", (proc.time() - ptm)["elapsed"]/60, "minutes"))
	res
}

plot_full_G2G <- function(res, associations, AA.scenarios, SNP.scenarios) {
	res_tidy = res %>% gather(AA, pvalue, -SNP, -Correction, factor_key = T)
	res_tidy$pvalue = -log10(as.numeric(res_tidy$pvalue))
	threshold = 0.05/(ncol(res)*nrow(res))
	#	res_tidy$SNP_num = rep(1:(length(res_tidy$pvalue)/length(levels(res_tidy$Correction))), length(levels(res_tidy$Correction))) 
	association_table = as.data.frame(do.call(rbind, mapply(function(AA, SNP){
		do.call(rbind, lapply(unlist(AA), function(aa){
			do.call(rbind, lapply(unlist(SNP), function(snp){c(snp, aa)}))}))}
		, AA.scenarios$id, AA.scenarios$associated_SNPs)))
	colnames(association_table) <- c("SNP", "AA")
	association_table$associated= T
	res_tidy$SNP <- as.factor(res_tidy$SNP)
	res_tidy$AA <- as.factor(res_tidy$AA)
	res_tidy = full_join(as.data.frame(association_table), res_tidy, by = c("SNP", "AA"))
}

plot_collapsed_G2G <- function(res, associations, AA.scenarios, SNP.scenarios) {
	res_tidy = res %>% gather(AA, pvalue, -SNP, -Correction, factor_key = T)
	res_tidy$pvalue = -log10(as.numeric(res_tidy$pvalue))
	threshold = 0.05/(ncol(res)*nrow(res))
	#	res_tidy$SNP_num = rep(1:(length(res_tidy$pvalue)/length(levels(res_tidy$Correction))), length(levels(res_tidy$Correction))) 
	association_table = as.data.frame(do.call(rbind, mapply(function(AA, SNP){
		if(!is.na(unlist(SNP))) {
			do.call(rbind, lapply(unlist(AA), function(aa){
				do.call(rbind, lapply(unlist(SNP), function(snp){c(snp, aa)}))}))}},
		AA.scenarios$id, AA.scenarios$associated_SNP_tag)	
	))
	
	colnames(association_table) <- c("SNP", "AA")
	
	association_table$associated <- T
	res_tidy = right_join(as.data.frame(association_table), res_tidy, by = c("SNP", "AA"))
	
	res_tidy$associated[is.na(res_tidy$associated)] <- F
	
	res_tidy$SNP <- as.factor(res_tidy$SNP)
	res_tidy$AA <- as.factor(res_tidy$AA)
	#mapply(function(aa_src, snp_src){}, res_tidy$AA, res_tidy$SNP)
	
	lp = lapply(levels(res_tidy$Correction), function(correction) {
		rt = filter(res_tidy, Correction==correction)
		#rt = filter(res_tidy, Correction=="With strain PC")
		#rt = filter(res_tidy, Correction=="Without correction")
		#rt = filter(res_tidy, Correction=="With both groups")
		
		rt$SNP = as.numeric(rt$SNP)
		rt = arrange(rt,associated)
		p <- ggplot(rt, aes(SNP, pvalue, color=associated))
		p + geom_point() + scale_colour_manual(values =c("black", "red")) + geom_hline(yintercept = -log10(threshold)) + labs(title = correction, x = "SNP")+
			theme(axis.text = element_text(size=24), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36)) +
			scale_y_continuous(limits = c(0, 45))
		ggsave(filename = paste0("file","-",correction,".png"), width = 20, height = 14)
	})}

plot_G2G_by_tag <- function(res, associations, AA.scenarios, SNP.scenarios) {
	#save(res, associations, AA.scenarios, SNP.scenarios, file = "FullScenartio.RData")
	res_tidy = res %>% gather(AA, pvalue, -SNP, -Correction, factor_key = T)
	res_tidy = do.call(rbind, mapply(function(tag,id) {
		data.frame(SNP.tag = as.factor(tag), filter(res_tidy, SNP %in% unlist(id)))}
		, SNP.scenarios$id_tag, SNP.scenarios$id, SIMPLIFY = F))
	res_tidy = do.call(rbind, mapply(function(tag,id) {
		data.frame(AA.tag = as.factor(tag), filter(res_tidy, AA %in% unlist(id)))}
		, AA.scenarios$id_tag, AA.scenarios$id, SIMPLIFY = F))
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
		p <- ggplot(res_tidy, aes(AA.tag, colour=SNP.tag, pvalue, fill=Correction))
		p + geom_boxplot(outlier.color = "grey") + 
			geom_hline(yintercept = -log10(threshold), colour="red") + 
			labs(title = paste("pvalue for full model, with corrections:", paste("", collapse = " ")), y = "-log10(pval)", x="Scenario") + 
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
	
	plot_with_correction_no_ass(save = T)}
#plot_with_correction("WO_correction", "W_human_groups", "W_human_PC", save = T)
#plot_with_correction("WO_correction", "W_strain_group", "W_strain_PC", save = T)
#plot_with_correction("WO_correction", "W_both_groups", "W_both_PC", save = T)
#lapply(levels(res_tidy$Correction), plot_with_correction, save = T)


plot_G2G_exp_by_tag <- function(res, associations, AA.scenarios, SNP.scenarios) {
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

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	library(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
										 ncol = cols, nrow = ceiling(numPlots/cols))
	}
	
	if (numPlots==1) {
		print(plots[[1]])
		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
																			layout.pos.col = matchidx$col))
		}
	}
}