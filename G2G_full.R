source("G2G.R")

#@`...`: G2G_conf()
parse_G2G_config <- function(study_design, ...) {
	set_env("AA", 0)
	set_env("SNP", 0)
	#Turn it to calls list to set last AA_id and last SNP_id after last G2G_conf call
	s = list(...)
	AA.scenarios = rbind_all(lapply(s, function(ele) ele$AA.scenarios))
	SNP.scenarios = rbind_all(lapply(s, function(ele) ele$SNP.scenarios))
	SNP.data = get_SNP(study_design, SNP.scenarios)
	AA.data = get_AA(study_design, AA.scenarios, SNP.data)
	list(`AA.scenarios` = AA.scenarios,`AA.data` = AA.data, `SNP.data` = SNP.data, `SNP.scenarios` = SNP.scenarios)}

#Here receive last AA_id, last SNP_id, turn replicate to for loops to set last id after last AA or SNP call.
#@`...`: AA() | SNP() | association()
G2G_conf<- function(...,bio_tag, replicate = 1) {
	calls = match.call(expand.dots = FALSE)$`...`
	
	#Acutal working one
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
	AA.scenarios = rbind_all(lapply(unlist(res,recursive =F), function(res) res$AA.scenarios))
	SNP.scenarios = rbind_all(lapply(unlist(res,recursive =F), function(res) res$SNP.scenarios))
	
	res = list(`AA.scenarios` = AA.scenarios, `SNP.scenarios` = SNP.scenarios)	}

#@`...`: AA() | SNP()
association <- function(...) {
	res = list(...)
	AA.scenarios = rbind_all(lapply(res, function(ele) ele$AA.scenarios))
	SNP.scenarios = rbind_all(lapply(res, function(ele) ele$SNP.scenarios))
	AA.scenarios$associated_SNPs = list(unlist(SNP.scenarios$id))
	AA.scenarios$associated_SNP_tag = list(unique(unlist(SNP.scenarios$bio_tag)))
	list(`AA.scenarios` = AA.scenarios, `SNP.scenarios` = SNP.scenarios)	}

AA <- function(size, stratified = NA, partial_strat = NA, fst_strat=NA, biased = NA, partial_bias = NA, fst_bias=NA, associated_strains = NA, associated_populations =NA, beta=NA, bio_tag=NA) {
	list(`AA.scenarios` = do.call(rbind, lapply(fst_strat, function(fst_strat) {
		do.call(rbind, lapply(fst_bias, function(fst_bias) {
			do.call(rbind, lapply(beta, function(beta) {
				id =  get_id("AA",size)
				id_tag = generate_id_tag(fst_strat, partial_strat, fst_bias, partial_bias, beta)
				#bio_tag = if(is.character(bio_tag)) setNames(list(id), bio_tag) else generate_biological_tag(size, bio_tag,id)
				data_frame(Stratified = list(stratified), Biased = list(biased), Partial_Stratification = list(partial_strat), Partial_Bias = list(partial_bias),
									 Associated_Strains= list(associated_strains), Associated_Populations = list(associated_populations), beta, `fst_strat` = fst_strat, `fst_bias` = fst_bias, size, id = list(id), id_tag, bio_tag)}))}))})))}

SNP <- function(size, stratified = NA, partial_strat = NA, fst_strat=NA, biased = NA, partial_bias = NA, fst_bias=NA, bio_tag=NA) {
	list(`SNP.scenarios` = do.call(rbind, lapply(fst_strat, function(fst_strat) {
		do.call(rbind, lapply(fst_bias, function(fst_bias) {
			id =  get_id("SNP",size)
			id_tag = generate_id_tag(fst_strat, partial_strat, fst_bias, partial_bias)
			#bio_tag = if(is.character(bio_tag)) setNames(list(id), bio_tag)  else generate_biological_tag(size, bio_tag, id)
			data_frame(Stratified = list(stratified), Biased = list(biased), Partial_Stratification = list(partial_strat), Partial_Bias = list(partial_bias), `fst_strat` = fst_strat, `fst_bias`= fst_bias, size, id = list(id), id_tag, bio_tag)}))})))}

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
		res = setNames(list(id[1:bio_tag_range[num]]),	paste0(id[1],"_BT_",num,"_",bio_tag_range[num]))
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

analyse_G2G <- function(data, WO_correction = F, W_human_group = F, W_strain_group = F, W_both_groups = F, W_human_PC = F, W_strain_PC = F, W_both_PC = F, W_non_linear_PC =F, analyse="logistic", nb_cpu) {
	#SNP.data = data$SNP.data
	#AA.data = data$AA.data
	#SNP.scenarios = data$SNP.scenarios
	#AA.scenarios = data$AA.scenarios
	#rm(data)
	attach(data)
	
	if(trace) print(paste0(Sys.time()," : Computing PC"))
	
	SNP_PC = if(W_human_PC || W_both_PC) {
		SNP_PC = prcomp(SNP.data, .scale = F)
		SNP_PC$x[,1:ifelse(ncol(SNP_PC$x)<5, ncol(SNP_PC$x), 5)]}
	
	AA_PC = if(W_strain_PC || W_both_PC){
		AA_PC = prcomp(AA.data, .scale = F)
		AA_PC$x[,1:ifelse(ncol(AA_PC$x)<5, ncol(AA_PC$x), 5)]}
	
	AA_NL_PC = if(W_non_linear_PC){
		AA_NL_PC = homals(AA.data, ndims = 5)
		matrix(unlist(AA_NL_PC$loadings), nrow = nrow(SNP), ncol = 5)} 
	
	save(AA_PC, SNP_PC, file ="PC-savestates.RData")
	
	if(trace) print(paste(Sys.time(),": Computing ", analyse, "with ",  nb_cpu, " CPU(s)"))
	ptm <- proc.time()
	
	logistic_analyse <- function() {
		cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
		
		analyse_AA <- function(Y) {
			Filter(length, list(
				`Without correction` = if(WO_correction == T) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X)))[,4][2]),
				`With human group` = if(W_human_group == T) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+study_design[,"Population"])))[,4][2]),
				`With strain group` =  if(W_strain_group == T) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+study_design[,"Strain"])))[,4][2]),
				`With both groups` = if(W_both_groups == T) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+study_design[,"Population"]+study_design[,"Strain"])))[,4][2]),
				`With human PC` = if(W_human_PC == T) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+SNP_PC)))[,4][2]),
				`With strain PC` = if(W_strain_PC == T) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+AA_PC)))[,4][2]),
				`With both PC` = if(W_both_PC == T) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+SNP_PC+AA_PC)))[,4][2]),
				`With non linear PC` = if(W_non_linear_PC == T) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+AA_NL_PC)))[,4][2])))}
		
		res = parApply(cl,AA.data, 2, analyse_AA)
		#nb_covariates = sum(WO_correction, W_human_group, W_strain_group, W_both_groups, W_human_PC, W_strain_PC, W_both_PC, W_non_linear_PC)
		SNPcol = rep(1:ncol(SNP.data), length(res[[1]]))
		bio_tag_col = rep(unlist(mapply(function(tag, size) rep(tag,size), SNP.scenarios$bio_tag, SNP.scenarios$size)), length(res[[1]]))
		CorrectionCol = as.factor(rep(names(res[[1]]), each = ncol(SNP.data)))
		
		res = do.call(cbind, lapply(names(res), function(aa_id) {
			do.call(rbind, unname(lapply(res[[aa_id]], function(SNP) {
				setNames(data.frame(unname(SNP)), aa_id)})))}))
		detach(data)
		res = cbind(`SNP` = SNPcol, `Correction` =  CorrectionCol, `Tag` = bio_tag_col, res)}
	
	SKAT_analyse <- function() {
		skat_it = {
			if(analyse == "skat-LW") function(HO) lapply(SNP_batch, function(batch) SKAT(SNP.data[,batch], HO, kernel = "linear.weighted")$p.value)
			else if(analyse  == "skato-LW") function(HO) lapply(SNP_batch, function(batch) SKAT(SNP.data[,batch], HO, kernel = "linear.weighted", method = "optimal.adj")$p.value)
			else if(analyse  == "skat-L") function(HO) lapply(SNP_batch, function(batch) SKAT(SNP.data[,batch], HO, kernel = "linear")$p.value)
			else if(analyse  == "skato-L") function(HO) lapply(SNP_batch, function(batch) SKAT(SNP.data[,batch], HO, kernel = "linear", method = "optimal.adj")$p.value)
			else stop(paste(analyse, " is incorrect value for analyse parameter"))
		}
		
		analyse_AA <- function(Y) {
			Filter(length, list(
				##Add  data=study_design ?
				`Without correction` = if(WO_correction == T)	skat_it(SKAT_Null_Model(Y~1, out_type = "D")),
				`With human group` = if(W_human_group == T) skat_it(SKAT_Null_Model(Y ~ study_design[,"Population"], out_type = "D")),
				`With strain group` =  if(W_strain_group == T) skat_it(SKAT_Null_Model(Y ~ study_design[,"Strain"], out_type = "D")),
				`With both groups` = if(W_both_groups == T) skat_it(SKAT_Null_Model(Y ~ study_design[,"Population"]+study_design[,"Strain"], out_type = "D")),
				`With human PC` = if(W_human_PC == T) skat_it(SKAT_Null_Model(Y ~ SNP_PC, out_type = "D")),
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
		
		#SNP = flip_SNP(SNP.data)
		bio_tag_key = unique(SNP.scenarios$bio_tag)
		SNP_batch = setNames(lapply(bio_tag_key, function(tag) unlist(filter(SNP.scenarios, bio_tag == tag)$id)), bio_tag_key)
		
		cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
		res = parApply(cl,AA.data, 2, analyse_AA) 		#res = apply(AA, 2, analyse_AA)
		
		tag = rep(names(SNP_batch), length(res[[1]]))
		correction_col = as.factor(rep(names(res[[1]]),  each = length(SNP_batch)))
		res = as.data.frame(do.call(cbind,lapply(res, function(aa) {
			do.call(rbind, lapply(aa, function(correction) {
				do.call(rbind, lapply(correction, function(SNP_tag) {
					SNP_tag}))}))})))
		colnames(res) <- colnames(AA.data)
		rownames(res) <- NULL
		res = as.data.frame(cbind(`Tag` = tag, `Correction` = correction_col, res))}
	
	GT_analyse <- function() {
		
		analyse_AA <- function(Y) {
			Filter(length, list(
				##Add  data=study_design ?
				`Without correction` = if(WO_correction == T)	lapply(SNP_batch, function(batch) as.numeric(result(gt(Y, SNP.data[,batch]))["p-value"])),
				`With human group` = if(W_human_group == T) lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ study_design[,"Population"], SNP.data[,batch]))["p-value"])),
				`With strain group` =  if(W_strain_group == T) lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ study_design[,"Strain"], SNP.data[,batch]))["p-value"])),
				`With both groups` = if(W_both_groups == T) lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ study_design[,"Population"]+study_design[,"Strain"], SNP.data[,batch]))["p-value"])),
				`With human PC` = if(W_human_PC == T) lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ SNP_PC, SNP.data[,batch]))["p-value"])),
				`With strain PC` = if(W_strain_PC == T) lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ AA_PC, SNP.data[,batch]))["p-value"])),
				`With both PC` = if(W_both_PC == T) lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ SNP_PC + AA_PC, SNP.data[,batch]))["p-value"])),
				`With non linear PC` = if(W_non_linear_PC == T) lapply(SNP_batch, function(batch) as.numeric(result(gt(Y + AA_NL_PC, SNP.data[,batch]))["p-value"]))))}
		
		bio_tag_key = unique(SNP.scenarios$bio_tag)
		SNP_batch = setNames(lapply(bio_tag_key, function(tag) unlist(filter(SNP.scenarios, bio_tag == tag)$id)), bio_tag_key)
		
		cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
		res = parApply(cl,AA.data, 2, analyse_AA) 		#res = apply(AA, 2, analyse_AA)
		
		tag = rep(names(SNP_batch), length(res[[1]]))
		correction_col = as.factor(rep(names(res[[1]]),  each = length(SNP_batch)))
		res = as.data.frame(do.call(cbind,lapply(res, function(aa) {
			do.call(rbind, lapply(aa, function(correction) {
				do.call(rbind, lapply(correction, function(SNP_tag) {
					SNP_tag}))}))})))
		colnames(res) <- colnames(AA.data)
		rownames(res) <- NULL
		res = as.data.frame(cbind(`Tag` = tag, `Correction` = correction_col, res))}
	
	res = {if(analyse == "logistic") logistic_analyse()
		else if(grepl("skat", analyse) ) SKAT_analyse()
		else if (analyse == "gt") GT_analyse()}
	#save(AA_PC, SNP_PC, res, AA.data,SNP.data, AA.scenarios, SNP.scenarios, file ="res-savestates.RData")
	if(trace) print(paste(Sys.time(),": Analysis took ", (proc.time() - ptm)["elapsed"]/60, "minutes"))
	detach(data)
	res}

plot_collapsed_G2G <- function(res, SNP.scenarios, AA.scenarios, analyse, file_tag="") {
	
	if(analyse == "logistic") {
		threshold = 0.05/((ncol(res)-3)*(nrow(res)/length(levels(res$Correction))))
		res = res %>% gather(AA, pvalue, -SNP, -Correction, -Tag, convert = T)
		res$pvalue = -log10(res$pvalue)
		
		association_table =	as.data.frame(do.call(rbind, mapply(function(AA, SNP){
			if(!is.null(unlist(SNP))) {
				do.call(rbind, lapply(unlist(AA), function(AA){
					do.call(rbind, lapply(unlist(SNP), function(SNP){data.frame(SNP, AA, `associated` = T)}))}))}},
			AA.scenarios$id, AA.scenarios$associated_SNPs)))
		
		res = right_join(as.data.frame(association_table), res, by = c("SNP", "AA"))
		res$associated[is.na(res$associated)] <- F
		
		lapply(levels(res$Correction), function(correction) {
			rt = filter(res, Correction==correction)
			rt = arrange(rt,associated)
			
			p <- ggplot(rt, aes(SNP, pvalue, shape=associated))
			p + geom_point(aes(color=Tag)) + geom_hline(yintercept = -log10(threshold)) + labs(title = correction, x = "SNP")+
				theme(axis.text = element_text(size=24), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36))# +
			#scale_y_continuous(limits = c(0, 15))
			ggsave(filename = paste0(getwd(), "/../gen-data/",file_tag,"-",analyse, "-",correction,".png"), width = 20, height = 14)})}
	
	else if(analyse == "skat-L" | analyse == "skato-L"|analyse == "skato-LW"| analyse == "skat-LW" | analyse == "gt") {
		threshold = 0.05/((ncol(res)-2)*(nrow(res)/length(levels(res$Correction))))
		res = res %>% gather(AA, pvalue, -Tag, -Correction, convert = T)
		res$pvalue = -log10(res$pvalue)
		
		association_table = as.data.frame(do.call(rbind, 
																							mapply(function(AA_id, associated_SNP_tag){
																								if(!is.null(unlist(associated_SNP_tag))) {
																									do.call(rbind, lapply(unlist(AA_id), function(AA){
																										do.call(rbind, lapply(unlist(associated_SNP_tag), function(Tag){ data.frame(`Tag` = factor(Tag, levels = levels(res$Tag)), AA, `associated` = T)}))}))}},
																								AA.scenarios$id, AA.scenarios$associated_SNP_tag)))
		
		res = right_join(association_table, res, by = c("Tag", "AA"))
		res$associated[is.na(res$associated)] <- F
		
		lapply(levels(res$Correction), function(correction) {
			rt = filter(res, Correction==correction)
			rt = arrange(rt,associated)
			
			p <- ggplot(rt, aes(Tag, pvalue, color=associated))
			p + geom_point() + scale_colour_manual(values =c("black", "red")) + geom_hline(yintercept = -log10(threshold)) + labs(title = correction, x = "SNP")+
				theme(axis.text = element_text(size=12, angle = 90, hjust = 1), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36)) #+
			#scale_y_continuous(limits = c(0, 15))
			ggsave(filename = paste0(getwd(), "/../gen-data/",file_tag,"-",analyse, "-",correction,".png"), width = 20, height = 14)})}}

plot_G2G_by_tag <- function(res, associations, AA.scenarios, SNP.scenarios) {
	#save(res, associations, AA.scenarios, SNP.scenarios, file = "FullScenartio.RData")
	threshold = 0.05/(ncol(res)*nrow(res))
	res_tidy = res %>% gather(AA, pvalue, -SNP, -Correction, factor_key = T)
	res_tidy = do.call(rbind, mapply(function(tag,id) {
		data.frame(SNP.tag = as.factor(tag), filter(res_tidy, SNP %in% unlist(id)))}
		, SNP.scenarios$id_tag, SNP.scenarios$id, SIMPLIFY = F))
	res_tidy = do.call(rbind, mapply(function(tag,id) {
		data.frame(AA.tag = as.factor(tag), filter(res_tidy, AA %in% unlist(id)))}
		, AA.scenarios$id_tag, AA.scenarios$id, SIMPLIFY = F))
	res_tidy$pvalue = -log10(res_tidy$pvalue)
	
	plot_with_correction <- function(...,save = F) {
		correction=c(...)
		p <- ggplot(filter(res_tidy, Correction %in% correction), aes(AA.tag, pvalue, colour=SNP.tag, fill=Correction))
		p + geom_boxplot(outlier.color = "grey") + 
			geom_hline(yintercept = -log10(threshold), colour="red") + 
			labs(title = paste("pvalue for full model, with corrections:", paste(correction, collapse = " ")), y = "-log10(pval)", x="covariate") + 
			scale_y_continuous(limits = c(0, 15)) +
			geom_point(aes(y = pval_ass))
		if(save) ggsave(filename = paste0(getwd(),paste0(correction, collapse = "_"), Sys.time(), ".png"))}
	
	plot_with_correction_no_ass <- function(save = F) {
		p <- ggplot(res_tidy, aes(AA.tag, colour=SNP.tag, pvalue, fill=Correction))
		p + geom_boxplot(outlier.color = "grey") + 
			geom_hline(yintercept = -log10(threshold), colour="red") + 
			labs(title = paste("pvalue for full model, with corrections:", paste("", collapse = " ")), y = "-log10(pval)", x="Scenario") + 
			scale_y_continuous(limits = c(0, 15))
		if(save) ggsave(filename = paste0(getwd(),paste0(collapse = "_"), Sys.time(), ".png"))}}

analyse_and_plot <- function(all_data, WO_correction=F,W_human_group=F,W_strain_group=F,W_both_groups=F,W_human_PC=F,W_strain_PC=F,W_both_PC=F,W_non_linear_PC=F) {
	mapply(function(data, file_tag) {
		analyse = c("logistic", "skat-LW", "skat-L", "skato-LW", "skato-L", "gt")
		lapply(analyse, function(analyse) {
			res = analyse_G2G(data, WO_correction, W_human_group, W_strain_group, W_both_groups, W_human_PC, W_strain_PC, W_both_PC, W_non_linear_PC, analyse, nb_cpu)
			plot_collapsed_G2G(res, data$SNP.scenarios, data$AA.scenarios, analyse, file_tag)	})
		
	}, all_data, names(all_data), SIMPLIFY = F)}