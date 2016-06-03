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
G2G_conf<- function(...,bio_tag=NA, replicate = 1) {
	calls = match.call(expand.dots = FALSE)$`...`
	
	#Should be cleaned, first time using functionnal patern
	res = lapply(paste0(bio_tag,"_", 1:replicate), function(bio_tag) {
		lapply(calls, function(call){
			if(call[1] == "SNP()" || call[1] == "AA()") {
				call$bio_tag = bio_tag
				eval(call)}
			else if(call[1] == "association()"){
				call_mod = lapply(call[2:length(call)], function(ass_call){
					if(ass_call[1] == "SNP()" || ass_call[1] == "AA()") {
						ass_call$bio_tag = bio_tag
						ass_call}
					else(ass_call)})
				do.call(association, call_mod)}
			else {stop(paste0(call[1], " is not a valid function call"))}})})
	##Could merge scenario having same id_tag and different biotag, need to recalculate 'Size', and merge 'id'
	AA.scenarios = rbind_all(lapply(unlist(res,recursive =F), function(res) res$AA.scenarios))
	SNP.scenarios = rbind_all(lapply(unlist(res,recursive =F), function(res) res$SNP.scenarios))
	
	res = list(`AA.scenarios` = AA.scenarios, `SNP.scenarios` = SNP.scenarios)	}

#@`...`: AA() | SNP()
association <- function(..., replicate = 1) {
	saved_eval = substitute(list(...))
	res = replicate(replicate, {
		res = eval(saved_eval)
		AA.scenarios = rbind_all(lapply(res, function(ele) ele$AA.scenarios))
		SNP.scenarios = rbind_all(lapply(res, function(ele) ele$SNP.scenarios))
		#Here make a thirs object, association table instead of putting it into AA.scenario
		AA.scenarios$associated_SNPs = list(unlist(SNP.scenarios$id))
		AA.scenarios$associated_SNP_tag = list(unique(unlist(SNP.scenarios$bio_tag)))
		list(`AA.scenarios` = AA.scenarios, `SNP.scenarios` = SNP.scenarios)}, simplify = F)
	list(`AA.scenarios` = rbind_all(lapply(res, function(ele) ele$AA.scenarios)), `SNP.scenarios` = rbind_all(lapply(res, function(ele) ele$SNP.scenarios)))}

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

analyse_G2G <- function(data, study_design, WO_correction = F, W_human_group = F, W_strain_group = F, W_both_groups = F, W_human_PC = F, W_strain_PC = F, W_both_PC = F, W_strain_groups_human_PC = F, W_non_linear_PC =F, analyse, nb_cpu) {
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
	
	#save(AA_PC, SNP_PC, file ="PC-savestates.RData")
	
	res = lapply(analyse, function(analyse) {
		
		if(trace) print(paste(Sys.time(),": Computing ", analyse, "with ",  nb_cpu, " CPU(s)"))
		ptm <- proc.time()
		
		logistic_analyse <- function() {
			if(trace) print(paste(Sys.time(),": doing logistic regression"))
			analyse_AA <- function(Y) {
				Filter(length, list(
					`Without correction` = if(WO_correction) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X)))[,4][2]),
					`With human group` = if(W_human_group) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+study_design[,"Population"])))[,4][2]),
					`With strain group` =  if(W_strain_group) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+study_design[,"Strain"])))[,4][2]),
					`With both groups` = if(W_both_groups) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+study_design[,"Population"]+study_design[,"Strain"])))[,4][2]),
					`With human PC` = if(W_human_PC) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+SNP_PC)))[,4][2]),
					`With strain PC` = if(W_strain_PC) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+AA_PC)))[,4][2]),
					`With both PC` = if(W_both_PC) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+SNP_PC+AA_PC)))[,4][2]),
					`W_strain_groups_human_PC` = if(W_strain_groups_human_PC) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+SNP_PC+study_design[,"Strain"])))[,4][2]),
					`With non linear PC` = if(W_non_linear_PC) apply(SNP.data, 2, function(X) coef(summary(glm(Y~X+AA_NL_PC+study_design[,"Population"])))[,4][2])))}
			
			cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
			res = parApply(cl,AA.data, 2, analyse_AA)
			stopCluster(cl)
			
			SNPcol = rep(1:ncol(SNP.data), length(res[[1]]))
			SNP_Tag = rep(unlist(mapply(function(tag, size) rep(tag,size), SNP.scenarios$bio_tag, SNP.scenarios$size)), length(res[[1]]))
			CorrectionCol = as.factor(rep(names(res[[1]]), each = ncol(SNP.data)))
			
			res = do.call(cbind, lapply(names(res), function(aa_id) {
				do.call(rbind, unname(lapply(res[[aa_id]], function(SNP) {
					setNames(data.frame(unname(SNP)), aa_id)})))}))
			res = cbind(`SNP_Tag` = SNP_Tag, `Correction` =  CorrectionCol, res)}
		
		SKAT_analyse <- function() {
			if(trace) print(paste(Sys.time(),": doing SKAT"))
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
					`Without correction` = if(WO_correction)	skat_it(SKAT_Null_Model(Y~1, out_type = "D")),
					`With human group` = if(W_human_group) skat_it(SKAT_Null_Model(Y ~ study_design[,"Population"], out_type = "D")),
					`With strain group` =  if(W_strain_group) skat_it(SKAT_Null_Model(Y ~ study_design[,"Strain"], out_type = "D")),
					`With both groups` = if(W_both_groups) skat_it(SKAT_Null_Model(Y ~ study_design[,"Population"]+study_design[,"Strain"], out_type = "D")),
					`With human PC` = if(W_human_PC) skat_it(SKAT_Null_Model(Y ~ SNP_PC, out_type = "D")),
					`With strain PC` = if(W_strain_PC) skat_it(SKAT_Null_Model(Y ~ AA_PC, out_type = "D")),
					`With both PC` = if(W_both_PC) skat_it(SKAT_Null_Model(Y ~ AA_PC + SNP_PC, out_type = "D")),
					`With non linear PC` = if(W_non_linear_PC) skat_it(SKAT_Null_Model(Y ~ AA_NL_PC, out_type = "D"))))}
			
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
			if(trace) print(paste(Sys.time(),": doing gt"))
			analyse_AA <- function(Y) {
				Filter(length, list(
					##Add  data=study_design ?
					`Without correction` = if(WO_correction)	lapply(SNP_batch, function(batch) as.numeric(result(gt(Y, SNP.data[,batch]))["p-value"])),
					`With human group` = if(W_human_group) lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ study_design[,"Population"], SNP.data[,batch]))["p-value"])),
					`With strain group` =  if(W_strain_group) lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ study_design[,"Strain"], SNP.data[,batch]))["p-value"])),
					`With both groups` = if(W_both_groups) lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ study_design[,"Population"]+study_design[,"Strain"], SNP.data[,batch]))["p-value"])),
					`With human PC` = if(W_human_PC) lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ SNP_PC, SNP.data[,batch]))["p-value"])),
					`With strain PC` = if(W_strain_PC) lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ AA_PC, SNP.data[,batch]))["p-value"])),
					`With both PC` = if(W_both_PC) lapply(SNP_batch, function(batch) as.numeric(result(gt(Y ~ SNP_PC + AA_PC, SNP.data[,batch]))["p-value"])),
					`With non linear PC` = if(W_non_linear_PC) lapply(SNP_batch, function(batch) as.numeric(result(gt(Y + AA_NL_PC, SNP.data[,batch]))["p-value"]))))}
			
			bio_tag_key = unique(SNP.scenarios$bio_tag)
			SNP_batch = setNames(lapply(bio_tag_key, function(tag) unlist(filter(SNP.scenarios, bio_tag == tag)$id)), bio_tag_key)
			
			cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
			res = parApply(cl,AA.data, 2, analyse_AA) 		#res = apply(AA, 2, analyse_AA)
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
			analyse_AA <- function(Y) {
				lapply(SNP_batch_id, function(X) G2(AA.data[,Y], SNP.data[,X],nperm = 2500)$std_pval)
			}
			cl = makeCluster(nb_cpu, type = "FORK", outfile='outcluster.log')
			res = parLapply(cl,AA_batch_id, analyse_AA) 		#res = apply(AA, 2, analyse_AA)
			stopCluster(cl)
			res = do.call(cbind, lapply(res, function(AA) do.call(rbind, AA)))
			colnames(res) <- unique(AA.scenarios$bio_tag)
			cbind("SNP_Tag" = paste0("SNP_", unique(SNP.scenarios$bio_tag)), as.data.frame(res))}
		
		res = {if(analyse == "logistic") logistic_analyse()
			else if(grepl("skat", analyse) ) SKAT_analyse()
			else if (analyse == "gt") GT_analyse()
			else if (analyse == "G2") G2_analyse()
			else stop(paste(analyse, " is incorrect value for analyse parameter"))}
		#save(AA_PC, SNP_PC, res, AA.data,SNP.data, AA.scenarios, SNP.scenarios, file ="res-savestates.RData")
		if(trace) print(paste(Sys.time(),": Analysis took ", (proc.time() - ptm)["elapsed"]/60, "minutes"))
		res})
	
	detach(data)
	names(res) <- analyse
	res}

plot_collapsed_G2G <- function(res, AA.scenarios, analyse, file_tag="") {
	invisible(lapply(names(res), function(analyse) {
		if(analyse == "logistic") {
			res = res[[analyse]]
			res$SNP <- 1:nrow(res) 
			threshold = 0.05/((ncol(res)-3)*(nrow(res)/length(levels(res$Correction))))
			res = as.data.frame(res) %>% gather(AA, pvalue,-SNP, -Correction, -SNP_Tag, convert = T)
			res$pvalue = -log10(as.numeric(res$pvalue))
			
			association_table = get_association_AA_SNP(AA.scenarios)
			
			res = right_join(as.data.frame(association_table), res, by = c("SNP", "AA"))
			res$associated[is.na(res$associated)] <- F
			
			invisible(lapply(levels(res$Correction), function(correction) {
				colnames(res) <- gsub("SNP_Tag", "SNP_Gene", colnames(res))
				res = filter(res, Correction==correction)
				res = arrange(res,associated)
				p <- ggplot(res, aes(SNP, pvalue, shape=associated))
				p + geom_point(aes(color=SNP_Gene)) + geom_hline(yintercept = -log10(threshold), color = "red") + labs(title = "G2G Uncollapsed Logistic regression results", x = "SNP")
				ggsave(filename = paste0(getwd(), "/../gen-data/",file_tag,analyse, "-",correction,".png"), dpi = 300)}))}
		
		else if(analyse == "skat-L" | analyse == "skato-L"|analyse == "skato-LW"| analyse == "skat-LW" | analyse == "gt") {
			res = res[[analyse]]
			threshold = 0.05/((ncol(res)-2)*(nrow(res)/length(levels(res$Correction))))
			res = res %>% gather(AA, pvalue, -SNP_Tag, -Correction, convert = T)
			res$pvalue = -log10(res$pvalue)
			association_table = get_association_AA_SNP_tag(AA.scenarios, res)
			association_AA_tag = do.call(rbind, lapply(1:nrow(AA.scenarios), function(num) do.call(rbind, lapply(AA.scenarios[num,]$id, function(aa_id) data_frame("AA" = aa_id, "AA_Tag" = factor(AA.scenarios[num,]$bio_tag, levels = AA.scenarios$bio_tag) )))))
			res = right_join(association_table, res, by = c("SNP_Tag", "AA"))
			res = right_join(association_AA_tag, res, by =c("AA"))
			res$associated[is.na(res$associated)] <- F
			
			invisible(lapply(levels(res$Correction), function(correction) {
				colnames(res) <- gsub("AA_Tag", "AA_Position", colnames(res))
				colnames(res) <- gsub("SNP_Tag", "SNP_Gene", colnames(res))
				res = filter(res, Correction==correction)
				res = arrange(res,associated)
				p <- ggplot(res, aes(SNP_Gene, pvalue, color = AA_Position, size=associated))
				p + geom_point() + geom_hline(yintercept = -log10(threshold), colour = "red") + labs(title = paste0("Collapsed on SNP in human side with ",analyse), x = "SNP Gene")
				ggsave(filename = paste0(getwd(), "/../gen-data/",file_tag,analyse, "-",correction,".png"), dpi = 300)}))}
		
		else if(analyse == "G2"){
			res = res[[analyse]]
			threshold = 0.05/(ncol(res)*nrow(res))
			association_table = data_frame(`SNP_Tag` = res[,"SNP_Tag"], `AA_Tag` = factor(colnames(res[2:ncol(res)])), `associated`=T)
			res = res %>% gather(AA_Tag, pvalue, -SNP_Tag, factor_key = T)
			res$pvalue = -log10(res$pvalue)
			res = right_join(association_table, res, by = c("SNP_Tag", "AA_Tag"))
			#Do not like that
			res$AA_Tag = as.factor(res$AA_Tag)
			res$associated[is.na(res$associated)] <- F
			res = arrange(res,associated)
			colnames(res) <- gsub("AA_Tag", "AA_Position", colnames(res))
			colnames(res) <- gsub("SNP_Tag", "SNP_Gene", colnames(res))
			p <- ggplot(res, aes(SNP_Gene, pvalue, color = AA_Position, size=associated))
			p + geom_point() + geom_hline(yintercept = -log10(threshold), colour = "red") + labs(title = paste0("Collapsed on SNP in human side and AA in viral side with G2"), x = "SNP Gene")
			#			theme(axis.text = element_text(size=12, angle = 90, hjust = 1), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36)) +
			ggsave(filename = paste0(getwd(), "/../gen-data/",file_tag,analyse, ".png"), dpi = 300)}}))}

get_association_AA_SNP_tag <- function(AA.scenario,res) {
	as.data.frame(
		do.call(rbind, mapply(function(AA_id, associated_SNP_tag){
			if(!is.null(unlist(associated_SNP_tag))) {
				do.call(rbind, lapply(unlist(AA_id), function(AA){
					#do.call(rbind, lapply(unlist(associated_SNP_tag), function(Tag){ data.frame(`Tag` = factor(Tag), AA, `associated` = T)}))}))}},
					do.call(rbind, lapply(unlist(associated_SNP_tag), function(SNP_Tag){ data.frame(`SNP_Tag` = factor(SNP_Tag, levels = levels(res$SNP_Tag)), AA, `associated` = T)}))}))}},
			AA.scenario$id, AA.scenario$associated_SNP_tag, SIMPLIFY = F)))}

get_association_AA_SNP <- function(AA.scenario) {
	as.data.frame(do.call(rbind, mapply(function(AA, SNP){
		if(!is.null(unlist(SNP))) {
			do.call(rbind, lapply(unlist(AA), function(AA){
				do.call(rbind, lapply(unlist(SNP), function(SNP){data.frame(SNP, AA, `associated` = T)}))}))}},
		AA.scenario$id, AA.scenario$associated_SNPs, SIMPLIFY = F)))}

plot_collapsed_G2G_on_tags <- function(res, AA.scenarios, SNP.scenarios) {
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

plot_pvalue_by_methods <- function(res, AA.scenarios) {
	analyses = names(res)
	res_all = lapply(analyses, function(analyse){ 
		res = res[[analyse]]
		if(analyse == "logistic") {
			res$SNP <- 1:nrow(res) 
			res = as.data.frame(res) %>% gather(AA, pvalue, -SNP, -Correction, -SNP_Tag, convert = T)
			res$pvalue = -log10(as.numeric(res$pvalue))
			association_table = get_association_AA_SNP(AA.scenarios)
			res = inner_join(as.data.frame(association_table), res, by = c("SNP", "AA"))
			res$framework <- factor("logistic", level = analyses)
			select(res, pvalue, framework)}
		else if(analyse %in% c("skat-L", "skato-L", "gt")){
			res = res %>% gather(AA, pvalue, -SNP_Tag, -Correction, convert = T)
			res$pvalue = -log10(res$pvalue)
			association_table = get_association_AA_SNP_tag(AA.scenarios, res)
			res = inner_join(association_table, res, by = c("SNP_Tag", "AA"))
			res$framework <- factor(analyse, level = analyses)
			select(res, pvalue, framework)}
		else if(analyse == "G2"){
			##This is were the selection occurs, without going through association table
			res = data.frame("pvalue" = -log10(diag(as.matrix(res[2:ncol(res)]))))
			res$framework <- factor("G2", level = analyses)
			res}})
	
	thresholds = get_thresholds(res)
	res_all = do.call(rbind, res_all)
	levels(res_all$framework) <- paste("Method ", c("Logistic Regression", "Global Test", "SKAT", "SKATO", "G2" ))
	
	p <- ggplot(res_all,mapping = aes(framework,pvalue))
	p + geom_boxplot(outlier.color = "grey", aes(fill=framework)) + 
		geom_hline(data = thresholds, aes(yintercept = pvalue, colour = collapsing)) + 
		labs(title = "pvalue distribution for association in function of collapsing or framework", x = "Framework") + 
		theme(axis.text = element_blank())
	ggsave(filename = paste0(getwd(), "/../gen-data/","framework-pvalue-dist",".png"), dpi = 300)}

###Optional
plot_median_pvalue_trend <- function(res) {
	###To implement from prototype script to framework
	res_all_med = do.call(rbind, lapply(1:length(res_test), function(index){
		res_collapsed = lapply(c("skat-L", "skato-L", "gt"), function(analyse) {
			res = res_test[[index]]$res[[analyse]]
			res = res %>% gather(AA, pvalue, -Tag, -Correction, convert = T)
			res$pvalue = -log10(res$pvalue)
			association_table = get_association_AA_SNP_tag(AA.scenarios[[index]], res)
			res_joined = inner_join(association_table, res, by = c("Tag", "AA"))
			ass_num = res_test[[index]]$associated_SNP_nb
			pvalue <- median(res_joined$pvalue)
			data.frame("associated_num" = ass_num, pvalue, "framework" = factor(analyse, level = framework_level))})
		
		res_G2 = {
			res = as.data.frame(res_test[[index]]$res[["G2"]])
			ass_num = res_test[[index]]$associated_SNP_nb
			pvalue <- -log10(median(diag(as.matrix(res[2:ncol(res)]))))
			data.frame("associated_num" = ass_num, pvalue, "framework" = factor("G2", level = framework_level))}
		
		res_logistic = {
			res = res_test[[index]]$res[["logistic"]]
			res = as.data.frame(res) %>% gather(AA, pvalue, -SNP, -Correction, -Tag, convert = T)
			res$pvalue = -log10(as.numeric(res$pvalue))
			association_table = get_association_AA_SNP(AA.scenarios[[index]])
			res_joined = inner_join(as.data.frame(association_table), res, by = c("SNP", "AA"))
			res = data.frame()
			ass_num = res_test[[index]]$associated_SNP_nb
			pvalue <- median(res_joined$pvalue)
			data.frame("associated_num" = ass_num, pvalue, "framework" = factor("logistic", level = framework_level))}
		
		rbind(res_logistic, res_G2, do.call(rbind, res_collapsed))}))
	
	res_all_med = as.data.frame(res_all_med) 
	levels(res_all_med$framework) <- paste("Method ", c("Logistic Regression", "Global Test", "SKAT", "SKATO", "G2" ))
	
	thresholds = get_thresholds(res_test[[1]])
	p <- ggplot(res_all_med, mapping = aes(associated_num, pvalue))
	p + geom_smooth(aes(colour = framework)) + 
		geom_hline(data = thresholds, aes(yintercept = pvalue, colour = collapsing)) + 
		labs(title = "pvalue trend in function of associated SNP number, for the different collapsing or framework", x = "associated SNP") + 
		theme(axis.text = element_text(size=12, angle = 90, hjust = 1), axis.title=element_text(size=18,face="bold"), plot.title = element_text(size = 18))
	ggsave(filename = paste0(getwd(), "/../gen-data/","framework-trend",".png"), dpi = 300)}


get_thresholds <- function(res) {
	threshold_seter = res[["gt"]]
	threshold_1C = -log10(0.05/((ncol(threshold_seter)-2)*(nrow(threshold_seter)/length(levels(threshold_seter$Correction)))))
	
	threshold_seter = res[["G2"]]
	threshold_2C = -log10(0.05/((ncol(threshold_seter)-1)*(nrow(threshold_seter))))
	
	threshold_seter = res[["logistic"]]
	threshold_logistic = -log10(0.05/((ncol(threshold_seter)-3)*(nrow(threshold_seter)/length(levels(threshold_seter$Correction)))))
	gather(data_frame(`Threshold Logistic R` = threshold_logistic, `Threshold Collapsing SNPs` = threshold_1C, `Threshold Collapsing SNPs and AAPV` =threshold_2C), collapsing, pvalue)}
