###associated_SNPs used for G2G setup

##@G
to_pop_structure <- function(study_design) {
	do.call(rbind, lapply(levels(study_design$Population), function(population) do.call(rbind, lapply(levels(study_design$Strain), function(strain) data_frame(strain,population, `nb` = sum(study_design[,"Strain"] == strain & study_design[,"Population"] == population))))))}

##@G
to_study_design <- function(structure) {
	do.call(rbind, lapply(1:nrow(structure), function(pop_num) {
		do.call(rbind, lapply(1:ncol(structure), function(strain_num) {
			if(structure[strain_num,pop_num] > 0) data.frame(Strain = rep(chartr("123456789", "ABCDEFGHI", strain_num), structure[strain_num,pop_num]), Population = paste0("P",pop_num) ) }))})) }

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

##@G
get_AA <- function(study_design, AA.scenarios, SNP.data, associated_SNPs = NULL, generate_AAs = generate_AAs_custom) {
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
		AA.data = generate_AAs_custom(AA.freq, AA.scenario, study_design, associated_SNPs)
		colnames(AA.data) <- unlist(AA.scenario$id)
		AA.data}))}

##@G
get_SNP <- function(study_design, SNP.scenarios) {
	populations = levels(study_design$Population)
	strains = levels(study_design$Strain)
	do.call(cbind, lapply(1:nrow(SNP.scenarios), function(scenario_num) {
		SNP.scenario = SNP.scenarios[scenario_num,]
		SNP.freq = get_frequencies(SNP.scenario, populations, strains)
		SNP.data = generate_SNPs_for_G2G(SNP.freq, study_design)
		colnames(SNP.data) <- unlist(SNP.scenario$id)
		SNP.data}))}

##@G
get_frequencies <- function(scenario, main_subgroup, secondary_subgroup) {
	nb_main_subgroup = length(main_subgroup)
	nb_secondary_subgroup = length(secondary_subgroup)
	##Things comapred are on list, that is not natural!!!!
	strat = {if(is.na(scenario$Stratified) | scenario$Stratified == "no") rep(0, nb_main_subgroup) 
		else if(scenario$Stratified == "full" | scenario$Stratified == "yes") sample(1:nb_main_subgroup, nb_main_subgroup)
		else rank(unlist(scenario$Stratified))} ###Here we will use user defined order
	###Here we make a filter that will later remove unstratified group, this act on seconday_subgroup
	p_strat = if(!is.na(scenario$Partial_Stratification)) rep(strat, each = nb_secondary_subgroup) * rep(as.numeric(secondary_subgroup == scenario$Partial_Stratification), nb_secondary_subgroup)
	
	bias = {if(is.na(scenario$Biased) | scenario$Biased == "no") rep(0, nb_main_subgroup)
		else if(scenario$Biased == "full" | scenario$Biased == "yes") sample(1:nb_secondary_subgroup, nb_secondary_subgroup)
		else rank(unlist(scenario$Biased))} ###Here we will use user defined order
	###Here we make a filter that will later remove unbiased group, this act on seconday_subgroup
	p_bias = if(!is.na(scenario$Partial_Bias)) rep(bias, nb_main_subgroup) * rep(as.numeric(main_subgroup == scenario$Partial_Bias), each = nb_main_subgroup)
	
	#F_strat = if(length(unique(strat)) >1) sort(get_AF(RF, scenario$fst_strat, length(unique(strat))), decreasing = TRUE)[strat] else rep(get_AF(RF, 0.2,1), nb_main_subgroup)
	###Bias first then strat
	freq = if(!is.null(p_strat)) {
		t(replicate(scenario$size, {
			RF = runif(1, 0.1, 0.4)
			F_bias = if(length(unique(bias)) >1) sort(get_AF(RF, scenario$fst_bias, length(unique(bias))), decreasing = TRUE)[bias] else rep(RF, nb_secondary_subgroup)
			F_strat =  if(length(unique(strat)) >1) c(sapply(F_bias, function(RF) sort(get_AF(RF, scenario$fst_strat, length(unique(strat))), decreasing = TRUE)[strat])) else rep(F_bias, each = nb_main_subgroup)
			#F_strat[!is.null(p_bias) & p_bias == 0] = rep(F_bias, each = nb_secondary_subgroup)[p_bias == 0]
			F_strat[!is.null(p_strat) & p_strat == 0] = rep(F_bias, nb_main_subgroup)[p_strat == 0]
			F_strat}))
	}
	###By default strat first then bias
	else {
		t(replicate(scenario$size, {
			RF = runif(1, 0.1, 0.4)
			F_strat = if(length(unique(strat)) >1) sort(get_AF(RF, scenario$fst_strat, length(unique(strat))), decreasing = TRUE)[strat] else rep(RF, nb_main_subgroup)
			F_bias =  if(length(unique(bias)) >1) c(sapply(F_strat, function(RF) sort(get_AF(RF, scenario$fst_bias, length(unique(bias))), decreasing = TRUE)[bias])) else rep(F_strat, each = nb_secondary_subgroup)
			#F_bias[!is.null(p_strat) & p_strat == 0] = rep(F_strat, nb_main_subgroup)[p_strat == 0]
			F_bias[!is.null(p_bias) & p_bias == 0] = rep(F_strat, each = nb_secondary_subgroup)[p_bias == 0]
			F_bias}))
	}
	colnames(freq) <- paste0("F_", rep(main_subgroup, each = nb_secondary_subgroup),".", secondary_subgroup)
	freq}

#Get alternate frequency according to Nicholas model
get_AF <- function(allele, fst, nb = 1) {
	if(is.na(fst)) stop("Specify fst")
	s1 = allele*(1-fst)/fst
	s2 = (1-allele)*(1-fst)/fst
	f = rbeta(n = nb, shape1 = s1,shape2 = s2)
	#	f[f < 0.01] = 0.01
	#	f
}

generate_SNPs_for_G2G <- function(SNP.freq, study_design) {
	nb_strains = length(levels(study_design$Strain))
	nb_populations = length(levels(study_design$Population))
	#if(trace == TRUE) print(paste(Sys.time(),"Generating SNPs dsitribution with viral properties", sep=" : "))
	#nb = summarize(group_by(study_design,Population, Strain), nb = n())$nb
	nb = unname(unlist(to_pop_structure(study_design)[,"nb"]))
	SNP = unlist(lapply(1:nrow(SNP.freq), function(snp_num) {
		p = as.data.frame(lapply(SNP.freq[snp_num,], function(alternate_allele) {c((1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2)}))    
		unlist(lapply(1:nb_populations, function(population_num) {
			iterator = ((population_num - 1 ) * nb_strains) + 1
			unlist(lapply(iterator:(iterator + nb_strains - 1), function(selector) {
				if(nb[selector]>0) sample(c(0,1,2), size = nb[selector], prob = p[,selector], replace = TRUE) })) })) }))
	matrix(data = SNP ,nrow = nrow(study_design), ncol = nrow(SNP.freq), dimnames = list(rownames(study_design), rownames(SNP.freq)))}

generate_AAs_classic <- function(AA.freq, scenario, study_design, associated_SNPs) {
	strains = levels(study_design$Strain)
	populations = levels(study_design$Population)
	#if(trace == TRUE) print(paste(Sys.time(),"Generating Viral output for", scenario$id_tag, sep=" : "))
	res = sapply(1:nrow(AA.freq), function(aa_num) {
		unlist(lapply(populations, function(population) {
			unlist(lapply(strains, function(strain) {
				filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
				if(!is.null(associated_SNPs) && population %in% unlist(scenario$Associated_Populations) && strain %in% unlist(scenario$Associated_Strains)) {
					z = associated_SNPs[which(filter),,drop=F]%*%rep(scenario$beta, ncol(associated_SNPs))
					pr = 1/(1+exp(-z))
					unlist(lapply(pr, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))}
				else {
					AF = AA.freq[aa_num,paste0("F_",strain ,".",population)]
					sample(0:1, size = sum(filter), prob = c(1 - AF,AF), replace = TRUE)
				}})) })) })
	rownames(res) <- paste0(rownames(study_design),"_",study_design[,"Population"],"_",study_design[,"Strain"]) 
	res}

generate_AAs_custom <- function(AA.freq, scenario, study_design, associated_SNPs) {
	strains = levels(study_design$Strain)
	populations = levels(study_design$Population)
	#if(trace == TRUE) print(paste(Sys.time(),"Generating Viral output for", scenario$id_tag, sep=" : "))
	res = sapply(1:nrow(AA.freq), function(aa_num) {
		unlist(lapply(populations, function(population) {
			unlist(lapply(strains, function(strain) {
				filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
				if(sum(filter)>0) {
					Y_ass =	if(!is.null(associated_SNPs) && population %in% unlist(scenario$Associated_Populations) && strain %in% unlist(scenario$Associated_Strains)) {
						z = associated_SNPs[which(filter),,drop=F]%*%rep(scenario$beta, ncol(associated_SNPs))
						pr = 1/(1+exp(-z))
						unlist(lapply(pr -0.5, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))} 
					else numeric(sum(filter))
					AF = AA.freq[aa_num,paste0("F_",strain ,".",population)]
					Y = Y_ass + sample(0:1, size = sum(filter), prob = c(1 - AF,AF), replace = TRUE)
					Y[Y>1] = 1
					Y	}})) })) })
	rownames(res) <- paste0(rownames(study_design),"_",study_design[,"Population"],"_",study_design[,"Strain"]) 
	res}

generate_AAs_pg_new <-  function(AA.freq, scenario, study_design, associated_SNPs) {
	strains = levels(study_design$Strain)
	populations = levels(study_design$Population)
	#if(trace == TRUE) print(paste(Sys.time(),"Generating Viral output for", scenario$id_tag, sep=" : "))
	res = sapply(1:nrow(AA.freq), function(aa_num) {
		Y_sam = unlist(lapply(populations, function(population) {
			unlist(lapply(strains, function(strain) {
				filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
				AF = AA.freq[aa_num,paste0("F_",strain ,".",population)]
				sample(0:1, size = sum(filter), prob = c(1 - AF,AF), replace = TRUE)})) })) 
		
		Y_ass =	unlist(lapply(populations, function(population) {
			unlist(lapply(strains, function(strain) {
				filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
				if(!is.null(associated_SNPs) && population %in% unlist(scenario$Associated_Populations) && strain %in% unlist(scenario$Associated_Strains)) {
					z = associated_SNPs[which(filter),,drop=F]%*%rep(scenario$beta, ncol(associated_SNPs))
					pr = 1/(1+exp(-z))
					unlist(lapply(pr -0.5, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))}
				else numeric(sum(filter))})) }))
		
		Y = Y_sam + Y_ass
		Y[Y>1] = 1
		nb_to_draw = sum(Y[study_design[,"Population"] == "P1" & study_design[,"Strain"] == "A"]) - sum(Y[study_design[,"Population"] == "P2" & study_design[,"Strain"] == "A"])
		if(nb_to_draw>0) Y[sample(which(study_design[,"Population"] == "P2" & study_design[,"Strain"] == "A" & Y == 0), size =nb_to_draw)] = 1
		else  Y[sample(which(study_design[,"Population"] == "P1" & study_design[,"Strain"] == "A" & Y == 0), size =-nb_to_draw)] = 1
		Y})}#We never removed a Y coming from association, we add random 1 to the group that receive less association because of SNP stratification in host side. The goal is to have same amount of 1 in P1.A and P2.A, we different nature.

generate_AAs_pg_new_2 <-  function(AA.freq, scenario, study_design, associated_SNPs) {
	strains = levels(study_design$Strain)
	populations = levels(study_design$Population)
	#if(trace == TRUE) print(paste(Sys.time(),"Generating Viral output for", scenario$id_tag, sep=" : "))
	res = sapply(1:nrow(AA.freq), function(aa_num) {
		Y_sam = unlist(lapply(populations, function(population) {
			unlist(lapply(strains, function(strain) {
				filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
				AF = AA.freq[aa_num,paste0("F_",strain ,".",population)]
				sample(0:1, size = sum(filter), prob = c(1 - AF,AF), replace = TRUE)})) })) 
		
		Y_ass =	unlist(lapply(populations, function(population) {
			unlist(lapply(strains, function(strain) {
				filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
				if(!is.null(associated_SNPs) && population %in% unlist(scenario$Associated_Populations) && strain %in% unlist(scenario$Associated_Strains)) {
					z = associated_SNPs[which(filter),,drop=F]%*%rep(scenario$beta, ncol(associated_SNPs))
					pr = 1/(1+exp(-z))
					unlist(lapply(pr -0.5, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))}
				else numeric(sum(filter))})) }))
		
		nb_to_draw = sum(Y_ass[study_design[,"Population"] == "P1" & study_design[,"Strain"] == "A"]) - sum(Y_ass[study_design[,"Population"] == "P2" & study_design[,"Strain"] == "A"])
		if(nb_to_draw>0) Y_ass[sample(which(Y_ass[study_design[,"Population"] == "P2" & study_design[,"Strain"] == "A"] == 0), size = nb_to_draw)] = 1
		else  Y_ass[sample(which(Y_ass[study_design[,"Population"] == "P1" & study_design[,"Strain"] == "A"] == 0), size = -nb_to_draw)] = 1
		Y = Y_sam + Y_ass
		Y[Y>1] = 1
		Y})}#We never removed a Y coming from association, we add random 1 to the group that receive less association because of SNP stratification in host side. The goal is to have same amount of 1 in P1.A and P2.A, we different nature.

generate_AAs_pg_new_da <-  function(AA.freq, scenario, study_design, associated_SNPs) {
	strains = levels(study_design$Strain)
	populations = levels(study_design$Population)
	#if(trace == TRUE) print(paste(Sys.time(),"Generating Viral output for", scenario$id_tag, sep=" : "))
	res = sapply(1:nrow(AA.freq), function(aa_num) {
		Y_sam = unlist(lapply(populations, function(population) {
			unlist(lapply(strains, function(strain) {
				filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
				AF = AA.freq[aa_num,paste0("F_",strain ,".",population)]
				sample(0:1, size = sum(filter), prob = c(1 - AF,AF), replace = TRUE)})) })) 
		
		Y_ass =	unlist(lapply(populations, function(population) {
			unlist(lapply(strains, function(strain) {
				filter = study_design[,"Population"] == population & study_design[,"Strain"] == strain
				if(!is.null(associated_SNPs) && population %in% unlist(scenario$Associated_Populations) && strain %in% unlist(scenario$Associated_Strains)) {
					z = associated_SNPs[which(filter),,drop=F]%*%rep(scenario$beta, ncol(associated_SNPs))
					pr = 1/(1+exp(-z))
					unlist(lapply(pr -0.5, function(pri) sample(0:1, 1, prob = c(1-pri, pri))))}
				else numeric(sum(filter))})) }))
		
		Y = Y_sam + Y_ass
		Y[Y>1] = 1
		
		nb_to_draw = sum(Y[study_design[,"Population"] == "P1" & study_design[,"Strain"] == "A"]) - sum(Y[study_design[,"Population"] == "P2" & study_design[,"Strain"] == "A"])
		if(nb_to_draw>0) Y[sample(which(study_design[,"Population"] == "P2" & study_design[,"Strain"] == "A" & Y == 0), size =nb_to_draw)] = 1
		else  Y[sample(which(study_design[,"Population"] == "P1" & study_design[,"Strain"] == "A" & Y == 0), size =-nb_to_draw)] = 1
		
		nb_to_draw = sum(Y[study_design[,"Population"] == "P1" & study_design[,"Strain"] == "B"]) - sum(Y[study_design[,"Population"] == "P2" & study_design[,"Strain"] == "B"])
		if(nb_to_draw>0) Y[sample(which(study_design[,"Population"] == "P2" & study_design[,"Strain"] == "B" & Y == 0), size =nb_to_draw)] = 1
		else  Y[sample(which(study_design[,"Population"] == "P1" & study_design[,"Strain"] == "B" & Y == 0), size =-nb_to_draw)] = 1
		
		Y})}#We never removed a Y coming from association, we add random 1 to the group that receive less association because of SNP stratification in host side. The goal is to have same amount of 1 in P1.A and P2.A, we different nature.
