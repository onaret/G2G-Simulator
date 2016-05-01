#rm(list = ls())
#cat("\014")
#library(rjson)
library(dplyr)
library(ggplot2)
library(tidyr)
library(parallel)
library(SKAT)

source("G2G_full.R")
source("G2G_single.R")
#source("GWAS.R")
#source("summarize.R")
trace <- TRUE

study_design = get_study_design(sample_size = 5000,
																nb_strain = 2,
																nb_pop = 2)
nb_cpu = 30
mixstrat = c(0.2)
mixbeta = c(0.3)

##Biotag fs2_ss2

data_exp = parse_G2G_config(
	study_design,
	G2G_conf(
		association(
			SNP(1),
			SNP(5, stratified = "full", fst_strat = mixstrat),
			AA(
				1,
				associated_strains = "full",
				associated_populations = "full",
				beta = c(0.25))),
		AA(7),
		SNP(44),
		bio_tag = "fs2_ss1",
		replicate = 5))


data_full = parse_G2G_config(
	study_design,
	G2G_conf(
		association(
			SNP(45),
			AA(
				1,
				associated_strains = "full",
				associated_populations = "full",
				beta = c(0.1))),
		AA(7),
		bio_tag = "fs2_ss1",
		replicate = 20))

data_half = parse_G2G_config(
	study_design,
	G2G_conf(
		association(
			SNP(20),
			AA(
				1,
				associated_strains = "full",
				associated_populations = "full",
				beta = c(0.25))),
		SNP(25),
		AA(7),
		bio_tag = "fs2_ss1",
		replicate = 5))

data_base =	parse_G2G_config(
	study_design,
	G2G_conf(
		association(
			SNP(5),
			AA(
				1,
				associated_strains = "full",
				associated_populations = "full",
				beta = c(0.25))),
		SNP(42),
		AA(7),
		bio_tag = "fs2_ss1",
		replicate = 5))

data_base =	parse_G2G_config(
	study_design,
	G2G_conf(
		association(
			SNP(5, stratified = "full", fst_strat = mixstrat),
			AA(
				1,
				stratified = "full",
				fst_strat = mixstrat,
				associated_strains = "full",
				associated_populations = "full",
				beta = mixbeta)),
		SNP(1, stratified = "full", fst_strat = mixstrat),
		SNP(42),
		AA(1, stratified = "full", fst_strat = mixstrat),
		AA(7),
		bio_tag = "fs2_ss1",
		replicate = 1))

data_imp = 	 
	parse_G2G_config(
		study_design,
		G2G_conf(
			association(
				SNP(5, stratified = "full", fst_strat = mixstrat),
				AA(
					1,
					stratified = "full",
					fst_strat = mixstrat,
					associated_strains = "full",
					associated_populations = "full",
					beta = mixbeta
				)
			),
			
			SNP(1, stratified = "full", fst_strat = mixstrat),
			SNP(42),
			AA(1, stratified = "full", fst_strat = mixstrat),
			AA(7),
			bio_tag = "fs2_ss1",
			replicate = 5
		),
		G2G_conf(
			association(
				SNP(5, stratified = "full", fst_strat = mixstrat),
				AA(
					1,
					stratified = "full",
					fst_strat = mixstrat,
					associated_strains = "full",
					associated_populations = "full",
					beta = mixbeta
				)
			),
			
			SNP(1, stratified = "full", fst_strat = mixstrat),
			SNP(42),
			AA(1, stratified = "full", fst_strat = mixstrat),
			AA(7),
			bio_tag = "fs2_ss_bis",
			replicate = 5
		)
	)



data = {
	##Big Neutral
	if (TRUE) {
		parse_G2G_config(study_design,
										 ###Generic stratification pattern
										 G2G_conf(AA(650), bio_tag = c('min' = 5, 'max' = 20)),
										 G2G_conf(SNP(10000), bio_tag = c('min' = 30, 'max' = 60)))
	}
	
	##Association min
	else if (FALSE) {
		parse_G2G_config(
			study_design,
			G2G_conf(
				association(
					SNP(5, stratified = "full", fst_strat = mixstrat),
					AA(
						1,
						stratified = "full",
						fst_strat = mixstrat,
						associated_strains = "full",
						associated_populations = "full",
						beta = 0.5
					)
				),
				SNP(1, stratified = "full", fst_strat = mixstrat),
				SNP(42),
				AA(1, stratified = "full", fst_strat = mixstrat),
				AA(7),
				bio_tag = "fs2_ss1",
				replicate = 1
			)
		)
	}
	
	##Association
	else if (FALSE) {
		parse_G2G_config(
			study_design,
			G2G_conf(
				association(
					SNP(5, stratified = "full", fst_strat = mixstrat),
					AA(
						1,
						stratified = "full",
						fst_strat = mixstrat,
						associated_strains = "full",
						associated_populations = "full",
						beta = mixbeta
					)
				),
				
				SNP(1, stratified = "full", fst_strat = mixstrat),
				SNP(42),
				AA(1, stratified = "full", fst_strat = mixstrat),
				AA(7),
				bio_tag = "fs2_ss1",
				replicate = 5
			)
		)
	}
	
	else if (FALSE) {
		parse_G2G_config(
			study_design,
			###Generic stratification pattern
			G2G_conf(
				SNP(1, stratified = "full", fst_strat = mixstrat),
				SNP(42),
				bio_tag = "contain_strat_SNP",
				replicate = 250
			),
			
			#      G2G_conf(SNP(10000), bio_tag = c('min'=30, 'max'=60)),
			G2G_conf(SNP(42), bio_tag = "homogenous_tag", replicate = 238),
			
			G2G_conf(
				AA(1, stratified = "full", fst_strat = mixstrat),
				AA(8),
				bio_tag = "contain_strat_AA",
				replicate = 50
			),
			
			#      G2G_conf(AA(650), bio_tag = c('min'=5, 'max'=20)),
			G2G_conf(AA(8), bio_tag = "homogenous_tag", replicate = 81),
			
			##C1 biotag
			G2G_conf(
				SNP(
					1,
					biased = "full",
					partial_bias = "P1",
					fst_bias = mixstrat
				),
				SNP(42),
				AA(
					1,
					stratified = "full",
					fst_strat = mixstrat,
					partial_bias = "P1",
					fst_bias = mixstrat
				),
				AA(1, stratified = "full", fst_strat = mixstrat),
				AA(7),
				bio_tag = "C1",
				replicate = 5
			),
			
			##FS1 biotag
			G2G_conf(
				SNP(1, stratified = "full", fst_strat = mixstrat),
				SNP(42),
				AA(1, biased = "full", fst_bias = mixstrat),
				AA(1, stratified = "full", fst_strat = mixstrat),
				AA(7),
				bio_tag = "FS1",
				replicate = 5
			),
			
			##SC1 biotag
			G2G_conf(
				SNP(
					1,
					stratified = c("P1", "P2"),
					fst_strat = mixstrat
				),
				SNP(42),
				AA(
					1,
					stratified = c("A", "B"),
					fst_strat = mixstrat,
					biased = c("P1", "P2"),
					fst_bias = mixstrat
				),
				AA(1, stratified = "full", fst_strat = mixstrat),
				AA(7),
				bio_tag = "SC1",
				replicate = 5
			),
			
			##FS2_SS1
			G2G_conf(
				association(
					SNP(5, stratified = "full", fst_strat = mixstrat),
					AA(
						1,
						stratified = "full",
						fst_strat = mixstrat,
						associated_strains = "full",
						associated_populations = "full",
						beta = mixbeta
					)
				),
				SNP(1, stratified = "full", fst_strat = mixstrat),
				SNP(42),
				AA(1, stratified = "full", fst_strat = mixstrat),
				AA(7),
				bio_tag = "fs2_ss1",
				replicate = 5
			),
			
			##FS2_SS1_T1
			G2G_conf(
				association(
					SNP(5, stratified = "full", fst_strat = mixstrat),
					AA(
						2,
						stratified = "full",
						fst_strat = mixstrat,
						associated_strains = "full",
						associated_populations = "full",
						beta = mixbeta
					)
				),
				SNP(1, stratified = "full", fst_strat = mixstrat),
				SNP(42),
				AA(1, stratified = "full", fst_strat = mixstrat),
				AA(7),
				bio_tag = "fs2_ss1_t1",
				replicate = 5
			),
			
			###fs2_ss1_t2
			G2G_conf(
				association(
					SNP(5, stratified = "full", fst_strat = mixstrat),
					AA(
						3,
						stratified = "full",
						fst_strat = mixstrat,
						associated_strains = "full",
						associated_populations = "full",
						beta = mixbeta
					)
				),
				SNP(1, stratified = "full", fst_strat = mixstrat),
				SNP(42),
				AA(1, stratified = "full", fst_strat = mixstrat),
				AA(7),
				bio_tag = "fs2_ss1_t2",
				replicate = 5
			),
			
			##Biotag fs2_ss2
			G2G_conf(
				association(
					SNP(5, stratified = "full", fst_strat = mixstrat),
					AA(
						1,
						stratified = "full",
						fst_strat = mixstrat,
						associated_strains = "full",
						associated_populations = "half",
						beta = mixbeta
					)
				),
				SNP(1, stratified = "full", fst_strat = mixstrat),
				SNP(42),
				AA(1, stratified = "full", fst_strat = mixstrat),
				AA(7),
				bio_tag = "fs2_ss2",
				replicate = 5
			),
			
			##BiotagFS3_SS3
			G2G_conf(
				association(
					SNP(5, stratified = "full", fst_strat = mixstrat),
					AA(
						1,
						stratified = "full",
						fst_strat = mixstrat,
						associated_strains = "half",
						associated_populations = "full",
						beta = mixbeta
					)
				),
				SNP(1, stratified = "full", fst_strat = mixstrat),
				SNP(42),
				AA(1, stratified = "full", fst_strat = mixstrat),
				AA(7),
				bio_tag = "fs3_ss3",
				replicate = 5
			),
			
			##BiotagFS4
			G2G_conf(
				association(
					SNP(5, stratified = "full", fst_strat = mixstrat),
					AA(
						1,
						biased = "full",
						fst_bias = mixstrat,
						associated_strains = "half",
						associated_populations = "full",
						beta = mixbeta
					)
				),
				SNP(1, stratified = "full", fst_strat = mixstrat),
				SNP(42),
				AA(1, stratified = "full", fst_strat = mixstrat),
				AA(7),
				bio_tag = "fs4",
				replicate = 5
			),
			
			##BiotagSS4
			G2G_conf(
				association(
					SNP(5, stratified = "full", fst_strat = mixstrat),
					AA(
						1,
						associated_strains = "half",
						associated_populations = "half",
						beta = mixbeta
					)
				),
				SNP(1, stratified = "full", fst_strat = mixstrat),
				SNP(42),
				AA(1, stratified = "full", fst_strat = mixstrat),
				AA(7),
				bio_tag = "ss4",
				replicate = 5
			),
			
			##BiotagSC1B
			G2G_conf(
				association(
					SNP(
						5,
						stratified = c("P1", "P2"),
						fst_strat = mixstrat
					),
					AA(
						1,
						stratified = c("A", "B"),
						fst_strat = mixstrat,
						biased = c("P1", "P2"),
						fst_bias = mixstrat,
						associated_strains = "A",
						associated_populations = "full",
						beta = mixbeta
					)
				),
				SNP(1, stratified = "full", fst_strat = mixstrat),
				SNP(42),
				AA(1, stratified = "full", fst_strat = mixstrat),
				AA(7),
				bio_tag = "SC1b",
				replicate = 5
			)
		)
	}
}

data = data_full
data = data_half
data = data_exp
data = data_base
data=data_imp


AA.data = data$AA.data
SNP.data = data$SNP.data
AA.scenarios = data$AA.scenarios
SNP.scenarios = data$SNP.scenarios
rm(data)
gc()


WO_correction = T
#W_human_PC = T
#W_strain_PC = T
#W_both_PC = T

W_human_PC = F
W_strain_PC = F
W_both_PC = F


W_both_groups = T
W_human_group = F
W_strain_group = F
W_non_linear_PC = F

analyse = "logistic"
res_log = analyse_G2G(SNP.data, AA.data, study_design, SNP.scenarios, AA.scenarios, WO_correction,W_human_group,W_strain_group,W_both_groups,W_human_PC,W_strain_PC,W_both_PC,W_non_linear_PC,analyse,nb_cpu)
plot_collapsed_G2G(res_log, SNP.scenarios, AA.scenarios, analyse, file_tag)

analyse = "skat"
res_skat = analyse_G2G(SNP.data, AA.data, study_design, SNP.scenarios, AA.scenarios, WO_correction,W_human_group,W_strain_group,W_both_groups,W_human_PC,W_strain_PC,W_both_PC,W_non_linear_PC,analyse,nb_cpu)
plot_collapsed_G2G(res_skat, SNP.scenarios, AA.scenarios, analyse, file_tag)

analyse = "skato"
res_skato = analyse_G2G(SNP.data, AA.data, study_design, SNP.scenarios, AA.scenarios, WO_correction,W_human_group,W_strain_group,W_both_groups,W_human_PC,W_strain_PC,W_both_PC,W_non_linear_PC,analyse,nb_cpu)
plot_collapsed_G2G(res_skato, SNP.scenarios, AA.scenarios, analyse, file_tag)



close_me <-function() {
	
	mapply(function(associated_SNP_nb, no_associated_SNP_nb) {
		#		print(paste0("Associated : ", associated_SNP_nb))
		#		associated_SNP_nb <<- associated_SNP_nb
		#		no_associated_SNP_nb <<- no_associated_SNP_nb
		data = eval(substitute(parse_G2G_config(
			study_design,
			G2G_conf(
				association(
					SNP(associated_SNP_nb),
					AA(
						1,
						associated_strains = "full",
						associated_populations = "full",
						beta = c(0.1))),
				AA(10),
				SNP(no_associated_SNP_nb),
				bio_tag = "fs2_ss1",
				replicate = 5),
			G2G_conf(AA(200), bio_tag = c('min' = 5, 'max' = 20)),
			G2G_conf(SNP(1000), bio_tag = c('min' = 30, 'max' = 60))
			
		), list(`associated_SNP_nb` = associated_SNP_nb, `no_associated_SNP_nb` = no_associated_SNP_nb) ) )
		
		AA.data = data$AA.data
		SNP.data = data$SNP.data
		AA.scenarios = data$AA.scenarios
		SNP.scenarios = data$SNP.scenarios
		rm(data)
		
		
		WO_correction = T
		W_human_PC = F
		W_strain_PC = F
		W_both_PC = F
		W_both_groups = F
		W_human_group = F
		W_strain_group = F
		W_non_linear_PC = F
		
		file_tag = paste("ass_SNP_", associated_SNP_nb)
		analyse = "logistic"
		res_log = analyse_G2G(SNP.data, AA.data, study_design, SNP.scenarios, AA.scenarios, WO_correction,W_human_group,W_strain_group,W_both_groups,W_human_PC,W_strain_PC,W_both_PC,W_non_linear_PC,analyse,nb_cpu)
		plot_collapsed_G2G(res_log, SNP.scenarios, AA.scenarios, analyse, file_tag)
		
		analyse = "skat"
		res_skat = analyse_G2G(SNP.data, AA.data, study_design, SNP.scenarios, AA.scenarios, WO_correction,W_human_group,W_strain_group,W_both_groups,W_human_PC,W_strain_PC,W_both_PC,W_non_linear_PC,analyse,nb_cpu)
		plot_collapsed_G2G(res_skat, SNP.scenarios, AA.scenarios, analyse, file_tag)
		
		analyse = "skato"
		res_skato = analyse_G2G(SNP.data, AA.data, study_design, SNP.scenarios, AA.scenarios, WO_correction,W_human_group,W_strain_group,W_both_groups,W_human_PC,W_strain_PC,W_both_PC,W_non_linear_PC,analyse,nb_cpu)
		plot_collapsed_G2G(res_skato, SNP.scenarios, AA.scenarios, analyse, file_tag)
	}
	
	, c(1,5,10,35,40), rev(c(1,5,10,35,40)))
	
	
}
close_me()