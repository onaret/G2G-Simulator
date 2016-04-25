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

study_design = get_study_design(sample_size=5000, nb_strain=2, nb_pop=2)
nb_cpu = 30
mixstrat = c(0.2)
mixbeta = c(0.25)

##Biotag fs2_ss2

data = {
	if(TRUE){
		parse_G2G_config(study_design,	
										 ###Generic stratification pattern
										 G2G_conf(AA(650), bio_tag = c('min'=5, 'max'=20)),
										 G2G_conf(SNP(10000), bio_tag = c('min'=30, 'max'=60)))}}

data =  
	parse_G2G_config(study_design,	
									 G2G_conf(
									 	association(
									 		SNP(5, stratified = "full", fst_strat = mixstrat), 
									 		AA(1, stratified = "full", fst_strat = mixstrat, associated_strains = "full", associated_populations = "full", beta=mixbeta)),
									 	SNP(1, stratified = "full", fst_strat=mixstrat),
									 	SNP(42),
									 	AA(1, stratified = "full", fst_strat = mixstrat),
									 	AA(7), bio_tag="fs2_ss1", replicate = 5))

data = {
	if(TRUE){
		parse_G2G_config(study_design,	
										 ###Generic stratification pattern
										 G2G_conf(
										 	SNP(1, stratified = "full", fst_strat=mixstrat),
										 	SNP(42), bio_tag = "contain_strat_SNP", replicate = 250),
										 
										 #		G2G_conf(SNP(10000), bio_tag = c('min'=30, 'max'=60)),
										 G2G_conf(SNP(42), bio_tag = "homogenous_tag", replicate = 238),
										 
										 G2G_conf(
										 	AA(1, stratified = "full", fst_strat = mixstrat),
										 	AA(8), bio_tag = "contain_strat_AA", replicate = 50),
										 
										 #		G2G_conf(AA(650), bio_tag = c('min'=5, 'max'=20)),
										 G2G_conf(AA(8), bio_tag = "homogenous_tag", replicate = 81),
										 
										 ##C1 biotag
										 G2G_conf(
										 	SNP(1, biased="full", partial_bias="P1", fst_bias = mixstrat),
										 	SNP(42),
										 	AA(1, stratified = "full", fst_strat = mixstrat, partial_bias="P1", fst_bias = mixstrat),
										 	AA(1, stratified = "full", fst_strat = mixstrat),
										 	AA(7), bio_tag="C1", replicate = 5),
										 
										 ##FS1 biotag
										 G2G_conf(
										 	SNP(1, stratified = "full", fst_strat=mixstrat),
										 	SNP(42),
										 	AA(1, biased = "full", fst_bias = mixstrat),
										 	AA(1, stratified = "full", fst_strat = mixstrat),
										 	AA(7), bio_tag="FS1", replicate = 5),		
										 
										 ##SC1 biotag
										 G2G_conf(
										 	SNP(1, stratified = c("P1","P2"), fst_strat=mixstrat),
										 	SNP(42),
										 	AA(1,  stratified = c("A","B"), fst_strat = mixstrat, biased = c("P1","P2"), fst_bias = mixstrat),
										 	AA(1, stratified = "full", fst_strat = mixstrat),
										 	AA(7), bio_tag="SC1", replicate = 5),
										 
										 ##FS2_SS1
										 G2G_conf(
										 	association(
										 		SNP(5, stratified = "full", fst_strat = mixstrat), 
										 		AA(1, stratified = "full", fst_strat = mixstrat, associated_strains = "full", associated_populations = "full", beta=mixbeta)),
										 	SNP(1, stratified = "full", fst_strat=mixstrat),
										 	SNP(42),
										 	AA(1, stratified = "full", fst_strat = mixstrat),
										 	AA(7), bio_tag="fs2_ss1", replicate = 5),
										 
										 ##FS2_SS1_T1
										 G2G_conf(
										 	association(
										 		SNP(5, stratified = "full", fst_strat = mixstrat), 
										 		AA(2, stratified = "full", fst_strat = mixstrat, associated_strains = "full", associated_populations = "full", beta=mixbeta)),
										 	SNP(1, stratified = "full", fst_strat=mixstrat),
										 	SNP(42),
										 	AA(1, stratified = "full", fst_strat = mixstrat),
										 	AA(7), bio_tag="fs2_ss1_t1", replicate = 5),
										 
										 ###fs2_ss1_t2
										 G2G_conf(
										 	association(
										 		SNP(5, stratified = "full", fst_strat = mixstrat), 
										 		AA(3, stratified = "full", fst_strat = mixstrat, associated_strains = "full", associated_populations = "full", beta=mixbeta)),
										 	SNP(1, stratified = "full", fst_strat=mixstrat),
										 	SNP(42),
										 	AA(1, stratified = "full", fst_strat = mixstrat),
										 	AA(7),bio_tag="fs2_ss1_t2", replicate = 5),
										 
										 ##Biotag fs2_ss2
										 G2G_conf(
										 	association(
										 		SNP(5, stratified = "full", fst_strat = mixstrat),
										 		AA(1, stratified = "full", fst_strat = mixstrat, associated_strains = "full", associated_populations = "half", beta=mixbeta)),
										 	SNP(1, stratified = "full", fst_strat=mixstrat),
										 	SNP(42),
										 	AA(1, stratified = "full", fst_strat = mixstrat),
										 	AA(7), bio_tag="fs2_ss2", replicate = 5),
										 
										 ##BiotagFS3_SS3
										 G2G_conf(
										 	association(
										 		SNP(5, stratified = "full", fst_strat = mixstrat), 
										 		AA(1, stratified = "full", fst_strat = mixstrat, associated_strains = "half", associated_populations = "full", beta=mixbeta)),
										 	SNP(1, stratified = "full", fst_strat=mixstrat),
										 	SNP(42),
										 	AA(1, stratified = "full", fst_strat = mixstrat),
										 	AA(7), bio_tag="fs3_ss3", replicate = 5),	
										 
										 ##BiotagFS4
										 G2G_conf(
										 	association(
										 		SNP(5, stratified = "full", fst_strat = mixstrat), 
										 		AA(1, biased = "full", fst_bias = mixstrat, associated_strains = "half", associated_populations = "full", beta = mixbeta)),
										 	SNP(1, stratified = "full", fst_strat=mixstrat),
										 	SNP(42),
										 	AA(1, stratified = "full", fst_strat = mixstrat),
										 	AA(7), bio_tag="fs4", replicate = 5),
										 
										 ##BiotagSS4
										 G2G_conf(
										 	association(
										 		SNP(5, stratified = "full", fst_strat = mixstrat),
										 		AA(1, associated_strains = "half", associated_populations = "half", beta = mixbeta)),
										 	SNP(1, stratified = "full", fst_strat=mixstrat),
										 	SNP(42),
										 	AA(1, stratified = "full", fst_strat = mixstrat),
										 	AA(7), bio_tag="ss4", replicate = 5),
										 
										 ##BiotagSC1B
										 G2G_conf(
										 	association(
										 		SNP(5, stratified = c("P1","P2"), fst_strat=mixstrat), 
										 		AA(1, stratified = c("A","B"), fst_strat = mixstrat, biased = c("P1","P2"), fst_bias = mixstrat, associated_strains = "A", associated_populations = "full", beta = mixbeta)),
										 	SNP(1, stratified = "full", fst_strat=mixstrat),
										 	SNP(42),
										 	AA(1, stratified = "full", fst_strat = mixstrat),
										 	AA(7), bio_tag="SC1b", replicate = 5)
		)}
	
	############
	else if(FALSE){
		parse_G2G_config(study_design,
										 G2G_conf(
										 	association(
										 		SNP(5, stratified = "full", fst_strat = mixstrat),
										 		AA(1, stratified = "full", fst_strat = mixstrat, associated_strains = "full", associated_populations = "half", beta=mixbeta)),
										 	SNP(1, stratified = "full", fst_strat=mixstrat),
										 	SNP(42),
										 	AA(1, stratified = "full", fst_strat = mixstrat),
										 	AA(7), bio_tag="fs2_ss2", replicate = 5))}
	
	else if(FALSE) {
		parse_G2G_config(study_design,	
										 add_SNP(
										 	SNP_conf(10, biased="full", partial_bias="P1", fst_bias = mixstrat, bio_tag="C1"),
										 	SNP_conf(10, stratified = "full", fst_strat=mixstrat, bio_tag="FS1"),
										 	SNP_conf(10, stratified = c("P1","P2"), fst_strat=mixstrat, bio_tag="SC1"),
										 	SNP_conf(250, stratified = "full", fst_strat=mixstrat),
										 	SNP_conf(10000)
										 ),
										 add_AA(
										 	AA_conf(10, stratified = "full", fst_strat = mixstrat, partial_bias="P1", fst_bias = mixstrat, bio_tag="C1"),
										 	AA_conf(10, biased = "full", fst_bias = mixstrat, bio_tag="FS1"),
										 	AA_conf(10, stratified = c("A","B"), fst_strat = mixstrat, biased = c("P1","P2"), fst_bias = mixstrat, bio_tag="SC1"),
										 	AA_conf(1000),
										 	AA_conf(50, stratified = "full", fst_strat = mixstrat)),
										 association(
										 	add_SNP(SNP_conf(1, stratified = "full", fst_strat = mixstrat, bio_tag="fs2_ss1")), 
										 	add_AA(AA_conf(1, stratified = "full", fst_strat = mixstrat, associated_strains = "full", associated_populations = "full", beta=mixbeta, bio_tag="fs2_ss1")), replicate=5),
										 association(
										 	add_SNP(SNP_conf(5, stratified = "full", fst_strat = mixstrat, bio_tag="fs2_ss1_t1")), 
										 	add_AA(AA_conf(5, stratified = "full", fst_strat = mixstrat, associated_strains = "full", associated_populations = "full", beta=mixbeta, bio_tag="fs2_ss1_t1"))),
										 association(
										 	add_SNP(SNP_conf(5, stratified = "full", fst_strat = mixstrat, bio_tag="fs2_ss1_t2")), 
										 	add_AA(AA_conf(1, stratified = "full", fst_strat = mixstrat, associated_strains = "full", associated_populations = "full", beta=mixbeta, bio_tag="fs2_ss1_t2"))),
										 association(
										 	add_SNP(SNP_conf(1, stratified = "full", fst_strat = mixstrat, bio_tag="fs2_ss2")),
										 	add_AA(AA_conf(1, stratified = "full", fst_strat = mixstrat, associated_strains = "full", associated_populations = "half", beta=mixbeta, bio_tag="fs2_ss2")), replicate=5),
										 association(
										 	add_SNP(SNP_conf(1, stratified = "full", fst_strat = mixstrat, bio_tag="fs3_ss3")), 
										 	add_AA(AA_conf(1, stratified = "full", fst_strat = mixstrat, associated_strains = "half", associated_populations = "full", beta=mixbeta, bio_tag="fs3_ss3")), replicate=5),
										 association(
										 	add_SNP(SNP_conf(1, stratified = "full", fst_strat = mixstrat, bio_tag="fs4")), 
										 	add_AA(AA_conf(1, biased = "full", fst_bias = mixstrat, associated_strains = "half", associated_populations = "full", beta = mixbeta, bio_tag="fs4")), replicate=5),
										 association(
										 	add_SNP(SNP_conf(1, stratified = "full", fst_strat = mixstrat, bio_tag="ss4")),
										 	add_AA(AA_conf(1, associated_strains = "half", associated_populations = "half", beta = mixbeta, bio_tag="ss4")), replicate=5),
										 association(
										 	add_SNP(SNP_conf(10, stratified = c("P1","P2"), fst_strat=mixstrat, bio_tag="SC1b")), 
										 	add_AA(AA_conf(10, stratified = c("A","B"), fst_strat = mixstrat, biased = c("P1","P2"), fst_bias = mixstrat, associated_strains = "A", associated_populations = "full", beta = mixbeta, bio_tag="SC1b"))))}
	
	else if(FALSE){
		parse_G2G_config(study_design,	
										 association(
										 	AA= add_AA(AA_conf(5, associated_strains = "full", associated_populations = "full", beta = 0.1, bio_tag ="AA2a")),
										 	SNP = add_SNP(SNP_conf(1, stratified = "full", fst_strat=c(0.2, 0.02), bio_tag ="SNP2a1"),
										 								SNP_conf(3, stratified = "full", fst_strat=0.2, bio_tag ="SNP2a2")), replicate = 5),
										 add_AA(AA_conf(50)),
										 add_SNP(SNP_conf(100)))}
	
	else if(FALSE){
		parse_G2G_config(study_design,
										 add_SNP(SNP_conf(100)),
										 add_AA(AA_conf(50)))}}

AA = data$AA.data
SNP = data$SNP.data
AA.scenarios = data$AA.scenarios
SNP.scenarios = data$SNP.scenarios
rm(data)
gc()


WO_correction = T
W_human_PC = T
W_strain_PC = T
W_both_PC = T
W_both_groups = T
W_human_group = F
W_strain_group = F
W_non_linear_PC = F

analyse = "skat"

res = analyse_G2G(SNP, AA, study_design, WO_correction = T, W_human_PC = T, W_strain_PC = T, W_both_PC = T, W_both_groups = T, analyse = "skat",nb_cpu)
plot_G2G(res, data.associations, data.AA.scenarios, data.SNP.scenarios)
#save(AA, SNP, AA_PC, SNP_PC, AA.scenarios, SNP.scenarios, file = "savestates.RData")
if(FALSE) {
	t2 =parse_G2G_config(study_design, 
											 association(
											 	add_AA(AA_conf(1, associated_strains = "full", associated_populations = "full", beta = 0.25, bio_tag ="SNP2b_aa")),
											 	add_SNP(SNP_conf(1, stratified = "full", fst_strat=0.2, bio_tag ="SNP2b")), replicate = 2))
	
	data = parse_G2G_config(study_design,
													add_AA(25, fst_bias=0.02, biased = "full"),
													add_AA(25, fst_bias=0.2, biased = "full"),
													add_AA(25, fst_strat=0.02, stratified = "full"),
													add_AA(25, fst_strat=0.2, stratified = "full"),
													add_AA(25, fst_bias=0.02, fst_strat=0.2, stratified = "full", biased = "full"),
													add_AA(25, fst_bias=0.2, fst_strat=0.02, stratified = "full", biased = "full"),
													add_AA(25, fst_bias=0.02, fst_strat=0.2, stratified = "P1", biased = "A"),
													add_AA(25, fst_bias=0.2, fst_strat=0.02, stratified = "P2", biased = "A"),
													add_AA(1, associated_strains = "full", associated_populations = "full", beta = 0.1, associated_SNPs = add_SNP(10, stratified = "full", fst_strat=0.2, bio_tag ="SNP1a")),
													add_AA(1, associated_strains = "full", associated_populations = "full", beta = 0.25, associated_SNPs = add_SNP(10, stratified = "full", fst_strat=0.2, bio_tag ="SNP1b")),
													add_AA(1, associated_strains = "full", associated_populations = "full", beta = 0.5, associated_SNPs = add_SNP(10, stratified = "full", fst_strat=0.2, bio_tag ="SNP1c")),
													add_AA(1, associated_strains = "full", associated_populations = "full", beta = 0.1, associated_SNPs = add_SNP(1, stratified = "full", fst_strat=0.2, bio_tag ="SNP2a"), replicate = 9),
													add_AA(1, associated_strains = "full", associated_populations = "full", beta = 0.25, associated_SNPs = add_SNP(1, stratified = "full", fst_strat=0.2, bio_tag ="SNP2b"), replicate = 9),
													add_AA(1, associated_strains = "full", associated_populations = "full", beta = 0.5, associated_SNPs = add_SNP(1, stratified = "full", fst_strat=0.2, bio_tag ="SNP2c"), replicate = 9),
													add_SNP(500, stratified = "full", fst_strat=0.1),
													add_SNP(10000),
													add_AA(1000))
	
	
	##TODO: check that AA is generated by effect size sum (this is a new function)
	data = parse_G2G_config(study_design,	
													association(
														add_AA(AA_conf(1, associated_strains = "full", associated_populations = "full", beta = 0.1, bio_tag ="AA2a")),
														add_SNP(SNP_conf(1, stratified = "full", fst_strat=0.2, bio_tag ="SNP2a1"),
																		SNP_conf(3, stratified = "full", fst_strat=0.2, bio_tag ="SNP2a2")), replicate = 9),
													association(
														add_AA(AA_conf(1, associated_strains = "full", associated_populations = "full", beta = 0.1, bio_tag ="SNP1a")),
														add_SNP(SNP_conf(10, stratified = "full", fst_strat=0.2, bio_tag ="SNP1a"))),
													association(
														add_AA(AA_conf(1, associated_strains = "full", associated_populations = "full", beta = 0.25, bio_tag ="SNP1b")),
														add_SNP(SNP_conf(10, stratified = "full", fst_strat=0.2, bio_tag ="SNP1b"))),
													association(
														add_AA(AA_conf(1, associated_strains = "full", associated_populations = "full", beta = 0.5, bio_tag ="SNP1c")),
														add_SNP(SNP_conf(10, stratified = "full", fst_strat=0.2, bio_tag ="SNP1c"))),
													
													
													
													association(
														add_AA(AA_conf(1, associated_strains = "full", associated_populations = "full", beta = 0.25, bio_tag ="SNP2b_aa")),
														add_SNP(SNP_conf(1, stratified = "full", fst_strat=0.2, bio_tag ="SNP2b")), replicate = 9),
													association(
														add_AA(AA_conf(1, associated_strains = "full", associated_populations = "full", beta = 0.5)),
														add_SNP(SNP_conf(1, stratified = "full", fst_strat=0.2, bio_tag ="SNP2c")), replicate = 9),
													add_SNP(
														SNP_conf(500, stratified = "full", fst_strat=0.1),
														SNP_conf(10000)),
													add_AA(AA_conf(1000)))}

if(FALSE) {
	AA = data$AA
	SNP = data$SNP
	AA.scenarios = data$AA.scenarios
	SNP.scenarios = data$SNP.scenarios
	associations = data$associations
	rm(data)
	gc()
	res = analyse_G2G(SNP, AA, study_design, nb_cpu)
	plot_G2G(res, data.associations, data.AA.scenarios, data.SNP.scenarios)	}

