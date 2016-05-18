#rm(list = ls())
cat("\014")  
#library(rjson)
library(dplyr)
library(ggplot2)
library(tidyr)
library(parallel)

source("G2G.R")
source("GWAS.R")
trace <- TRUE

###Generate SNPs
if(FALSE) {
	
	####Constants
	C0 = generate_population_for_GWAS(list(`P1` = c(`case` = 1250, `control` = 1250), `P2` = c(`case` = 1250, `control`  = 1250)))
	C1 = generate_population_for_GWAS(list(`P1` = c(`case` = 200, `control` = 400), `P2` = c(`case` = 400, `control`  = 200)))
	C2 = generate_population_for_GWAS(list(`P1` = c(`case` = 400, `control` = 200), `P2` = c(`case` = 200, `control` = 400)))
	C3 = generate_population_for_GWAS(list(`P1` = c(`case` = 300, `control` = 0), `P2` = c(`case` = 300, `control` = 600)))
	C4 = generate_population_for_GWAS(list(`P1` = c(`case` = 300, `control` = 200), `P2` = c(`case` = 200, `control` = 100), `P3` = c(`case` = 100, `control` = 300)))
	C5 = generate_population_for_GWAS(list(`P1` = c(`case` = 200, `control` = 0), `P2` = c(`case` = 400, `control` = 200), `P3` = c(`case` = 0, `control` = 400)))
	AllPop = generate_population_for_GWAS(list(`C1` = C1, `C2`= C2, `C3` = C3, `C4` = C4, `C5` = C5))
	
	GWAS_scenario(C1, 50000, 0.05, causal_NS = seq(1,2, by = 0.05), causal_S = seq(1,2, by = 0.05), 0.2)
	res_GWAS = GWAS_scenario(C1, 50000, 0.05, causal_NS = rep(1.3, 100), causal_S = rep(1.3, 100), 0.2)
	
	res_GWAS_coreval = GWAS_scenario(populations = C1, causal_NS = rep(1.3, 1000), causal_S = rep(1.3, 1000), fst_strat = 0.2)
	
	gwas_s_wo = res_GWAS$pvalues$WO$pval[50101:50200]
	gwas_s_wgr = res_GWAS$pvalues$W_simulated_PC$pval[50101:50200]
	gwas_ns_wo = res_GWAS$pvalues$WO_correction$pval[50001:50100]
	gwas_ns_wgr = res_GWAS$pvalues$W_simulated_PC$pval[50001:50100]
	
	gwas_ns_wo = res_GWAS_coreval$pvalues$WO$pval[1:1000]
	gwas_ns_wgr = res_GWAS_coreval$pvalues$W_simulated_PC$pval[1:1000]
	gwas_s_wo = res_GWAS_coreval$pvalues$WO_correction$pval[1001:2000]
	gwas_s_wgr = res_GWAS_coreval$pvalues$W_simulated_PC$pval[1001:2000]
	
	
	
	#Stratified with group and non stratified with group should be equivalent ideally and they are
	t.test(gwas_s_wgr, gwas_ns_wgr)
	
	#Stratified with group and non stratified without group should be equivalent ideally and they are!
	t.test(gwas_s_wgr, gwas_ns_wo)
	#This reflect a correction on effect size from the groups.
	
	#Proof of overcorection
	t.test(gwas_ns_wo, gwas_ns_wgr)
	
	cl = makeCluster(30, type = "FORK", outfile='outcluster.log')
	res = parLapply(cl,1:30, function(nop){
		res_g = GWAS_scenario(populations = C1, causal_NS = rep(1.3, 10000), causal_S = rep(1.3, 10000), fst_strat = 0.2)
		gwas_ns_wo = res_g$pvalues$WO$pval[1:10000]
		gwas_ns_wgr = res_g$pvalues$W_simulated_PC$pval[1:10000]
		gwas_s_wo = res_g$pvalues$WO_correction$pval[10001:20000]
		gwas_s_wgr = res_g$pvalues$W_simulated_PC$pval[10001:20000]
		c(nop, 
			`gwas_s_wgr-gwas_ns_wgr` = t.test(gwas_s_wgr, gwas_ns_wgr)$p.value, 
			`gwas_s_wgr-gwas_ns_wo` =t.test(gwas_s_wgr, gwas_ns_wo)$p.value,
			`gwas_ns_wo-gwas_ns_wgr` = t.test(gwas_ns_wo, gwas_ns_wgr)$p.value)})  
	stopCluster(cl)
	res_merged = do.call(rbind, res)
	
	cl = makeCluster(30, type = "FORK", outfile='outcluster.log')
	res = clusterCall(cl, GWAS_scenario, `populations` = C1, `causal_NS` = rep(1.3, 10000), `causal_S` = rep(1.3, 10000), `fst_strat` = 0.2)
	stopCluster(cl)
	
	res_p = lapply(res, function(res_g){
		gwas_ns_wo = res_g$pvalues$WO$pval[1:10000]
		gwas_ns_wgr = res_g$pvalues$W_simulated_PC$pval[1:10000]
		gwas_s_wo = res_g$pvalues$WO_correction$pval[10001:20000]
		gwas_s_wgr = res_g$pvalues$W_simulated_PC$pval[10001:20000]
		c(`gwas_s_wgr-gwas_ns_wgr` = t.test(gwas_s_wgr, gwas_ns_wgr)$p.value, 
			`gwas_s_wgr-gwas_ns_wo` = t.test(gwas_s_wgr, gwas_ns_wo)$p.value,
			`gwas_ns_wo-gwas_ns_wgr` = t.test(gwas_ns_wo, gwas_ns_wgr)$p.value)})  
	res_p_merged = do.call(rbind, res_p)
	
	load( "G_Study_SNP_strat_correction4(best).RData")
	cl = makeCluster(25, type = "FORK", outfile='outcluster.log')
	res = parLapply(cl, seq(0, 0.002,length.out = 100), function (strat_rate){
		GWAS_scenario(C1, 50000, strat_rate, causal_S = NULL, seq(1,2, by = 0.05), 0.2)
	})
	nb_pc = 5
	analyse_FP_in_function_of_s_rate(res)
	
	plot_GWAS(res[[1]]$pvalues, res[[1]]$SNP_params, title = paste("neutral rate: ", res[1]$params$neutral_S_rate)) 
	plot_GWAS(res[[2]]$pvalues, res[[2]]$SNP_params, title = paste("neutral rate: ", res[2]$params$neutral_S_rate)) 
	plot_GWAS(res[[3]]$pvalues, res[[3]]$SNP_params, title = paste("neutral rate: ", res[3]$params$neutral_S_rate)) 
	plot_GWAS(res[[4]]$pvalues, res[[4]]$SNP_params, title = paste("neutral rate: ", res[4]$params$neutral_S_rate)) 
	
	load("G_Study_SNP_strat_correction4(best+PC).RData")
	nb_SNP = 50000
	nb_pc = 10
	cl = makeCluster(25, type = "FORK", outfile='outcluster.log')
	res = parLapply(cl, seq(0, 0.002,length.out = 100), function (strat_rate){
		GWAS_scenario(C1, nb_SNP, strat_rate, causal_S = NULL, seq(1,2, by = 0.05), 0.2,nb_pc = nb_pc)
	})
	
	analyse_FP_in_function_of_s_rate(res)
	
	load("5_from_0_to_0.01_strat_10PC_fst0.2.RData")
	nb_SNP = 50000
	nb_pc = 10
	cl = makeCluster(25, type = "FORK", outfile='outcluster.log')
	res = parLapply(cl, seq(0, 0.01,length.out = 500), function (strat_rate){
		GWAS_scenario(C1, nb_SNP, strat_rate, causal_S = NULL, seq(1,2, by = 0.05), 0.2,nb_pc = nb_pc)
	})
	analyse_FP_in_function_of_s_rate(res)
	res_half = lapply(1:100, function(ext) res[[ext]])
	analyse_FP_in_function_of_s_rate(res_half)
	
	####TO do diff PC
	nb_SNP = 50000
	cl = makeCluster(25, type = "FORK", outfile='outcluster.log')
	res = parLapply(cl, seq(1, 100,length.out = 100), function (nb_pc){
		GWAS_scenario(C1, neutral = nb_SNP, neutral_S_rate = 0.0005,  causal_S = NULL, causal_NS = NULL, 0.2,nb_pc = nb_pc)
	})
}