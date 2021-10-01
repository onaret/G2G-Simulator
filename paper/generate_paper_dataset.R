#IMPORTANT: First set the working directory to the path of the G2G-Simulator folder
setwd("/home/onaret/G2G-Simulator")
source("G2G_simulator.R")

###Here there causal SNP are stratified with PG
G2G_results = analyse_G2G(
  G2G_data = get_G2G_data(
    study_design = get_study_design(list(
      `P1` = c(`A` = 1500, `B` = 1000), 
      `P2` = c(`A` = 1000, `B` = 1500))),
    G2G_conf(
      association(
        AA(
          size=1,
          stratified = c("A","B"),
          fst_strat = 0.2,
          biased = c("P1","P2"),
          fst_bias = 0.01,
          beta = 0.3,
          bio_tag = "Asso_Stratified_Biased_AA_PG2"),
        SNP(
          size=1,
          stratified = c("P2","P1"),
          biased = c("B","A"),
          fst_strat = 0.2,
          fst_bias = 0.016,
          bio_tag = "Stratified_biased_SNP"),
        replicate = 100),
      association(
        AA(
          size=1,
          beta = 0.3,
          bio_tag = "Asso_Unstratified_AA"),
        SNP(
          size=1,
          bio_tag = "Unstratified_SNP"),
        replicate = 100),
      AA(
        size=100,
        stratified = c("A","B"), 
        biased = c("P1","P2"),
        fst_strat = 0.2,
        fst_bias = 0.005, 
        bio_tag = "Stratified_biased_AA"),
      SNP(
        size=100,
        stratified = c("P1","P2"),
        biased = c("A","B"),
        fst_strat = 0.2,
        fst_bias = 0.016, 
        bio_tag = "Stratified_biased_SNP"),
      AA(
        size=100,
        stratified = c("A","B"), 
        fst_strat = 0.2,
        bio_tag = "Stratified_AA"),
      SNP(
        size=10000,
        stratified = c("P1","P2"),
        fst_strat = 0.2,
        bio_tag = "Stratified_SNP"),
      SNP(
        size = 40000,
        bio_tag = "Unstratified_SNP"))), 
  correction = get_correction(WO_correction = T, W_host_PC = T, W_pathogen_group = T, W_pathogen_groups_host_PC = T), 
  nb_cpu = 47)

save(G2G_results, file ="G2G_results.RData")