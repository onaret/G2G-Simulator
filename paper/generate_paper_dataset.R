source("../G2G_simulator.R")

###Here there causal SNP are straified with PG
scenario5_V5 = analyse_G2G(
  parse_G2G_config(
    to_study_design(list(
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
          associated_strains = "full",
          associated_populations = "full",
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
          associated_strains = "full",
          associated_populations = "full",
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
  correction(WO_correction = T, W_human_PC = T, W_strain_group = T, W_strain_groups_human_PC = T), 
  analyse(logistic = T), 
  nb_cpu = 47)
save(scenario5_V5, file ="scenario5_V5.RData")