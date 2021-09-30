source("G2G_simulator.R")

###Option 1: Full G2G simulation
study_design = get_study_design(list(
  `P1` = c(`A` = 1125, `B` = 625), 
  `P2` = c(`A` = 1125, `B`  = 2150)))

G2G_conf = get_G2G_conf(
  association(
    SNP(1),
    AA(1,
       stratified = "full",
       fst_strat = 0.2,
       beta = 0.5),
    replicate = 100),
  SNP(100, stratified = "full", fst_strat = 0.2),
  SNP(800), bio_tag = 'test')

G2G_data =	get_G2G_data(
	study_design,
	G2G_conf)

res = analyse_G2G(G2G_data = G2G_data, 
                  correction = get_correction(WO_correction = T, W_host_PC = T, W_pathogen_PC = T, W_both_PC = T, W_pathogen_group = T, W_pathogen_groups_host_PC = T),
                  nb_cpu = 20)

###Option 2: Single G2G pattern analysis

#Neutral scenarios
study_design = to_study_design(list(`P1` = c(`A` = 2500, `B` = 0), `P2` = c(`A` = 0, `B`  = 2500)))
neutral_A = get_G2G_setup(1000, s_stratified = c("P1","P2"), a_stratified = c("A","B"))
test_G2G_setup(study_design, neutral_A, fst_host_strat = 0.2, fst_pathogen_strat = 0.2, tag = "neutral_A")

study_design = to_study_design(list(`P1` = c(`A` = 1500, `B` = 1000), `P2` = c(`A` = 1000, `B`  = 1500)))
neutral_B = get_G2G_setup(1000, s_stratified = c("P1","P2"), a_stratified = c("A","B"))
test_G2G_setup(study_design, neutral_B, fst_host_strat = 0.2, fst_pathogen_strat = 0.2, tag = "neutral_B")

study_design = to_study_design(list(`P1` = c(`A` = 1250, `B` = 1250), `P2` = c(`A` = 0, `B` = 2500)))
neutral_C = get_G2G_setup(1000, s_stratified = c("P1","P2"), s_biased = c("A","B"), s_partial_bias = c("P1"), a_stratified = c("A","B"))
test_G2G_setup(study_design, neutral_C, fst_host_strat = 0.2, fst_host_bias = 0.2, fst_pathogen_strat = 0.2, tag = "neutral_C")

study_design = to_study_design(list(`P1` = c(`A` = 1500, `B` = 1000), `P2` = c(`A` = 1000, `B`  = 1500)))
neutral_E = get_G2G_setup(100, s_stratified = c("P1","P2"), a_stratified = c("A","B"))
test_G2G_setup(study_design, neutral_E, fst_host_strat = 0.2, fst_pathogen_strat = 0.2, tag = "neutral_E")

#Causal scenarios
study_design = to_study_design(list(`P1` = c(`A` = 2500), `P2` = c(`A` = 2500)))
causal_A = get_G2G_setup(1000, s_stratified = "full", associated_strains = "full", associated_populations = "full", beta = 0.25)
test_G2G_setup(study_design, causal_A, fst_host_strat = 0.2, tag = "causal_A")

study_design = to_study_design(list(`P1` = c(`A` = 5000)))
causal_B = get_G2G_setup(1000, associated_strains = "full", associated_populations = "full", beta = 0.25)
test_G2G_setup(study_design, causal_B,  fst_host_strat = 0.2, tag = "causal_B")

study_design = to_study_design(list(`P1` = c(`A` = 1500, `B` = 1000), `P2` = c(`A` = 1000, `B`  = 1500)))
causal_C = get_G2G_setup(1000, s_stratified = c("P2","P1"), s_biased = c("B","A"), a_stratified = c("A","B"), a_biased = c("P1","P2"), associated_strains = "full", associated_populations = "full", beta = 0.25)
test_G2G_setup(study_design, causal_C, fst_host_strat = 0.2, fst_pathogen_strat = 0.2, fst_host_bias = 0.01, fst_pathogen_bias = 0.01, tag = "causal_C")


###Option 3: GWAS analysis
C1 = generate_population_for_GWAS(list(`P1` = c(`case` = 200, `control` = 400), `P2` = c(`case` = 400, `control`  = 200)))
C2 = generate_population_for_GWAS(list(`P1` = c(`case` = 400, `control` = 200), `P2` = c(`case` = 200, `control` = 400)))
C3 = generate_population_for_GWAS(list(`P1` = c(`case` = 300, `control` = 1), `P2` = c(`case` = 300, `control` = 600)))
C4 = generate_population_for_GWAS(list(`P1` = c(`case` = 300, `control` = 200), `P2` = c(`case` = 200, `control` = 100), `P3` = c(`case` = 100, `control` = 300)))
C5 = generate_population_for_GWAS(list(`P1` = c(`case` = 200, `control` = 1), `P2` = c(`case` = 400, `control` = 200), `P3` = c(`case` = 1, `control` = 400)))

res_C1 = GWAS_scenario(C1, neutral = 100000, neutral_S_rate = 0.05, causal_NS = seq(1,2, by = 0.05), causal_S = seq(1,2, by = 0.05), fst_strat = 0.2)
res_C2 = GWAS_scenario(C2, neutral = 100000, neutral_S_rate = 0.05, causal_NS = seq(1,2, by = 0.05), causal_S = seq(1,2, by = 0.05), fst_strat = 0.2)
res_C3 = GWAS_scenario(C3, neutral = 100000, neutral_S_rate = 0.05, causal_NS = seq(1,2, by = 0.05), causal_S = seq(1,2, by = 0.05), fst_strat = 0.2)
res_C4 = GWAS_scenario(C4, neutral = 100000, neutral_S_rate = 0.05, causal_NS = seq(1,2, by = 0.05), causal_S = seq(1,2, by = 0.05), fst_strat = 0.2)
res_C5 = GWAS_scenario(C5, neutral = 100000, neutral_S_rate = 0.05, causal_NS = seq(1,2, by = 0.05), causal_S = seq(1,2, by = 0.05), fst_strat = 0.2)

plot_GWAS_manhattan(res_C1)
plot_GWAS_manhattan(res_C2)
plot_GWAS_manhattan(res_C3)
plot_GWAS_manhattan(res_C4)
plot_GWAS_manhattan(res_C5)

plot_GWAS_QQ(res_C1)
plot_GWAS_QQ(res_C2)
plot_GWAS_QQ(res_C3)
plot_GWAS_QQ(res_C4)
plot_GWAS_QQ(res_C5)
