#IMPORTANT: First set the working directory to the path of the G2G-Simulator folder
setwd("/home/onaret/G2G-Simulator")

source("G2G_simulator.R")

library(dplyr)
library(ggplot2)
library(tidyr)
library(parallel)

# Here you can load a data you regenerated using 'G2G-Simulator/paper/generate_paper_dataset.R' or 
# the results in the paper itself 'G2G-simulator-generated-data.RData'available at https://zenodo.org/record/1100977

####IMPORTANT TEMPORARY NOTE: If you use this dataset, please run the following lines of code before.
# A new updated version will be uploaded with a new DOI for the paper before publishing
# Why this lines? A leftover of some tests that are not discussed anymore in the paper need to be removed
scenario5_V5$data$SNP.scenarios = scenario5_V5$data$SNP.scenarios[-(1:100),]
scenario5_V5$data$AA.scenarios = scenario5_V5$data$AA.scenarios[-(1:100),]
scenario5_V5$data$AA.scenarios$associated_SNPs[1:200] = as.list(1:200)

scenario5_V5$data$AA.scenarios$id[1:200] = as.list(1:200)
scenario5_V5$data$AA.scenarios$id[[201]] = c(201:300)
scenario5_V5$data$AA.scenarios$id[[202]] = c(301:400)

scenario5_V5$data$SNP.scenarios$id[1:200] = as.list(1:200)
scenario5_V5$data$SNP.scenarios$id[[201]] = c(201:300)
scenario5_V5$data$SNP.scenarios$id[[202]] = c(301:10300)
scenario5_V5$data$SNP.scenarios$id[[203]] = c(10301:50300)

scenario5_V5$results$logistic = scenario5_V5$results$logistic[-(1:100),-(4:103)]
scenario5_V5$results$logistic = filter(scenario5_V5$results$logistic, !SNP %in% c(1:100))
scenario5_V5$AA.scenarios$bio_tag = gsub("_1","", scenario5_V5$AA.scenarios$bio_tag)
scenario5_V5$SNP.scenarios$bio_tag = gsub("_1","", scenario5_V5$SNP.scenarios$bio_tag)
scenario5_V5$results$logistic$SNP = scenario5_V5$results$logistic$SNP - 100
names(scenario5_V5$results$logistic) <- c("SNP", "SNP_Tag", "Correction", 1:400)

G2G_generated = scenario5_V5
####END

load("G2G-Simulator/data/G2G-simulator-generated-data.RData")

####Make the results plot friendly
results = G2G_generated$results$logistic

threshold = 0.05/((ncol(results)-3)*(nrow(results)/length(levels(results$Correction))))
results = as.data.frame(results) %>% gather(AA, pvalue,-SNP, -Correction, -SNP_Tag, convert = T)

SNP_nb = do.call(sum, as.list(G2G_generated$data$SNP.scenarios))
results = cbind(
  `AA_biotag` = unlist(mapply(function(size, bio_tag) 
    rep(bio_tag, size * SNP_nb * length(levels(results$Correction))), 
    G2G_generated$data$AA.scenarios$size, 
    G2G_generated$data$AA.scenarios$bio_tag)), 
  results)

results$pvalue = -log10(as.numeric(results$pvalue))
association_table = get_association_AA_SNP(G2G_generated$data$AA.scenarios)
results = right_join(as.data.frame(association_table), results, by = c("SNP", "AA"))
results$associated[is.na(results$associated)] <- F

#Aesthetic purpose
results$SNP_Tag = gsub("_1","",results$SNP_Tag)
results$AA_biotag = gsub("_1","",results$AA_biotag)

results$cross_tag = paste0(results$AA_biotag, "|", results$SNP_Tag)
analyse = "logistic"

###Here there is 18 kinds of couples of SNP and AA, some of them show similar signal patten
#If one or the other side is not stratified, there is no false signal
results$cross_tag = gsub("^Asso_Unstratified_AA.*", "Dummy", results$cross_tag)
results$cross_tag = gsub(".*\\|Unstratified_SNP", "Dummy", results$cross_tag)

#Whether it is associated or not is not something we want to see as a different label on plot
#this information will be displayed by dot shape from the 'results$associated' boolean
results$cross_tag = gsub("^Asso_Stratified_Biased_AA_PG2\\|Stratified_SNP", "Stratified_biased_AA\\|Stratified_SNP", results$cross_tag)
results$cross_tag = gsub("^Asso_Stratified_Biased_AA_PG2\\|Stratified_biased_SNP", "Stratified_biased_AA\\|Stratified_biased_SNP", results$cross_tag)

######Boxplot plot 
boxplot_subres <- function(results,title) {
  levels(results$Correction)[1] = "Correcting on both sides"
  levels(results$Correction)[2] = "Correcting on host side"
  levels(results$Correction)[3] = "Correcting on pathogen side"
  levels(results$Correction)[4] = "Without correction"
  cbp <- c("#CC79A7","#E69F00","#0072B2", "#009E73")
  results$Correction = factor(results$Correction,
                        levels = c("Without correction",
                                   "Correcting on host side",
                                   "Correcting on pathogen side",
                                   "Correcting on both sides"), ordered = TRUE)
  p <- ggplot(results, aes(Correction, pvalue))
  p + geom_boxplot(outlier.color = "grey", aes(fill=Correction)) +
    geom_hline(yintercept = -log10(threshold), color = "black", linetype = 2) + 
    labs(title = paste0('boxplot - ', title), y = "Association test pvalue on negative log scale", x="Correction") + 
    scale_fill_manual(values = cbp) +
    coord_cartesian(ylim = c(0, 15)) +
    theme(axis.text = element_blank()) 
  ggsave(filename = paste0('gen-data/boxplot - ', title,".png"), dpi = 300)}

boxplot_subres(filter(results, cross_tag == "Stratified_AA|Stratified_SNP"), "for stratified SNPs and pathogen variants")
boxplot_subres(filter(results, cross_tag == "Stratified_AA|Stratified_biased_SNP"), "for stratified biased SNPs and pathogen variants")
boxplot_subres(filter(results, cross_tag == "Stratified_biased_AA|Stratified_SNP"), "for stratified SNPs and biased pathogen variants")
boxplot_subres(filter(results, cross_tag == "Stratified_biased_AA|Stratified_biased_SNP" & associated == FALSE), "for stratified biased SNPs and biased pathogen variants")
boxplot_subres(filter(results, cross_tag == "Stratified_biased_AA|Stratified_biased_SNP" & associated == TRUE), "for power gain, stratified biased SNPs and biased pathogen variants")
boxplot_subres(filter(results, cross_tag == "Dummy" & associated == TRUE), "for non stratified causal association")
boxplot_subres(filter(results, cross_tag == "Dummy" & associated == FALSE), "for non stratified non causal association")

######General Manhattan plot 
manhattan_plot_subres <- function(results) {
  cl = makeCluster(4, type = "FORK")
  #Spread SNPs along x axis to make it looks better
  results$SNP = rep(sample(1:max(results$SNP), max(results$SNP), replace = F), (max(results$AA))*length(levels(results$Correction)))
  invisible(parLapply(cl, levels(results$Correction), function(correction) {
    results = select(filter(results, Correction==correction), SNP, pvalue, associated, cross_tag)
    cbPalette <- c("#999999", "#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#D55E00")
    ylim = if(correction != "Without correction") c(0, 35) else c(0,90)
    ggsave(
      plot = ggplot(results, aes(SNP, pvalue)) + 
        geom_point(aes(color = cross_tag, shape = associated), size = 3, data = subset(results, (cross_tag != "Dummy"))) + 
        geom_point(aes(color = cross_tag), size = 3, data = subset(results, (cross_tag == "Dummy" & !associated))) + 
        geom_point(aes(color = cross_tag, shape = associated), size = 3, data = rbind(subset(results, (cross_tag == "Dummy"))[1,], subset(results, associated))) +
        geom_hline(yintercept = -log10(threshold), color = "black", linetype = 2) + 
        labs(title = paste0("Manhattan plot - ", correction), x = "SNP", y = "AA association score on negative log scale") + 
        scale_color_manual(values = cbPalette) +
        coord_cartesian(ylim = ylim),
      filename = paste0(getwd(), "/gen-data/manhattan -", correction,".png"), dpi = 300, 
      width = 30, height = 20, units = "cm")}))
  stopCluster(cl)}
manhattan_plot_subres(results)

######Web-app specific Manhattan plot
spec = unique(results$cross_tag)
spec = spec[c(4,1,3,2,5,6)]
plotter = data.frame(`spec` = spec, 
                     `color` = c("#D55E00","#56B4E9","#009E73","#CC79A7","#E69F00","#F0E442"), 
                     `tag`= c("sH_sbP", "sbH_sbP_PG", "no_spurious", "sbH_sbP_FP", "sbH_sP", "sH_sP"))

manhattan_plot_webapp <- function(plotter, results, cor_tag, nb_cpu_outter_loop, nb_cpu_inner_loop) {
  #Spread SNPs along x axis to make it looks better
  results$SNP = rep(sample(1:max(results$SNP), max(results$SNP), replace = F), (max(results$AA))*length(levels(results$Correction)))
  levels(results$Correction) <- c("host_strains", "host", "strains", "no_cor")
  cl1 = makeCluster(nb_cpu_outter_loop, type = "FORK")
  invisible(parLapply(cl1, levels(results$Correction), function(cor_tag) {
    results = filter(results, Correction == cor_tag)
    ylim = if(cor_tag != "no_cor") c(0, 35) else c(0,90)
    cl2 = makeCluster(nb_cpu_inner_loop, type = "FORK")
    invisible(parLapply(cl2, 1:nrow(plotter), function(idx) {
      plot = ggplot(results, aes(SNP, pvalue)) + 
        geom_point(aes(shape = associated), size = 3, color = "#999999", 
                   data = subset(results, (cross_tag != plotter[idx,'spec']))) +
        geom_point(aes(shape = associated), size = 3, color = plotter[idx,'color'], 
                   data = subset(results, (cross_tag == plotter[idx,'spec']))) +
        geom_hline(yintercept = -log10(threshold), color = "black", linetype = 2) + 
        coord_cartesian(ylim = ylim) +
        theme(plot.title = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.position='none')
      ggsave(filename = paste0("gen-data/manhattan-webapp - ",plotter[idx,'tag'],"_",cor_tag,".png"), dpi = 300) }))
    stopCluster(cl2) }))
  stopCluster(cl1) }

manhattan_plot_webapp(plotter, results, "no_cor", nb_cpu_outter_loop = 1, nb_cpu_inner_loop = 6)
manhattan_plot_webapp(plotter, results, "host", nb_cpu_outter_loop = 1, nb_cpu_inner_loop = 6)
manhattan_plot_webapp(plotter, results, "strains", nb_cpu_outter_loop = 1, nb_cpu_inner_loop = 6)
manhattan_plot_webapp(plotter, results, "host_strains", nb_cpu_outter_loop = 1, nb_cpu_inner_loop = 6)