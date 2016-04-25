parse_pvalues<- function(pvalues, threshold) {
  apply(pvalues, 2 , function(pval) {
    signiff = c()
    signiff[threshold < pval] <- 1
    signiff[threshold*0.1 < pval & pval <= threshold] <- 2
    signiff[threshold*0.01 < pval & pval <= threshold*0.1] <- 3
    signiff[threshold*0.001 < pval & pval <= threshold*0.01] <- 4
    signiff[pval <= threshold*0.001] <- 5
    data.frame(pval, `signiff` = factor(signiff, level=1:5, label=c("ns","*","**","***","+")), `pval_diff` = -log10(pvalues$WO_correction/pval) )})}

regen_signiff <- function(pvalues, threshold) {
  parse_pvalues(data.frame(`WO_PC` = pvalues[,"WO_PC.pval"], `W_simulated_PC` = pvalues[,"W_simulated_PC.pval"], `W_computed_PC` = pvalues[,"W_computed_PC.pval"], row.names = rownames(pvalues)), threshold)}

summary_sim <- function(pvalues, params, filter) {
  mapply(function(result, name) {
    FP_neutral = rownames(params[result$signiff != "ns" & !params[,filter],])
    FN_causal = rownames(params[result$signiff == "ns" & params[,filter],])
    z = qnorm(result[params[,filter] == FALSE,]$pval/2)
    lambda = median(z^2,na.rm = T)/0.456
    names = c("FP_neutral", "FN_causal")
    sum = c(length(FP_neutral), length(FN_causal))
    ratio = (sum/c(sum(!params[,filter]), sum(params[,filter]) ) )*100
    ratio[is.nan(ratio)] = 0
    total = data.frame(sum, ratio, row.names = names)
    power_gain = median(-log10(pvalues$WO_correction$pval/result$pval), na.rm = T)
    res_full = list(`SNP` = c(FP_neutral, FN_causal), `sum`=setNames(sum,names),`ratio`=setNames(ratio, names), `lambda` = lambda)
    res = list(`FP_sum`=sum[1],`FP_ratio`=ratio[1], `FN_sum`=sum[2],`FN_ratio`=ratio[2], `lambda` = lambda, `power_gain` = power_gain)
    Filter(length, res)}, pvalues, names(pvalues))}

###@GWR : GWAS result
plot_GWAS <- function(GWR, save = TRUE, file = paste0(Sys.time(),".png")) {
  if(trace == TRUE) ifelse((save == FALSE), print(paste(Sys.time(),"Ploting...", sep=" : ")), print(paste(Sys.time(),"Writting plots", sep=" : ")))
  tag = c("Without correction", "With human groups", "With viral groups", "With human and viral groups")
  tag = names(GWR$pvalues)
  mapply(function(main_name, condition) {
    result = cbind(`Stratified` = GWR$SNP_params$Stratified, `Causal` = GWR$SNP_params$Causal,`R` = GWR$SNP_params$R, GWR$pvalues[[condition]])
    result$SNP_type = mapply(function(stratified, causal) {paste(if(stratified) "Stratified" else "Non-Stratified", if(causal) "Causal" else "Non-Causal")}, result$Stratified, result$Causal)
    result$pval = -log10(result$pval)
    result$SNP_num = 1:length(result$pval)
    p <- ggplot(result, aes(SNP_num, pval))
    p + geom_point(aes(colour = SNP_type)) + scale_colour_manual(values =c("orange", "black", "red", "grey")) + geom_hline(yintercept = -log10(threshold)) + labs(title = main_name, x = "SNP")+
      theme(axis.text = element_text(size=24), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36)) +
      scale_y_continuous(limits = c(0, 45))
    if(save == TRUE) ggsave(filename = paste0(file,"-",main_name,".png"), width = 20, height = 14)
  },tag, names(GWR$pvalues))}

###@G2SR : G2G setup results
plot_G2G_setup <- function(G2SR, save = FALSE, out = TRUE, file = Sys.time()) {
  if(trace == TRUE) ifelse((save == FALSE), print(paste(Sys.time(),"Ploting...", sep=" : ")), print(paste(Sys.time(),"Writting plots", sep=" : "))) 
  G2GSR = -log10(select(as.data.frame(res$GWR$pvalues), ends_with("pval")))
  G2GSRD = select(as.data.frame(res$GWR$pvalues), ends_with("pval_diff"))
    mapply(function(res, name) {
      colnames(res) <- c("Without correction", "With human groups", "With viral groups", "With human and viral groups")
      res = melt(res, measure.vars = 1:4)
      p <- ggplot(res, aes(variable, value))
      p <- if(out == TRUE) {
        # compute lower and upper whiskers
        ylim1 = boxplot.stats(res$value)$stats[c(1, 5)]
        # scale y limits based on ylim1
        p + coord_cartesian(ylim = ylim1*1.05)
      } else p
      p + geom_boxplot(outlier.color = "grey", aes(fill=variable)) + geom_hline(yintercept = -log10(threshold), colour="red") + 
        labs(title = paste(title, "pvalues with different covariates"), y = "-log10(pval)", x="covariate") + 
        theme(axis.text.x=element_blank(), axis.text.y = element_text(size=24), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36)) + 
        #if(!is.null(lim)) scale_y_continuous(limits = lim) +
        guides(fill=FALSE)
      if(save == TRUE) ggsave(filename = paste0(name, file), width = 14, height = 18)
    }, list(G2GSR, G2GSRD), c("G2GR", "G2GRD"))
  }

####@pvector: pvalue vector
qq <- function(pvector, title="Quantile-quantile plot of p-values", spartan=F) {
  y = -log10(sort(pvector,decreasing=F))
  x = -log10( 1:length(y)/length(y) )
  p <- ggplot(data.frame(x,y), aes(x, y))
  p + geom_smooth() + labs(title="QQ-plot", x = "Expected  -Log10(pval)", y = "Observed -Log10(pval)") + geom_abline(intercept = 0, colour ="red") + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))}

###@LGWR : List of GWAS results
analyse_FP_in_function_of_s_rate <- function(LGWR, save = TRUE, file = paste0(Sys.time(),".png", title="") ) {
  fst_strat = LGWR[[1]]$params$fst_strat
  nb_pc = LGWR[[1]]$params$nb_pc
  nb_SNP = nrow(LGWR[[1]]$SNP_params)
  data_to_plot = do.call(rbind, lapply(LGWR, function(r) data.frame(`S_rate` = r$params$neutral_S_rate, `FP` = r$summary_sim["FP_sum",])))  
  data_to_plot = data_to_plot  %>% gather(S_rate, FP)
  colnames(data_to_plot) <- c("S_rate", "Correction", "FP")
  data_to_plot$FP <- data_to_plot$FP/nb_SNP
  p <- ggplot(data = data_to_plot, aes(S_rate, FP, colour=Correction))
  p + geom_smooth() +  
    labs(title = paste("FP in function of stratified SNPs number with", fst_strat, "Fst and", nb_pc,"PCs"), y = "False positive per SNP", x=paste("Stratification rate in percentage for", nb_SNP,"SNP"))
 if(save == TRUE) ggsave(filename = file)}

write_G2G_setup <- function(res, tag = NULL, output_dir = getwd()) {
  end = Sys.time()
  dir.create(output_dir, showWarnings = FALSE, recursive = FALSE, mode = "0777")
  trace_plot(res$pvalues, save = TRUE, file = paste0(output_dir, tag, "-plots-(", end, ")"))
  plot_G2G_setup(res, out = TRUE, save = TRUE, file = paste0(output_dir, tag, "-boxplots_out-(", end, ").png"))
  plot_G2G_setup(res, out = FALSE, save = TRUE, file = paste0(output_dir, tag, "-boxplots-(", end, ").png"))
  res$SNP_params[,"Associated_Populations"] = vapply(res$SNP_params[,"Associated_Populations"], paste, collapse = ", ", character(1L))
  res$SNP_params[,"Associated_Strains"] = vapply(res$SNP_params[,"Associated_Strains"], paste, collapse = ", ", character(1L))
  write.table(x = cbind(res$SNP_params,res$pvalues), file = paste0(output_dir, tag, "-SNP_params-pvals-(", end, ").csv"))
  write.table(x = res$study_design, file = paste0(output_dir, tag, "-study_design-(", end, ").csv"))
  write.table(x = res$params, file = paste0(output_dir, tag, "-params-(", end, ").csv"))
  write.table(x = res$scenario, file = paste0(output_dir, tag, "-scenario-(", end, ").csv"))
  write.table(x =res$summary_sim, file = paste0(output_dir, tag, "-summary-(", end, ").csv"))}