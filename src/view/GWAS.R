#parse GWAS result into a friendly format for ggplot
parse_pvalues<- function(pvalues, threshold) {
  pvalues = pvalues[complete.cases(pvalues),, drop = F]
  apply(pvalues, 2 , function(pval) {
    signiff = c()
    signiff[threshold < pval] <- 1
    signiff[threshold*0.1 < pval & pval <= threshold] <- 2
    signiff[threshold*0.01 < pval & pval <= threshold*0.1] <- 3
    signiff[threshold*0.001 < pval & pval <= threshold*0.01] <- 4
    signiff[pval <= threshold*0.001] <- 5
    data.frame(pval, `signiff` = factor(signiff, level=1:5, label=c("ns","*","**","***","+")), `pval_diff` = -log10(pvalues$WO_correction/pval) )})}

summary_sim <- function(res, params,filter="Causal") {
  mapply(function(result, name) {
    FP_neutral = rownames(params[result$signiff != "ns" & !params[,filter],])
    FN_causal = rownames(params[result$signiff == "ns" & params[,filter],])
    chisq <- qchisq(1-result[params[,filter] == FALSE,]$pval,1)
    lambda = median(chisq, na.rm = T)/qchisq(0.5,1)
    names = c("FP_neutral", "FN_causal")
    sum = c(length(FP_neutral), length(FN_causal))
    ratio = (sum/c(sum(!params[,filter]), sum(params[,filter]) ) )*100
    ratio[is.nan(ratio)] = 0
    total = data.frame(sum, ratio, row.names = names)
    power_gain = median(-log10(res$WO_correction$pval/result$pval), na.rm = T)
    res_full = list(`SNP` = c(FP_neutral, FN_causal), `sum`=setNames(sum,names),`ratio`=setNames(ratio, names), `lambda` = lambda)
    res = list(`FP_sum`=sum[1],`FP_ratio`=ratio[1], `FN_sum`=sum[2],`FN_ratio`=ratio[2], `lambda` = lambda, `power_gain` = power_gain)
    Filter(length, res)}, res, names(res))}

###@GWR : GWAS result
plot_GWAS_manhattan <- function(GWR, save = TRUE, file = paste0("gen-data/GWAS-simulation-MH")) {
  if(trace == TRUE) ifelse((save == T), print(paste(Sys.time(),"Saving plot in", paste0("gen-data/"), sep=" : "))) else print(paste(Sys.time(),"Plotting"))
  GWR_name = deparse(substitute(GWR))
  threshold = -log10(0.05/length(GWR$pvalues[[1]]$pval))
  tag = c("Without correction", "With human groups", "With 5PCs")
  names(GWR$pvalues) <- tag
  invisible(mapply(function(main_name, condition) {
    result = cbind(`Stratified` = GWR$SNP_params$Stratified, `Causal` = GWR$SNP_params$Causal,`R` = GWR$SNP_params$R, GWR$pvalues[[condition]])
    result$SNP_type = mapply(function(stratified, causal) {paste(if(stratified) "Stratified" else "Non-Stratified", if(causal) "Causal" else "Non-Causal")}, result$Stratified, result$Causal)
    result$pval = -log10(result$pval)
    result$SNP_num = 1:length(result$pval)
    p <- ggplot(result, aes(SNP_num, pval))
    p + geom_point(aes(colour = SNP_type)) + 
      scale_colour_manual(values =c("orange", "black", "red", "grey")) + 
      geom_hline(yintercept = threshold,colour="red") + labs(title = main_name, x = "SNP") +
      scale_y_continuous(limits = c(0, 45))
    #print(p)
    if(save == TRUE) ggsave(filename = paste0(file,"-",GWR_name,"-",main_name, ".png"), dpi = 300)
  },tag, names(GWR$pvalues)))}

####@pvector: pvalue vector
plot_GWAS_QQ <- function(GWR, save = TRUE, file = paste0("gen-data/GWAS-simulation-QQ")) {
  GWR_name = deparse(substitute(GWR))
  tag = c("Without correction", "With human groups", "With 5PCs")
  names(GWR$pvalues) <- tag
  invisible(mapply(function(main_name, condition) {
    y = -log10(sort(GWR$pvalues[[condition]]$pval, decreasing=F))
    x = -log10( 1:length(y)/length(y) )
    p <- ggplot(data.frame(x,y), aes(x, y))
    p + geom_point() + 
      labs(title=main_name, x = "Expected  -Log10(pval)", y = "Observed -Log10(pval)") + 
      geom_abline(intercept = 0, colour ="red") + 
      scale_x_continuous(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0))
    if(save == TRUE) ggsave(filename = paste0(file,"-",GWR_name,"-",main_name, ".png"), dpi = 300)
  },tag, names(GWR$pvalues)))}

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
    labs(title = paste("FP in function of stratified SNPs number with", fst_strat, "Fst and", nb_pc,"PCs"), 
         y = "False positive per SNP", x=paste("Stratification rate in percentage for", nb_SNP,"SNP"))
  if(save == TRUE) ggsave(filename = file)}