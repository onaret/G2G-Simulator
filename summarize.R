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

summary_sim <- function(pvalues, aa.params) {
  mapply(function(result, name) {
    FP_neutral = rownames(aa.params[result$signiff != "ns" & !aa.params[,"Associated"],])
    FN_causal = rownames(aa.params[result$signiff == "ns" & aa.params[,"Associated"],])
    z = qnorm(result[aa.params[,"Associated"] == FALSE,]$pval/2)
    lambda = median(z^2)/0.456
    names = c("FP_neutral", "FN_causal")
    sum = c(length(FP_neutral), length(FN_causal))
    ratio = (sum/c(sum(!aa.params[,"Associated"]), sum(aa.params[,"Associated"]) ) )*100
    ratio[is.nan(ratio)] = 0
    total = data.frame(sum, ratio, row.names = names)
    power_gain = median(-log10(pvalues$WO_correction$pval/result$pval), na.rm = T)
    res_full = list(`SNP` = c(FP_neutral, FN_causal), `sum`=setNames(sum,names),`ratio`=setNames(ratio, names), `lambda` = lambda)
    res = list(`FP_sum`=sum[1],`FP_ratio`=ratio[1], `FN_sum`=sum[2],`FN_ratio`=ratio[2], `lambda` = lambda, `power_gain` = power_gain)
    print(paste0("For ", name))
    print(paste0(" -Median power gain from uncorrected GLM (-log10(pval)) on FN Causal is ", power_gain[1]))
    print(paste0(" -Ratio of FN causal/Total causal ", ratio[1]))
    #    print(paste0(" -Power Gain on FN Causal is ", power_gain[1], collapse=", "))
    #    print(paste0(" -Ratio is ", ratio, collapse=", "))
    Filter(length, res)}, pvalues, names(pvalues))}

trace_plot <- function(pvalues, save = TRUE, file = paste0(Sys.time(),".png")) {
  if(trace == TRUE) ifelse((save == FALSE), print(paste(Sys.time(),"Ploting...", sep=" : ")), print(paste(Sys.time(),"Writting plots", sep=" : ")))
  mapply(function(main_name, condition) {
    #qq(condition$pval,name)
    result = select(as.data.frame(pvalues), starts_with(condition))
    colnames(result) <- c("pval", "signiff", "dpval")
    result$pval = -log10(result$pval)
    result$SNP = 1:length(result$pval)
    p <- ggplot(result, aes(SNP, pval))
    p + geom_point(aes(colour = signiff)) + scale_colour_manual(values =c("ns"="black", "*"="yellow", "**"="orange", "***"="pink", "+"="red")) + geom_hline(yintercept = -log10(threshold)) + labs(title = main_name, x = "SNP")
    if(save == TRUE) ggsave(filename = paste0(file,"-",main_name,".png"), width = 20, height = 14)
  },c("Without correction", "With human groups", "With ciral groups", "With human and viral groups"), names(pvalues))}

trace_boxplot <- function(res, save = TRUE, out = TRUE, file = paste0(Sys.time(),".png"), prefix="") {
  if(trace == TRUE) ifelse((save == FALSE), print(paste(Sys.time(),"Ploting...", sep=" : ")), print(paste(Sys.time(),"Writting plots", sep=" : "))) 
  #colnames(res) <- c("Without correction", "With human groups", "With viral groups", "With human and viral groups")
  res = melt(res,measure.vars = 1:ncol(res))
  p <- ggplot(res, aes(variable, value))
  p <- if(out == TRUE) {
    # compute lower and upper whiskers
    ylim1 = boxplot.stats(res$value)$stats[c(1, 5)]
    # scale y limits based on ylim1
    p + coord_cartesian(ylim = ylim1*1.05)
  } else p
  p + geom_boxplot(outlier.color = "grey", aes(fill=variable)) + geom_hline(yintercept = -log10(threshold), colour="red") + 
    labs(title = "pvalues with different covariates", y = "-log10(pval)", x="covariate") + 
    theme(axis.text.x=element_blank(), axis.text.y = element_text(size=24), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36)) + 
    if(save == TRUE) ggsave(filename = paste0(file,prefix), width = 14, height = 18)}

trace_boxplotlim <- function(res, save = TRUE, out = TRUE, file = paste0(Sys.time(),".png")) {
  if(trace == TRUE) ifelse((save == FALSE), print(paste(Sys.time(),"Ploting...", sep=" : ")), print(paste(Sys.time(),"Writting plots", sep=" : "))) 
  colnames(res) <- c("Without correction", "With human groups", "With viral groups", "With human and viral groups")
  res = melt(res,measure.vars = 1:4)
  p <- ggplot(res, aes(variable, value))
  p <- if(out == TRUE) {
    # compute lower and upper whiskers
    ylim1 = boxplot.stats(res$value)$stats[c(1, 5)]
    # scale y limits based on ylim1
    p + coord_cartesian(ylim = ylim1*1.05)
  } else p
  p + geom_boxplot(outlier.color = "grey", aes(fill=variable)) + geom_hline(yintercept = -log10(threshold), colour="red") + 
    labs(title = "pvalues with different covariates", y = "-log10(pval)", x="covariate") + 
    theme(axis.text.x=element_blank(), axis.text.y = element_text(size=24), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36)) + 
    scale_y_continuous(limits = c(0, 15)) +
    guides(fill=FALSE)  
  if(save == TRUE) ggsave(filename = file, width = 14, height = 18)}

qq <- function(pvector, title="Quantile-quantile plot of p-values", spartan=F) {
  y = -log10(sort(pvector,decreasing=F))
  x = -log10( 1:length(y)/length(y) )
  p <- ggplot(data.frame(x,y), aes(x, y))
  p + geom_smooth() + labs(title="QQ-plot", x = "Expected  -Log10(pval)", y = "Observed -Log10(pval)") + geom_abline(intercept = 0, colour ="red") + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))}

write <- function(res, tag = NULL, output_dir = getwd()) {
  end = Sys.time()
  dir.create(output_dir, showWarnings = FALSE, recursive = FALSE, mode = "0777")
  trace_plot(res$pvalues, save = TRUE, file = paste0(output_dir, tag, "-plots-(", end, ")"))
  trace_boxplot(-log10(select(as.data.frame(res$pvalues), ends_with("pval"))), save = TRUE, file = paste0(output_dir, tag, "-boxplots-(", end, ").png"))
  trace_boxplotlim(-log10(select(as.data.frame(res$pvalues), ends_with("pval"))), save = TRUE, file = paste0(output_dir, tag, "-boxplots-lim-(", end, ").png"))
  trace_boxplot(select(as.data.frame(res$pvalues), ends_with("pval_diff")), save = TRUE, file = paste0(output_dir, tag, "-boxplots_pvaldiff-(", end, ").png"))
  trace_boxplot(-log10(select(as.data.frame(res$pvalues), ends_with("pval"))),out = FALSE, save = TRUE, file = paste0(output_dir, tag, "-boxplots-out-(", end, ").png"))
  trace_boxplot(select(as.data.frame(res$pvalues), ends_with("pval_diff")),out = FALSE, save = TRUE, file = paste0(output_dir, tag, "-boxplots_pvaldiff-out-(", end, ").png"))
  #  trace_plot(res$pvalues,res$scenarios, save = TRUE, file = paste0(output_dir, tag, "-plots-(", end, ").png"))
  #write.table(x = res$pvalues, file = paste0(output_dir, tag, "-pvalues-(", end, ").csv"))
  #write.table(x = res$SNP_params, file = paste0(output_dir, tag, "-SNP_params-(", end, ").csv"))
  res$SNP_params[,"Associated_Populations"] = vapply(res$SNP_params[,"Associated_Populations"], paste, collapse = ", ", character(1L))
  res$SNP_params[,"Associated_Strains"] = vapply(res$SNP_params[,"Associated_Strains"], paste, collapse = ", ", character(1L))
  write.table(x = cbind(res$SNP_params,res$pvalues), file = paste0(output_dir, tag, "-SNP_params-pvals-(", end, ").csv"))
  write.table(x = res$study_design, file = paste0(output_dir, tag, "-study_design-(", end, ").csv"))
  #write.table(x = res$params, file = paste0(output_dir, tag, "-params-(", end, ").csv"))
  #write.table(x = res$scenario, file = paste0(output_dir, tag, "-scenario-(", end, ").csv"))
  #cat(toJSON(res$summary_sim), file = paste0(output_dir, tag, "-summary-(", end, ").json"))
  write.table(x =res$summary_sim, file = paste0(output_dir, tag, "-summary-(", end, ").csv"))}