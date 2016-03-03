####Function
multiple_SNP_scenarios <- function(populations, neutral, neutral_S_rat, causal_S=c(), causal_NS=c(), repetitions) {
  #This condition is ugly
  cl <- makeCluster(CPU, outfile = "output.txt", type = "FORK")
  if(is.list(populations) && is.list(populations[[1]])) {
    parLapply(cl, 1:(repetitions*length(populations)), function(repetition, populations, names){
      pop_index = repetition%%length(populations) + 1
      print(paste(Sys.time(),": For population", names[pop_index], "executing repetition", trunc(repetition/length(populations))+1))
      write(SNP_scenario(populations[pop_index], neutral = neutral, neutral_S_rate, causal_S, causal_NS), repetition, names[pop_index])
    },populations, names(populations) )}
  else {
    lapply(1:repetitions, function(repetition){
      if(trace == TRUE) print(paste(Sys.time(),": Executing repetition", repetition))
      write(SNP_scenario(populations, neutral, neutral_S_rate, causal_S, causal_NS), repetition)})}
  if(trace == TRUE) print(paste(Sys.time(),"Done", sep=" : "))
  stopCluster(cl)}

SNP_scenario <- function(populations, neutral=0, neutral_S_rate=0, causal_S=c(), causal_NS=c()) {
  populations = generate_population_structure_for_CC(populations)
  SNP_params = generate_SNPs_frequencies(neutral, neutral_S_rate, causal_S, causal_NS, length(populations))
  SNPs = generate_SNPs(SNP_params, populations)
  pvalues = analyse_SNP(SNPs, populations)
  summary_sim = summary_sim(pvalues, SNP_params)
  trace_plot(pvalues)
  
  populations = rbind(populations, matrix(c(sum(populations[,"Case"]), sum(populations[,"Control"])),ncol = 2, nrow = 1, dimnames = list("Total") ) )
  list(#`SNPs` = SNPs,
    `study_design` = populations,
    `SNP_params` = SNP_params,
    `pvalues` = pvalues,
    `summary_sim` = summary_sim)}

get_stratified_SNPs_qtt <- function(stratification_rate, neutral_SNPs) {
  if(stratification_rate>1) stratification_rate = stratification_rate/100
  if(stratification_rate>100) stratification_rate = 1
  stratification_rate * neutral_SNPs}

generate_population_structure_for_CC <- function(populations) {
  populations = {
    if(!is.list(populations)) {
      possible_nb_pop = populations["min"]:populations["max"]
      size = populations["size"]
      populations = sample(c(possible_nb_pop), size=1, prob=rep(1/length(possible_nb_pop),length(possible_nb_pop)), replace = TRUE)
      populations = min(size, populations)
      populations = runif(populations,0,1)
      populations = round(populations/sum(populations) * size)
      populations[1] = populations[1] + size - sum(populations)
      t(sapply(populations, function(strat_size, pop_id) {
        case_nb = round(strat_size*runif(1,0,1))
        #case_nb = round(strat_size*abs(rnorm(5,0.5,0.2)))
        c(case_nb,(strat_size-case_nb))}))}
    else if(is.null(populations) && is.null(size)) stop("You must set size or population arguments...")
    else t(as.data.frame(populations))}
  colnames(populations) <- c("Case", "Control")
  rownames(populations) <- paste0("P", 1:nrow(populations))
  populations}

generate_SNPs_frequencies <- function(neutral_S, neutral_NS, causal_S, causal_NS, populations) {
  if(trace == TRUE) print(paste(Sys.time(),"Generating allele frequency", sep=" : "))
  AF = {
    if(neutral_S > 0) AF = sapply(1:nrow(populations), function(population) get_alternate_allele(replicate(neutral_S, runif(1,0.01,0.99))))
    if(neutral_NS > 0) AF = rbind(AF, replicate(nrow(populations), get_alternate_allele(replicate(neutral_NS, runif(1,0.4,0.5)))))
    if(length(causal_S) > 0) AF = rbind(AF, sapply(1:nrow(populations), function(population) get_alternate_allele(replicate(length(causal_S), runif(1,0.01,0.99)))))
    if(length(causal_NS) > 0) rbind(AF, replicate(nrow(populations), get_alternate_allele(replicate(length(causal_NS), runif(1,0.4,0.5)))))}
  
  AF = setNames(data.frame(AF), paste0("P", 1:nrow(populations), "_AF"))
  
  SNP_names = {
    if(neutral_S>0) SNP_names = sapply(1:neutral_S, function(snp_num)c(paste0('Snp_S_',snp_num,"_R0")))
    if(neutral_NS >0) SNP_names = c(SNP_names, sapply(1:neutral_NS, function(snp_num) paste0('Snp_NS_',snp_num,"_R0")))
    if(length(causal_S)>0) SNP_names = c(SNP_names, sapply(1:length(causal_S), function(snp_num)c(paste0('Snp_S_',snp_num,"_R", causal_S[snp_num]))))
    if(length(causal_NS)>0) c(SNP_names, sapply(1:length(causal_NS), function(snp_num)c(paste0('Snp_NS_',snp_num,"_R", causal_NS[snp_num]))))}
  
  data.frame(AF, `R`=c(rep(0,neutral_S + neutral_NS), causal_S, causal_NS), row.names = SNP_names)}
##TODO: parameter with numer to draw from beta distribution
get_alternate_allele <- function(reference_alleles) {
  sapply(reference_alleles, function(reference_allele) {
    s1 = reference_allele*(1-fcoeff)/fcoeff
    s2 = (1-reference_allele)*(1-fcoeff)/fcoeff
    rbeta(n = 1, shape1 = s1,shape2 = s2)})}

generate_SNPs <- function(SNP_params, populations) {
  if(trace == TRUE) print(paste(Sys.time(),"Generating SNPs dsitribution", sep=" : "))
  data = data = apply(SNP_params, 1, function(SNP_param){
    p = get_SNP_probability(SNP_param[1:length(populations)], SNP_param["R"], CC=TRUE)
    unlist(lapply(1:nrow(populations), function(population) {
      cases=sample(c(0,1,2), size = populations[population,"Case"], prob = p[1:3, population], replace = TRUE)
      controls=sample(c(0,1,2), size = populations[population,"Control"], prob = p[4:6, population], replace = TRUE)
      c(cases, controls)}))})
  rownames(data) <- get_samples_name(populations)
  data}

get_samples_name <- function(populations) {
  #Refactor this by vectorizing
  unlist(lapply(1:nrow(populations), function(population) {
    c(if(populations[population,"Case"] > 0) {paste('Case', 'P', population, 'sample', 1:populations[population,"Case"],sep='_')},
      if(populations[population,"Control"] > 0) {paste('Control', 'P', population, 'sample', 1:populations[population,"Control"],sep='_')})}))}

get_SNP_probability <- function(AFs, R = NA, CC = FALSE) {
  as.data.frame(lapply(AFs, function(alternate_allele) {
    if(CC == FALSE) {
      c((1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2)}
    else {
      c(`Case` =
          if(!is.na(R)) {
            p_associated = c((1-alternate_allele)^2, 2*as.numeric(R)*alternate_allele*(1-alternate_allele), as.numeric(R)^2 * alternate_allele^2 )
            p_associated/(sum(p_associated) + (1-alternate_allele^2))}
        else { 
          c( (1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2 )},
        `Control` = c( (1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2))}}))}

analyse_SNP <- function(SNPs, populations) {
  simulated_PC = rep(1:nrow(populations), populations[1:nrow(populations),"Control"] + populations[1:nrow(populations), "Case"])
  simulated_PC = chartr("123456789", "ABCDEFGHI", simulated_PC)
  y = unlist(lapply(1:nrow(populations), function(nb_pop) c( rep(1, populations[nb_pop,"Case"]), rep(0, populations[nb_pop,"Control"]))))
  if(trace == TRUE) print(paste(Sys.time(),"Computing PCA", sep=" : "))
  ca = prcomp(SNPs, scale. = FALSE) 
  if(trace == TRUE) print(paste(Sys.time(),"Computing logistic regression with computed PC", sep=" : "))
  W_computed_PC = apply(SNPs,2, function(SNP) coef(summary(glm(y~SNP+ca$x[,1:5])))[,4][2])
  if(trace == TRUE) print(paste(Sys.time(),"Computing logistic regression with simulated  PC", sep=" : "))
  W_simulated_PC = apply(SNPs,2, function(SNP) coef(summary(glm(y~SNP+simulated_PC)))[,4][2])
  if(trace == TRUE) print(paste(Sys.time(),"Computing logistic regression without PC", sep=" : "))
  WO_PC = apply(SNPs,2, function(SNP) coef(summary(glm(y~SNP)))[,4][2])
  parse_pvalues(data.frame(WO_PC, W_simulated_PC, W_computed_PC, row.names = colnames(SNPs)), threshold)}

parse_pvalues<- function(pvalues, threshold) {
  apply(pvalues, 2 , function(pval) {
    signiff = c()
    signiff[threshold < pval] <- 1
    signiff[threshold*0.1 < pval & pval <= threshold] <- 2
    signiff[threshold*0.01 < pval & pval <= threshold*0.1] <- 3
    signiff[threshold*0.001 < pval & pval <= threshold*0.01] <- 4
    signiff[pval <= threshold*0.001] <- 5
    data.frame(pval, `signiff` = factor(signiff, level=1:5, label=c("ns","*","**","***","+")), `pval_diff` = -log10(pvalues$WO_PC/pval) )})}

regen_signiff <- function(pvalues, threshold) {
  parse_pvalues(data.frame(`WO_PC` = pvalues[,"WO_PC.pval"], `W_simulated_PC` = pvalues[,"W_simulated_PC.pval"], `W_computed_PC` = pvalues[,"W_computed_PC.pval"], row.names = rownames(pvalues)), threshold)}

summary_sim <- function(pvalues, SNP_params) {
  mapply(function(result, name) {
    FP_neutral = rownames(SNP_params[result$signiff != "ns" & !SNP_params[,"Causal"],])
    FN_causal = rownames(SNP_params[result$signiff == "ns" & SNP_params[,"Causal"],])
    z = qnorm(result[SNP_params[,"Causal"] == FALSE,]$pval/2)
    lambda = median(z^2)/0.456
    names = c("FP_neutral", "FN_causal")
    sum = c(length(FP_neutral), length(FN_causal))
    ratio = (sum/c(sum(!SNP_params[,"Causal"]), sum(SNP_params[,"Causal"]) ) )*100
    ratio[is.nan(ratio)] = 0
    total = data.frame(sum, ratio, row.names = names)
    power_gain = median(-log10(pvalues$WO_PC$pval/result$pval), na.rm = T)
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

trace_boxplot <- function(res, save = TRUE, out = TRUE, file = paste0(Sys.time(),".png")) {
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
    guides(fill=FALSE)  
  if(save == TRUE) ggsave(filename = file, width = 14, height = 18)}

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
  write.table(x = res$params, file = paste0(output_dir, tag, "-params-(", end, ").csv"))
  write.table(x = res$scenario, file = paste0(output_dir, tag, "-scenario-(", end, ").csv"))
  #cat(toJSON(res$summary_sim), file = paste0(output_dir, tag, "-summary-(", end, ").json"))
  write.table(x =res$summary_sim, file = paste0(output_dir, tag, "-summary-(", end, ").csv"))}
