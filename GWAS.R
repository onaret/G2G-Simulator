####Function
multiple_GWAS_scenarios <- function(populations, neutral, neutral_S_rat, causal_S=c(), causal_NS=c(), repetitions) {
  #This condition is ugly
  cl <- makeCluster(CPU, outfile = "output.txt", type = "FORK")
  if(is.list(populations) && is.list(populations[[1]])) {
    parLapply(cl, 1:(repetitions*length(populations)), function(repetition, populations, names){
      pop_index = repetition%%length(populations) + 1
      print(paste(Sys.time(),": For population", names[pop_index], "executing repetition", trunc(repetition/length(populations))+1))
      write(GWAS_scenario(populations[pop_index], neutral = neutral, neutral_S_rate, causal_S, causal_NS), repetition, names[pop_index])
    },populations, names(populations) )}
  else {
    lapply(1:repetitions, function(repetition){
      if(trace == TRUE) print(paste(Sys.time(),": Executing repetition", repetition))
      write(GWAS_scenario(populations, neutral, neutral_S_rate, causal_S, causal_NS), repetition)})}
  if(trace == TRUE) print(paste(Sys.time(),"Done", sep=" : "))
  stopCluster(cl)}

GWAS_scenario <- function(populations, neutral=0, neutral_S_rate=0, causal_S=c(), causal_NS=c(), fst_strat, nb_pc=5, SNPs=NULL) {
  SNPs = if(is.null(SNPs)) generate_SNPs_for_GWAS(populations, neutral, neutral_S_rate, causal_S, causal_NS, fst_strat) else SNPs
  res = analyse_GWAS(SNPs$data, populations, nb_pc)
  threshold =  0.05/nrow(SNPs$params)
  pvalues = parse_pvalues(res, threshold)
  summary_sim = summary_sim(pvalues, SNPs$params, "Causal")
  populations = rbind(populations, matrix(c(sum(populations[,"Case"]), sum(populations[,"Control"])),ncol = 2, nrow = 1, dimnames = list("Total") ) )
  list(#`SNPs` = SNPs,
    `study_design` = populations,
    `SNP_params` = SNPs$params,
    `pvalues` = pvalues,
    `summary_sim` = summary_sim,
    `params` = list(`neutral_S_rate` = neutral_S_rate, `fst_strat` = fst_strat, `threshold` = threshold),
    `populations` = populations)}

generate_SNPs_for_GWAS <- function(populations, neutral=0, neutral_S_rate=0, causal_S=c(), causal_NS=c(), fst_strat) {
  neutral_S = neutral_S_rate * neutral
  SNP_params = get_SNP_frequencies(neutral - neutral_S, neutral_S, causal_S, causal_NS, populations, fst_strat)
  if(trace == TRUE) print(paste(Sys.time(),"Generating SNPs dsitribution", sep=" : "))
  SNPs = apply(SNP_params, 1, function(SNP_param){
    p = get_genotype_probability(SNP_param[1:nrow(populations)], SNP_param["R"])
    unlist(lapply(1:nrow(populations), function(population) {
      cases=sample(c(0,1,2), size = populations[population,"Case"], prob = p[1:3, population], replace = TRUE)
      controls=sample(c(0,1,2), size = populations[population,"Control"], prob = p[4:6, population], replace = TRUE)
      c(cases, controls)}))})
  rownames(SNPs) <- unlist(lapply(1:nrow(populations), function(population) {
    c(if(populations[population,"Case"] > 0) {paste('Case', 'P', population, 'sample', 1:populations[population,"Case"],sep='_')},
      if(populations[population,"Control"] > 0) {paste('Control', 'P', population, 'sample', 1:populations[population,"Control"],sep='_')})}))
  list(`data` = SNPs, `params` = SNP_params)}

generate_population_for_GWAS <- function(populations) {
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

to_study_design_for_GWAS<- function(populations) {
  do.call(rbind, lapply(1:nrow(populations), function(pop_num) {rbind(data.frame(CC = rep("Case", populations[pop_num,1]), Population = paste0("P",pop_num) ),data.frame(CC = rep("Control", populations[pop_num,2]), Population = paste0("P",pop_num) )) }))}

get_SNP_frequencies <- function(neutral_NS, neutral_S, causal_S, causal_NS, populations, fst_strat) {
  if(trace == TRUE) print(paste(Sys.time(),"Generating allele frequency", sep=" : "))
  # nb_pop = length(levels(populations$Population))
  nb_pop = ncol(populations)
  f_neutral_NS = if(neutral_NS >0) data.frame(t(data.frame(replicate(neutral_NS, rep(get_AF(allele = runif(1,0.4,0.5), fst = fst_strat), nb_pop),simplify = F))), Stratified = F, Causal = F, R= NA)
  f_neutral_S = if(neutral_S>0) data.frame(t(data.frame(replicate(neutral_S, get_AF(allele = runif(1,0.4,0.5), fst = fst_strat, nb = nb_pop),simplify = F))), Stratified = T, Causal = F, R= NA)
  f_causal_NS = if(length(causal_NS)>0) data.frame(t(data.frame(replicate(length(causal_NS), rep(get_AF(allele = runif(1,0.4,0.5), fst = fst_strat), nb_pop),simplify = F))), Stratified = F, Causal = T, R= causal_NS)
  f_causal_S = if(length(causal_S)>0) data.frame(t(data.frame(replicate(length(causal_S), get_AF(allele = runif(1,0.4,0.5), fst = fst_strat, nb = nb_pop),simplify = F))), Stratified = T, Causal = T, R= causal_S)
  SNP_params = rbind(f_neutral_NS, f_neutral_S, f_causal_NS, f_causal_S)
  colnames(SNP_params) <- c(paste0("P", 1:nb_pop), "Stratified", "Causal", "R")
  SNP_names = c(if(neutral_NS >0) paste0('Snp_NS_',1:neutral_NS,"_R0"),
                if(neutral_S>0) paste0('Snp_S_',1:neutral_S,"_R0"),
                if(length(causal_NS)>0) paste0('Snp_NS_',1:length(causal_NS),"_R_", causal_NS),
                if(length(causal_S)>0) paste0('Snp_S_',1:length(causal_S),"_R_", causal_S))
  rownames(SNP_params) <- SNP_names
  SNP_params}

get_genotype_probability <- function(AFs, R = NA) {
  as.data.frame(lapply(AFs, function(alternate_allele) {
    c(`Case` =
        if(!is.na(R)) {
          p_associated = c((1-alternate_allele)^2, 2*as.numeric(R)*alternate_allele*(1-alternate_allele), as.numeric(R)^2 * alternate_allele^2 )
          p_associated/(sum(p_associated) + (1-alternate_allele^2))}
      else { 
        c( (1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2 )},
      `Control` = c( (1-alternate_allele)^2, 2*alternate_allele*(1-alternate_allele), alternate_allele^2))}))}

analyse_GWAS <- function(SNPs, populations, nb_pc) {
  simulated_PC = rep(1:nrow(populations), populations[1:nrow(populations),"Control"] + populations[1:nrow(populations), "Case"])
  simulated_PC = chartr("123456789", "ABCDEFGHI", simulated_PC)
  y = unlist(lapply(1:nrow(populations), function(nb_pop) c( rep(1, populations[nb_pop,"Case"]), rep(0, populations[nb_pop,"Control"]))))
  if(trace == TRUE) print(paste(Sys.time(),"Computing PCA", sep=" : "))
  SNP_PC = prcomp(SNPs, scale. = FALSE)
  SNP_PC = SNP_PC$x[,1:ifelse(ncol(SNP_PC$x)<nb_pc, ncol(SNP_PC$x), nb_pc)]
  if(trace == TRUE) print(paste(Sys.time(),"Computing logistic regression with computed PC", sep=" : "))
  W_computed_PC = apply(SNPs,2, function(SNP) coef(summary(glm(y~SNP+SNP_PC)))[,4][2])
  if(trace == TRUE) print(paste(Sys.time(),"Computing logistic regression with simulated  PC", sep=" : "))
  W_simulated_PC = apply(SNPs,2, function(SNP) coef(summary(glm(y~SNP+simulated_PC)))[,4][2])
  if(trace == TRUE) print(paste(Sys.time(),"Computing logistic regression without PC", sep=" : "))
  WO_correction = apply(SNPs,2, function(SNP) coef(summary(glm(y~SNP)))[,4][2])
  data.frame(WO_correction, W_simulated_PC, W_computed_PC, row.names = colnames(SNPs))}

##@G
parse_pvalues<- function(pvalues, threshold) {
	apply(pvalues, 2 , function(pval) {
		signiff = c()
		signiff[threshold < pval] <- 1
		signiff[threshold*0.1 < pval & pval <= threshold] <- 2
		signiff[threshold*0.01 < pval & pval <= threshold*0.1] <- 3
		signiff[threshold*0.001 < pval & pval <= threshold*0.01] <- 4
		signiff[pval <= threshold*0.001] <- 5
		data.frame(pval, `signiff` = factor(signiff, level=1:5, label=c("ns","*","**","***","+")), `pval_diff` = -log10(pvalues$WO_correction/pval) )})}

##@G
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