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

#@SNPs: if you want to add already generated SNPs
GWAS_scenario <- function(populations, neutral=0, neutral_S_rate=0, causal_S=c(), causal_NS=c(), fst_strat, nb_pc=5, SNPs=NULL) {
  SNPs = if(is.null(SNPs)) generate_SNPs_for_GWAS(populations, neutral, neutral_S_rate, causal_S, causal_NS, fst_strat) else SNPs
  res = analyse_GWAS(SNPs$data, populations, nb_pc)
  threshold =  0.05/nrow(SNPs$params)
  pvalues = parse_pvalues(res, threshold)
  summary_sim = summary_sim(pvalues, SNPs$params)
  populations = rbind(populations, matrix(c(sum(populations[,"Case"]), sum(populations[,"Control"])),ncol = 2, nrow = 1, dimnames = list("Total") ) )
  list(#`SNPs` = SNPs,
    `study_design` = populations,
    `SNP_params` = SNPs$params,
    `pvalues` = pvalues,
    `summary_sim` = summary_sim,
    `params` = list(`neutral_S_rate` = neutral_S_rate, `fst_strat` = fst_strat, `threshold` = threshold),
    `populations` = populations)}

####Function
multiple_GWAS_scenarios <- function(populations, neutral, neutral_S_rat, causal_S=c(), causal_NS=c(), repetitions) {
  #This condition is ugly
  if(is.list(populations) && is.list(populations[[1]])) {
    cl <- makeCluster(CPU, outfile = "output.txt", type = "FORK")
    parLapply(cl, 1:(repetitions*length(populations)), function(repetition, populations, names){
      pop_index = repetition%%length(populations) + 1
      print(paste(Sys.time(),": For population", names[pop_index], "executing repetition", trunc(repetition/length(populations))+1))
      write(GWAS_scenario(populations[pop_index], neutral = neutral, neutral_S_rate, causal_S, causal_NS), repetition, names[pop_index])
    },populations, names(populations) )
    stopCluster(cl)}
  else {
    lapply(1:repetitions, function(repetition){
      if(trace == TRUE) print(paste(Sys.time(),": Executing repetition", repetition))
      write(GWAS_scenario(populations, neutral, neutral_S_rate, causal_S, causal_NS), repetition)})}
  if(trace == TRUE) print(paste(Sys.time(),"Done", sep=" : "))}

to_study_design_for_GWAS <- function(populations) {
  do.call(rbind, lapply(1:nrow(populations), function(pop_num) {rbind(data.frame(CC = rep("Case", populations[pop_num,1]), Population = paste0("P",pop_num) ),data.frame(CC = rep("Control", populations[pop_num,2]), Population = paste0("P",pop_num) )) }))}
