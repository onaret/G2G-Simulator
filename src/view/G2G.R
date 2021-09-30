###G2G single setup method
write_G2G_setup <- function(res, tag = NULL, output_dir = paste0(getwd(), "/")) {
  dir.create(paste(output_dir), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  #plot_G2G_setup(res, out = TRUE, save = TRUE, file = paste0(output_dir, tag, "-boxplots_out"))
  plot_G2G_setup(res, out = FALSE, save = TRUE, title = tag, file = paste0(output_dir, tag, "-boxplots"))
  #write.csv(x = summary(res$pvalues_shor),file = paste0(output_dir, tag, "-summary", end, ".csv"))
  write.table(x = res$pvalues_short, file = paste0(output_dir, tag, "-pvalues", ".csv"))}

#@G2SR : G2G test results
plot_G2G_setup <- function(G2SR, save = FALSE, out = TRUE, lim = 15, title ="", file = Sys.time()) {
  if(trace == TRUE) ifelse((save == FALSE), print(paste(Sys.time(),"Ploting...", sep=" : ")), print(paste(Sys.time(),"Writting plots", sep=" : "))) 
  G2GSR = as.data.frame(-log10(do.call(cbind, lapply(G2SR$pvalues, function(ele) ele$pval))))
  G2GSRD = as.data.frame(do.call(cbind, lapply(G2SR$pvalues, function(ele) ele$pval_diff)))
  t=mapply(function(res, name) {
    colnames(res) = gsub("WO_correction.pval", "Without correction", colnames(res))
    colnames(res) = gsub("W_human_groups.pval", "With human groups", colnames(res))
    colnames(res) = gsub("W_human_PCs.pval", "With human PCs", colnames(res))
    colnames(res) = gsub("W_strain_groups.pval", "With strain groups", colnames(res))
    colnames(res) = gsub("W_both.pval", "With both human and strain covariates", colnames(res))
    res = res %>% gather(factor_key = T)
    p <- ggplot(res, aes(key, value))
    p <- if(out == TRUE) {
      # compute lower and upper whiskers
      ylim1 = boxplot.stats(res$value)$stats[c(1, 5)]
      # scale y limits based on ylim1
      p + coord_cartesian(ylim = ylim1*1.05)
    } else p
    p + geom_boxplot(outlier.color = "grey", aes(fill=key)) + #geom_hline(yintercept = -log10(threshold), colour="red") + 
      labs(title = paste(name, title," with different covariates"), y = "-log10(pval)", x="covariate") + 
      theme(axis.text.x=element_blank())
    #if(save == TRUE) ggsave(filename = paste0(file,name,".png"), width = 14, height = 18)
    p <- ggplot(res, aes(key, value))
    p + geom_boxplot(outlier.color = "grey", aes(fill=key)) +
      labs(title = paste(name, title, " with different covariates"), y = "-log10(pval)", x="covariate") + 
      coord_cartesian(ylim = c(0, lim)) +
      theme(axis.text.x=element_blank())
    if(save == TRUE) ggsave(filename = paste0(file, name, "lim-",lim,".png"), dpi = 300)
  }, list(G2GSR, G2GSRD), c("_res", "_difference_from_correction"))}


###G2G full scenario methods
plot_collapsed_G2G <- function(res, AA.scenarios, analyse, file_tag="") {
  invisible(lapply(names(res), function(analyse) {
    if(analyse == "logistic") {
      res = res[[analyse]]
      threshold = 0.05/((ncol(res)-3)*(nrow(res)/length(levels(res$Correction))))
      res = as.data.frame(res) %>% gather(AA, pvalue,-SNP, -Correction, -SNP_Tag, convert = T)
      res$pvalue = -log10(as.numeric(res$pvalue))
      
      association_table = get_association_AA_SNP(AA.scenarios)
      
      res = right_join(as.data.frame(association_table), res, by = c("SNP", "AA"))
      res$associated[is.na(res$associated)] <- F
      
      invisible(lapply(levels(res$Correction), function(correction) {
        colnames(res) <- gsub("SNP_Tag", "SNP_Gene", colnames(res))
        res = filter(res, Correction==correction)
        res = arrange(res,associated)
        p <- ggplot(res, aes(SNP, pvalue, shape=associated))
        p + geom_point(aes(color=SNP_Gene)) + geom_hline(yintercept = -log10(threshold), color = "red") + labs(title = "G2G Uncollapsed Logistic regression results", x = "SNP")
        ggsave(filename = paste0(getwd(), "/gen-data/",file_tag,analyse, "-",correction,".png"), dpi = 300)}))}
    
    else if(analyse == "skat-L" | analyse == "skato-L"|analyse == "skato-LW"| analyse == "skat-LW" | analyse == "gt") {
      res = res[[analyse]]
      threshold = 0.05/((ncol(res)-2)*(nrow(res)/length(levels(res$Correction))))
      res = res %>% gather(AA, pvalue, -SNP_Tag, -Correction, convert = T)
      res$pvalue = -log10(res$pvalue)
      association_table = get_association_AA_SNP_tag(AA.scenarios, res)
      association_AA_tag = do.call(rbind, lapply(1:nrow(AA.scenarios), function(num) do.call(rbind, lapply(AA.scenarios[num,]$id, function(aa_id) data_frame("AA" = aa_id, "AA_Tag" = factor(AA.scenarios[num,]$bio_tag, levels = AA.scenarios$bio_tag) )))))
      res = right_join(association_table, res, by = c("SNP_Tag", "AA"))
      res = right_join(association_AA_tag, res, by =c("AA"))
      res$associated[is.na(res$associated)] <- F
      
      invisible(lapply(levels(res$Correction), function(correction) {
        colnames(res) <- gsub("AA_Tag", "AA_Position", colnames(res))
        colnames(res) <- gsub("SNP_Tag", "SNP_Gene", colnames(res))
        res = filter(res, Correction==correction)
        res = arrange(res,associated)
        p <- ggplot(res, aes(SNP_Gene, pvalue, color = AA_Position, size=associated))
        p + geom_point() + geom_hline(yintercept = -log10(threshold), colour = "red") + labs(title = paste0("Collapsed on SNP in human side with ",analyse), x = "SNP Gene")
        ggsave(filename = paste0(getwd(), "/gen-data/",file_tag,analyse, "-",correction,".png"), dpi = 300)}))}
    
    else if(analyse == "G2"){
      res = res[[analyse]]
      threshold = 0.05/(ncol(res)*nrow(res))
      association_table = data_frame(`SNP_Tag` = res[,"SNP_Tag"], `AA_Tag` = factor(colnames(res[2:ncol(res)])), `associated`=T)
      res = res %>% gather(AA_Tag, pvalue, -SNP_Tag, factor_key = T)
      res$pvalue = -log10(res$pvalue)
      res = right_join(association_table, res, by = c("SNP_Tag", "AA_Tag"))
      #Do not like that
      res$AA_Tag = as.factor(res$AA_Tag)
      res$associated[is.na(res$associated)] <- F
      res = arrange(res,associated)
      colnames(res) <- gsub("AA_Tag", "AA_Position", colnames(res))
      colnames(res) <- gsub("SNP_Tag", "SNP_Gene", colnames(res))
      p <- ggplot(res, aes(SNP_Gene, pvalue, color = AA_Position, size=associated))
      p + geom_point() + geom_hline(yintercept = -log10(threshold), colour = "red") + labs(title = paste0("Collapsed on SNP in human side and AA in viral side with G2"), x = "SNP Gene")
      #     theme(axis.text = element_text(size=12, angle = 90, hjust = 1), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36)) +
      ggsave(filename = paste0(getwd(), "/gen-data/",file_tag,analyse, ".png"), dpi = 300)}}))}

get_association_AA_SNP_tag <- function(AA.scenario,res) {
  as.data.frame(
    do.call(rbind, mapply(function(AA_id, associated_SNP_tag){
      if(!is.null(unlist(associated_SNP_tag))) {
        do.call(rbind, lapply(unlist(AA_id), function(AA){
          #do.call(rbind, lapply(unlist(associated_SNP_tag), function(Tag){ data.frame(`Tag` = factor(Tag), AA, `associated` = T)}))}))}},
          do.call(rbind, lapply(unlist(associated_SNP_tag), function(SNP_Tag){ data.frame(`SNP_Tag` = factor(SNP_Tag, levels = levels(res$SNP_Tag)), AA, `associated` = T)}))}))}},
      AA.scenario$id, AA.scenario$associated_SNP_tag, SIMPLIFY = F)))}

get_association_AA_SNP <- function(AA.scenario) {
  as.data.frame(do.call(rbind, mapply(function(AA, SNP){
    if(!is.null(unlist(SNP))) {
      do.call(rbind, lapply(unlist(AA), function(AA){
        do.call(rbind, lapply(unlist(SNP), function(SNP){data.frame(SNP, AA, `associated` = T)}))}))}},
    AA.scenario$id, AA.scenario$associated_SNPs, SIMPLIFY = F)))}

plot_G2G_on_tags <- function(res, AA.scenarios, SNP.scenarios) {
  threshold = 0.05/(ncol(res)*nrow(res))
  
  res=res$logistic
  
  res_tidy = res %>% gather(AA, pvalue, -SNP, -Correction, factor_key = T)
  res_tidy = do.call(rbind, mapply(function(tag,id) {
    data.frame(SNP.tag = as.factor(tag), filter(res_tidy, SNP %in% unlist(id)))}
    , SNP.scenarios$id_tag, SNP.scenarios$id, SIMPLIFY = F))
  res_tidy = do.call(rbind, mapply(function(tag,id) {
    data.frame(AA.tag = as.factor(tag), filter(res_tidy, AA %in% unlist(id)))}
    , AA.scenarios$id_tag, AA.scenarios$id, SIMPLIFY = F))
  res_tidy$pvalue = -log10(as.numeric(res_tidy$pvalue))
  
  plot_with_correction <- function(...,save = F) {
    correction=c(...)
    p <- ggplot(filter(res_tidy, Correction %in% correction), aes(AA.tag, pvalue, colour=SNP.tag, fill=Correction))
    p + geom_boxplot(outlier.color = "grey") + 
      geom_hline(yintercept = -log10(threshold), colour="red") + 
      labs(title = paste("pvalue for full model, with corrections:", paste(correction, collapse = " ")), y = "-log10(pval)", x="covariate") + 
      scale_y_continuous(limits = c(0, 15)) +
      geom_point(aes(y = pval_ass))
    if(save) ggsave(filename = paste0(getwd(),paste0(correction, collapse = "_"), Sys.time(), ".png"))}
  
  plot_with_correction_no_ass <- function(save = F) {
    p <- ggplot(res_tidy, aes(AA.tag, colour=SNP.tag, pvalue, fill=Correction))
    p + geom_boxplot(outlier.color = "grey") + 
      geom_hline(yintercept = -log10(threshold), colour="red") + 
      labs(title = paste("pvalue for full model, with corrections:", paste("", collapse = " ")), y = "-log10(pval)", x="Scenario") + 
      scale_y_continuous(limits = c(0, 15))
    if(save) ggsave(filename = paste0(getwd(),paste0(collapse = "_"), Sys.time(), ".png"))}}

plot_pvalue_by_methods <- function(res, AA.scenarios) {
  analyses = names(res)
  res_all = lapply(analyses, function(analyse){ 
    res = res[[analyse]]
    if(analyse == "logistic") {
      res = as.data.frame(res) %>% gather(AA, pvalue, -SNP, -Correction, -SNP_Tag, convert = T)
      res$pvalue = -log10(as.numeric(res$pvalue))
      association_table = get_association_AA_SNP(AA.scenarios)
      res = inner_join(as.data.frame(association_table), res, by = c("SNP", "AA"))
      res$framework <- factor("logistic", level = analyses)
      select(res, pvalue, framework)}
    else if(analyse %in% c("skat-L", "skato-L", "gt")){
      res = res %>% gather(AA, pvalue, -SNP_Tag, -Correction, convert = T)
      res$pvalue = -log10(res$pvalue)
      association_table = get_association_AA_SNP_tag(AA.scenarios, res)
      res = inner_join(association_table, res, by = c("SNP_Tag", "AA"))
      res$framework <- factor(analyse, level = analyses)
      select(res, pvalue, framework)}
    else if(analyse == "G2"){
      ##This is were the selection occurs, without going through association table
      res = data.frame("pvalue" = -log10(diag(as.matrix(res[2:ncol(res)]))))
      res$framework <- factor("G2", level = analyses)
      res}})
  
  thresholds = get_thresholds(res)
  res_all = do.call(rbind, res_all)
  levels(res_all$framework) <- paste("Method ", c("Logistic Regression", "Global Test", "SKAT", "SKATO", "G2" ))
  
  p <- ggplot(res_all,mapping = aes(framework,pvalue))
  p + geom_boxplot(outlier.color = "grey", aes(fill=framework)) + 
    geom_hline(data = thresholds, aes(yintercept = pvalue, colour = collapsing)) + 
    labs(title = "pvalue distribution for association in function of collapsing or framework", x = "Framework") + 
    theme(axis.text = element_blank())
  ggsave(filename = paste0(getwd(), "/gen-data/","framework-pvalue-dist",".png"), dpi = 300)}

###Optional
plot_median_pvalue_trend <- function(res) {
  ###To implement from prototype script to framework
  res_all_med = do.call(rbind, lapply(1:length(res_test), function(index){
    res_collapsed = lapply(c("skat-L", "skato-L", "gt"), function(analyse) {
      res = res_test[[index]]$res[[analyse]]
      res = res %>% gather(AA, pvalue, -Tag, -Correction, convert = T)
      res$pvalue = -log10(res$pvalue)
      association_table = get_association_AA_SNP_tag(AA.scenarios[[index]], res)
      res_joined = inner_join(association_table, res, by = c("Tag", "AA"))
      ass_num = res_test[[index]]$associated_SNP_nb
      pvalue <- median(res_joined$pvalue)
      data.frame("associated_num" = ass_num, pvalue, "framework" = factor(analyse, level = framework_level))})
    
    res_G2 = {
      res = as.data.frame(res_test[[index]]$res[["G2"]])
      ass_num = res_test[[index]]$associated_SNP_nb
      pvalue <- -log10(median(diag(as.matrix(res[2:ncol(res)]))))
      data.frame("associated_num" = ass_num, pvalue, "framework" = factor("G2", level = framework_level))}
    
    res_logistic = {
      res = res_test[[index]]$res[["logistic"]]
      res = as.data.frame(res) %>% gather(AA, pvalue, -SNP, -Correction, -Tag, convert = T)
      res$pvalue = -log10(as.numeric(res$pvalue))
      association_table = get_association_AA_SNP(AA.scenarios[[index]])
      res_joined = inner_join(as.data.frame(association_table), res, by = c("SNP", "AA"))
      res = data.frame()
      ass_num = res_test[[index]]$associated_SNP_nb
      pvalue <- median(res_joined$pvalue)
      data.frame("associated_num" = ass_num, pvalue, "framework" = factor("logistic", level = framework_level))}
    
    rbind(res_logistic, res_G2, do.call(rbind, res_collapsed))}))
  
  res_all_med = as.data.frame(res_all_med) 
  levels(res_all_med$framework) <- paste("Method ", c("Logistic Regression", "Global Test", "SKAT", "SKATO", "G2" ))
  
  thresholds = get_thresholds(res_test[[1]])
  p <- ggplot(res_all_med, mapping = aes(associated_num, pvalue))
  p + geom_smooth(aes(colour = framework)) + 
    geom_hline(data = thresholds, aes(yintercept = pvalue, colour = collapsing)) + 
    labs(title = "pvalue trend in function of associated SNP number, for the different collapsing or framework", x = "associated SNP") + 
    theme(axis.text = element_text(size=12, angle = 90, hjust = 1), axis.title=element_text(size=18,face="bold"), plot.title = element_text(size = 18))
  ggsave(filename = paste0(getwd(), "/gen-data/","framework-trend",".png"), dpi = 300)}

get_thresholds <- function(res) {
  threshold_seter = res[["gt"]]
  threshold_1C = -log10(0.05/((ncol(threshold_seter)-2)*(nrow(threshold_seter)/length(levels(threshold_seter$Correction)))))
  
  threshold_seter = res[["G2"]]
  threshold_2C = -log10(0.05/((ncol(threshold_seter)-1)*(nrow(threshold_seter))))
  
  threshold_seter = res[["logistic"]]
  threshold_logistic = -log10(0.05/((ncol(threshold_seter)-3)*(nrow(threshold_seter)/length(levels(threshold_seter$Correction)))))
  gather(data_frame(`Threshold Logistic R` = threshold_logistic, `Threshold Collapsing SNPs` = threshold_1C, `Threshold Collapsing SNPs and AAPV` =threshold_2C), collapsing, pvalue)}
