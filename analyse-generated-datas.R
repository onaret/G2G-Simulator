prefix = "Analysis/pop-structure/"

#prefix = "gen-data/"

prefix ="param-analysis/"

load_datas <- function(pattern) {
  sapply(list.files(prefix, pattern = pattern), function(file) {
    setNames(list(read.table(file = paste0(prefix,file))), file)})}

load_populations <- function() load_datas("populations")
load_SNP_struc <- function() load_datas("SNP_struc")
load_pvalues <- function() load_datas("pvals")
load_summary <- function() load_datas("summary")
load_SNP_freq <- function() load_datas("SNP_freq")

SNP_res = data.frame(unname(load_datas("11N12C4-SNPs_freq.csv")))
SNP_res_pval = data.frame(unname(load_datas("11N12C4-pvalues.csv")))
meltsnp = melt(SNP_res_pval,id.vars = c(1,2,3,4,5,6))
meltsnp = melt(SNP_res_pval,id.vars = c(1,2))

#colnames(SNP_res) <- c("P1", "P2", "P3", "R")


SNP_rem_WO = (SNP_res_pval[,1:2])
SNP_rem_WO[,1] = -log10(SNP_rem_WO[,1, drop = F])

SNP_rem_sim = (SNP_res_pval[,3:4])
SNP_rem_sim[,1] = -log10(SNP_rem_sim[,1, drop = F])

SNP_rem_comp = (SNP_res_pval[,5:6])
SNP_rem_comp[,1] = -log10(SNP_rem_comp[,1, drop = F])

#colnames(SNP_rem_WO) <- c("pvalue", "signiff")
#colnames(SNP_rem_sim) <- c("pvalue", "signiff")
#colnames(SNP_rem_comp) <- c("pvalue", "signiff")

colnames(SNP_rem_WO) <- c("pvalue")
colnames(SNP_rem_sim) <- c("pvalue")
colnames(SNP_rem_comp) <- c("pvalue")

#SNP_rem_WO$correction <- "Without correction"
#SNP_rem_sim$correction <- "Without human groups"
#SNP_rem_comp$correction <- "Without computed PC"

SNP_rem_WO$SNP_num <- 1:nrow(SNP_rem_WO)
SNP_rem_sim$SNP_num <- 1:nrow(SNP_rem_WO)
SNP_rem_comp$SNP_num <- 1:nrow(SNP_rem_WO)

SNP_rem_WO$SNP_type <- c(rep("Stratified", 5000), rep("Non-Stratified", 95000), rep("Causal Stratified", 21), rep("Causal Non-Stratified", 21))
SNP_rem_sim$SNP_type <- c(rep("Stratified", 5000), rep("Non-Stratified", 95000), rep("Causal Stratified", 21), rep("Causal Non-Stratified", 21))
SNP_rem_comp$SNP_type <- c(rep("Stratified", 5000), rep("Non-Stratified", 95000), rep("Causal Stratified", 21), rep("Causal Non-Stratified", 21))

threshold = 5e-7
main_name = "Without correction"
p <- ggplot(SNP_rem_WO, aes(SNP_num, pvalue))
p + geom_point(aes(colour = SNP_type)) + scale_colour_manual(values =c("orange", "red", "black", "yellow")) + geom_hline(yintercept = -log10(threshold)) + labs(title = main_name, x = "SNP")+
  theme(axis.text = element_text(size=24), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36)) +
  scale_y_continuous(limits = c(0, 45))

main_name = "With human groups"
p <- ggplot(SNP_rem_sim, aes(SNP_num, pvalue))
p + geom_point(aes(colour = SNP_type)) + scale_colour_manual(values =c("orange", "red", "black", "yellow")) + geom_hline(yintercept = -log10(threshold)) + labs(title = main_name, x = "SNP") +
  theme(axis.text = element_text(size=24), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36)) +
  scale_y_continuous(limits = c(0, 45))

main_name = "With computed PC"
p <- ggplot(SNP_rem_comp, aes(SNP_num, pvalue))
p + geom_point(aes(colour = SNP_type)) + scale_colour_manual(values =c("orange", "red", "black", "yellow")) + geom_hline(yintercept = -log10(threshold)) + labs(title = main_name, x = "SNP")+
  theme(axis.text = element_text(size=24), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36)) +
  scale_y_continuous(limits = c(0, 45))

SNP_rem_WO$pvector= sort(SNP_rem_WO$pvalue,decreasing=F)
SNP_rem_sim$pvector = sort(SNP_rem_sim$pvalue,decreasing=F)
SNP_rem_comp$pvector = sort(SNP_rem_comp$pvalue,decreasing=F)

SNP_rem_WO$y = sort(-log10( 1:length(SNP_rem_WO$pvector)/length(SNP_rem_WO$pvector) ),decreasing=F)
SNP_rem_sim$y = sort(-log10( 1:length(SNP_rem_WO$pvector)/length(SNP_rem_WO$pvector) ),decreasing=F)
SNP_rem_comp$y = sort(-log10( 1:length(SNP_rem_WO$pvector)/length(SNP_rem_WO$pvector) ),decreasing=F)

p <- ggplot(SNP_rem_WO, aes(y,pvector))
p + geom_point() + labs(title="QQ-plot: Without Correction", x = "Expected  -Log(pval)", y = "Observed -Log(pvalue)") + geom_abline(intercept = 0, colour ="red") + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_text(size=24), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36))+
  scale_y_continuous(limits = c(0, 40))

p <- ggplot(SNP_rem_sim, aes(y, pvector))
p + geom_point() + labs(title="QQ-plot: With human groups", x = "Expected  -Log(pval)", y = "Observed -Log(pvalue)") + geom_abline(intercept = 0, colour ="red") + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  theme(axis.text = element_text(size=24), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36))+
  scale_y_continuous(limits = c(0, 40))

p <- ggplot(SNP_rem_comp, aes(y, pvector))
p + geom_point() + labs(title="QQ-plot: With computed PC", x = "Expected  -Log(pval)", y = "Observed -Log(pvalue)") + geom_abline(intercept = 0, colour ="red") + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  theme(axis.text = element_text(size=24), axis.title=element_text(size=32,face="bold"), plot.title = element_text(size = 36))+
  scale_y_continuous(limits = c(0, 40))









  
t= load_datas("result_on_scenario_fs1_bias")
t = t(data.frame(read.table(file = paste0(prefix,"result_on_scenario_fs1_bias"))))
t = unname(t(t[[1]]))[1:40,]
x <- unname(t[,1])
y <- unname(-log10(t[,2]))
qplot(x,y, geom='smooth', span =0.5, main="Fst on pvalue", xlab = "Fst", ylab = "-Log10(pval)", margins = TRUE)

t = data.frame(t(tb))
colnames(t) <- c("Fst", "Without correction", "With human groups", "With viral groups", "With human and viral groups")
meltt = melt(t,id.vars = c(1))
meltt$value = -log10(meltt$value)
colnames(meltt)[2] <- "Correction"
f <- ggplot(meltt, aes(Fst, value))
f + geom_smooth(aes(colour = Correction)) + scale_x_continuous(limits = c(0, 0.4)) + 
  labs(title = "Fst parameter on median pvalue", y = "Median( -Log10(pvalue) )", x="Fst")


tb = read.table(file = paste0(prefix,"result_on_scenario_ss1_beta"))
t = data.frame(t(tb))
colnames(t) <- c("beta", "Without correction", "With human groups", "With viral groups", "With human and viral groups")
meltt = melt(t,id.vars = c(1))
meltt$value = -log10(meltt$value)
colnames(meltt)[2] <- "Correction"
f <- ggplot(meltt, aes(beta, value))
f + geom_smooth(aes(colour = Correction)) + #scale_x_continuous(limits = c(0, 0.4)) + 
  labs(title = "beta parameter on median pvalue", y = "Median( -Log10(pvalue) )", x="beta")

tbeta = fs2_beta
t = unname(t(t[[1]]))[1:40,]
x <- unname(t[,1])
y <- unname(-log10(t[,2]))
qplot(x,y, geom='smooth', span =0.5, main="Fst on pvalue", xlab = "Fst", ylab = "-Log10(pval)", margins = TRUE)

tb = read.table(file = paste0(prefix,"ss4_size"))
t = data.frame(t(tb))
rownames(t) <- NULL
colnames(t) <- c("size", "Without correction", "With human groups", "With viral groups", "With human and viral groups")
meltt = melt(t,id.vars = c(1))
meltt$value = -log10(meltt$value)
colnames(meltt)[2] <- "Correction"
f <- ggplot(meltt, aes(size, value))
f + geom_smooth(aes(colour = Correction)) + #scale_x_continuous(limits = c(0, 0.4)) + 
  labs(title = "size parameter on median pvalue", y = "Median( -Log10(pvalue) )", x="beta")

test <- function(x) {1/(1+exp(-x))}
f <- ggplot(data.frame(x = c(-12, 12)), aes(x))
f + stat_function(fun = test, colour = "red") +
  geom_hline(yintercept = 0, colour="blue", linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  geom_hline(yintercept = 1, colour="blue", linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = test(0.4*2), linetype = "dotdash", colour = "green") +
  geom_hline(yintercept = test(-0.4*2), linetype = "dotdash", colour = "green") +
  geom_vline(xintercept = 0.4*2, linetype = "dotdash", colour = "green") +
  geom_vline(xintercept = -0.4*2, linetype = "dotdash", colour = "green") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), plot.title = element_text(size = 26)) 


newfun <- function(x,k=1,y=0.8) {1/(1+exp(-k*(x-solution(y))))}
f <- ggplot(data.frame(x = c(-12, 12), aes(x)))
f + stat_function(fun = newfun, colour = "red") +
  geom_hline(yintercept = 0, colour="blue", linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  geom_hline(yintercept = 1, colour="blue", linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), plot.title = element_text(size = 26)) 


test2 <- function(x,y) {1/(1+exp(-x))}
f <- ggplot(data.frame(x = c(-12, 12)), aes(x))
f + stat_function(fun = newfun, colour = "red") +
  geom_hline(yintercept = 0, colour="blue", linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  geom_hline(yintercept = 1, colour="blue", linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), plot.title = element_text(size = 26)) 


problem <- function(x) {1/(1+exp(x))}
f <- ggplot(data.frame(x = c(-12, 12)), aes(x))
f + stat_function(fun = problem, colour = "red") +
  geom_hline(yintercept = 0, colour="blue", linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  geom_hline(yintercept = 1, colour="blue", linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), plot.title = element_text(size = 26)) 

solution <- function(x){log((1/x) - 1)}
f <- ggplot(data.frame(x = c(0, 1)), aes(x))
f + stat_function(fun = solution, colour = "red") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), plot.title = element_text(size = 26)) 



test2 <- function(x,y) {-x*log(x)}
f <- ggplot(data.frame(x = c(0, 1)), aes(x))
f + stat_function(fun = test2, colour = "red") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), plot.title = element_text(size = 26)) 

test3 <- function(x,y) {log(x)}
f <- ggplot(data.frame(x = c(0, 1)), aes(x))
f + stat_function(fun = test2, colour = "red") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), plot.title = element_text(size = 26)) 


derivatelogit <- function(x,y) {exp(x)/((exp(x) +1)^2)}
f <- ggplot(data.frame(x = c(-12, 12)), aes(x))
f + stat_function(fun = derivatelogit, colour = "red") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), plot.title = element_text(size = 26)) 


derivatelogitparm <- function(x,y) {exp(x+1/25)/(exp(x)+exp(1/25))^2}
f <- ggplot(data.frame(x = c(-12, 12)), aes(x))
f + stat_function(fun = derivatelogitparm, colour = "red") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), plot.title = element_text(size = 26)) 





f <- ggplot(data.frame(x = c(0, 10)), aes(x))
f + stat_function(fun = sin, colour = "red") +  stat_function(fun = cos, colour = "blue") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"), plot.title = element_text(size = 26)) + guides(fill=FALSE)  
