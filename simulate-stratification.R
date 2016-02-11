rm(list = ls())
cat("\014")  
library(rjson)
library(parallel)
library(pryr)
setwd("/home/onaret/workspace")
#setwd("/home/zod/Documents/Workspace/EPFL/simulate-stratification")

source("generate-SNP-datas.R")
source("generate-viral-datas.R")
source("analyse-generated-datas.R")

####Constants
C1 = list(`P1` = c(`case` = 200, `control` = 400), `P2` = c(`case` = 400, `control`  = 200))
C2 = list(`P1` = c(`case` = 400, `control` = 200), `P2` = c(`case` = 200, `control` = 400))
C3 = list(`P1` = c(`case` = 300, `control` = 0), `P2` = c(`case` = 300, `control` = 600))
C4 = list(`P1` = c(`case` = 300, `control` = 200), `P2` = c(`case` = 200, `control` = 100), `P3` = c(`case` = 100, `control` = 300))
C5 = list(`P1` = c(`case` = 200, `control` = 0), `P2` = c(`case` = 400, `control` = 200), `P3` = c(`case` = 0, `control` = 400))
AllPop = list(`C1` = C1, `C2`= C2, `C3` = C3, `C4` = C4, `C5` = C5)

AllPop = list(`C1` = C1, `C2`= C2, `C3` = C3, `C4` = C4, `C5` = C5)
rat = list(`C3` = C3, `C4` = C4, `C5` = C5,`C3` = C3, `C4` = C4, `C5` = C5,`C3` = C3, `C4` = C4, `C5` = C5)

fcoeff  = 0.01 ##### Wright's coefficient for inbreeding
#threshold = 3.63e--8
threshold = 5e-07
sequence = 11
CPU = trunc(detectCores() * (1/2))
######Scenarios
run_multiple_scenarios(populations = AllPop,
                      #populations = c(`min`=2, `max`=8, `size`=1200),
                      neutral = 100000, 
                      neutral_S_rate = 0.05, 
                      causal_NS = seq(1,2, by = 0.05),
                      causal_S = seq(1,2, by = 0.05),
                      repetitions = 200)

#populations = load_populations()
#SNP_struc = load_SNP_struc()
#pvalues = load_pvalues()
#summary = load_summary()
#SNP_freq = load_SNP_freq()
#test = regen_signiff(pvalues$`10N1-pvalues-(2016-02-08 14:54:40).csv.10N1-pvalues-(2016-02-08 14:54:40).csv`, 5e-06)
