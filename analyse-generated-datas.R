prefix = "gen-data/"

load_datas <- function(pattern) {
  sapply(list.files(prefix, pattern = pattern), function(file) {
    setNames(list(read.table(file = paste0(prefix,file))), file)})}

load_populations <- function() load_datas("populations")
load_SNP_struc <- function() load_datas("SNP_struc")
load_pvalues <- function() load_datas("pvalues")
load_summary <- function() load_datas("summary")
load_SNP_freq <- function() load_datas("SNP_freq")