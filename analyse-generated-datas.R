prefix = "gen-data/"

load_datas <- function(pattern) {
  lapply(list.files(prefix, pattern = pattern), function(file) {
    read.csv2(file = paste0(prefix,file))})}

load_populations <- function() load_datas("populations")
load_SNP_struc <- function() load_datas("SNP_struc")
load_pvalues <- function() load_datas("pvalues")
load_summary <- function() load_datas("summary")
load_SNP_freq <- function() load_datas("SNP_freq")