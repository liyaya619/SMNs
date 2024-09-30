source("step1_lib.R")

load("phegwas_result_list.Rdata")

for (i in 1:length(result_list)) {
  exposure <- result_list[[i]]
  outcome <- read_outcome_data("./second_malignant_neoplasm/dat_38/ GCST90079616 .tsv.gz", 
                               sep = "\t", 
                               snps = exposure$SNP, 
                               snp_col = "SNP",
                               effect_allele_col = "A2", 
                               other_allele_col = "A1",
                               pval_col = "P",
                               beta_col = "Effect", 
                               se_col = "SE", 
                               eaf_col = "FRQ", 
                               chr_col = "CHR")
}
