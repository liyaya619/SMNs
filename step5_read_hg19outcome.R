source("step1_lib.R")

overall_result_list <- readRDS("./Europe_overall_phegwas_result_list.rds")

# dir.create("./mr_outputs")
result_list <- overall_result_list

for (i in 1:length(result_list)) {
  exposure <- result_list[[i]]
  ##outcome data
  dir <- "./second_malignant_neoplasm"
  out <- list.files(dir,pattern = "37.tsv.gz")
  for (j in 1:length(out)) {
    outdat <- out[j]
    print(outdat)
    
    expname <- unique(exposure$id.exposure)
    outname <- trimws(substr(outdat,1,12))
    
    outcome<-tryCatch(
      outcome <- read_outcome_data(file.path(dir,outdat), 
                                   sep = "\t", 
                                   snps = exposure$SNP, 
                                   snp_col = "variant_id",
                                   effect_allele_col = "effect_allele", 
                                   other_allele_col = "other_allele",
                                   pval_col = "p_value",
                                   samplesize_col = "N",
                                   beta_col = "beta", 
                                   se_col = "standard_error", 
                                   eaf_col = "effect_allele_frequency", 
                                   chr_col = "chromosome"),
      error = function(e) {
        message(paste("Error occurred for outcome:", outname))
        NULL  # Return NULL when an error occurs
      }
    )
    
    if (!is.null(outcome)) {
      dat<-harmonise_data(exposure,outcome)
      
      mr.methods=c("mr_wald_ratio", "mr_ivw", "mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")
      mr.res <- mr(dat, method_list=mr.methods)
      
      threshold <- 0.05 ##定义p值范围
      target_column <- "pval"
      if (any(mr.res[[target_column]] < threshold)) {
        # 如果存在小于阈值的元素，则执行相应的操作
        print(paste("At least one value in column", target_column, "is less than", threshold))
        # 在这里你可以执行其他的操作
        ors <- generate_odds_ratios(mr.res)
        write.csv(ors, paste0("./mr_outputs/",expname,"_",outname,"_mr_results.csv"))
      }
    }
    
  }
}
