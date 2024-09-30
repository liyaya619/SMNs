rm(list = ls())

# https://htmlpreview.github.io/?https://raw.githubusercontent.com/MRCIEU/ukbb-gwas-analysis/master/docs/ldsc_clumped_analysis.html?token=AAOV6TBQXEXEPT7SUXXLWMC6DWP3O
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

library(dplyr)
library(tidyr)
library(ieugwasr)
library(TwoSampleMR)
# library(doParallel)
library(foreach)

load("available_outcomes.Rdata")

rfMRI <- subset(ao,grepl("ubm",id) & population=="European"& sample_size>1000)

# get IVs -----------------------------------------------------------------
# 定义每次循环处理的子集大小
subset_size <- 100

exp <- rfMRI$id

ubm_list <- list()  # Create an empty list to store the results


# 计算需要的循环次数
num_iterations <- ceiling(length(exp) / subset_size)

# 遍历每个子集
for (i in 1:num_iterations) {
  # i=38
  start_index <- (i - 1) * subset_size + 1  # 子集的起始索引
  end_index <- min(i * subset_size, length(exp))  # 子集的结束索引
  
  subset <- exp[start_index:end_index]  # 提取子集进行处理
  
  # 在此处进行子集处理
  for (j in subset) {
    print(j)
    exp_gwas <- tryCatch(
      TwoSampleMR::extract_instruments(outcomes = j, p1 = 5e-08, clump = TRUE,
                                       p2 = 5e-08, r2 = 0.01, kb = 10000,
                                       access_token = ieugwasr::check_access_token(),
                                       force_server = FALSE),
      error = function(e) {
        message(paste("Error occurred for GWAS:", j))
        NULL # Return NULL when an error occurs
      }
    )
    
    if (!is.null(exp_gwas)) {
      ubm_list[[j]] <- exp_gwas
    }
  }
  
}

# 保存结果列表
# print(names(ubm_list))
saveRDS(ubm_list, file = "Europe_ubm_result_list.rds")


