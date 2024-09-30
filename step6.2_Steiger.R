source("step1_lib.R")

overall_result_list <- readRDS("./Europe_overall_phegwas_result_list.rds")

file_path <- "./filtered_steiger"

if (!file.exists(file_path)) {
  # 创建文件
  file.create(file_path)
} else {
  # 文件已存在，跳过
}


for (i in 1:length(result_list)) {
  # i=3
  exposure <- result_list[[i]]
  ##outcome data
  dir <- "./second_malignant_neoplasm"
  out <- list.files(dir,pattern = "37.tsv.gz")
  for (j in 1:length(out)) {
    # j=1
    outdat <- out[j]
    print(outdat)
    
    expname <- unique(exposure$id.exposure)
    outname <- trimws(substr(outdat,1,12))
    
    outcome<-tryCatch(
      outcome <- read_outcome_data(file.path(dir,outdat), 
                                   sep = "\t", 
                                   phenotype_col = outname,
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
    
    if (!is.null(outcome) && !anyNA(exposure$eaf.exposure)) { ##如果存在NA值，则返回TRUE
      dat<-harmonise_data(exposure,outcome)
      out_steiger <- directionality_test(dat)
      
      steiger_filter <- as.data.frame(out_steiger)
      # 
      # threshold <- 0.05 ##定义p值范围
      target_column <- "correct_causal_direction"
      
      if (steiger_filter[[target_column]]== "TRUE") {
        # 如果存在小于阈值的元素，则执行相应的操作
        print(paste("MR Steiger test of directionality", target_column, "is TRUE"))
        
        # 在这里你可以执行其他的操作
        write.csv(steiger_filter, paste0("./filtered_steiger/",expname,"_",outname,"_mr_results.csv"))
      }
    }
    
  }
}


# 将阳性结果整合为一个表 -------------------------------------------------------------
source("step1_lib.R")


# 创建一个新文件夹来存放筛选后的文件
# dir.create("./filterered_mr_results")

# 列出文件夹中的所有 CSV 文件
dir <- "./filtered_steiger"


# 获取文件夹中的所有CSV文件名
csv_files <- list.files(dir,pattern = "*.csv")


# 设置文件夹路径
setwd("./filtered_steiger")

# 创建一个空的数据框来存储筛选后的结果
filtered_df <- data.frame()


# 遍历每个CSV文件
for (file in csv_files) {
  # 读取CSV文件
  data <- read.csv(file)
  
  # 筛选特定行，假设特定条件是某一列的值为特定值
  keyword_column <- "correct_causal_direction" # 你要筛选的列名
  keywords <- "TRUE" # 修改为你要筛选的特定值
  
  filtered_data <- data[data[[keyword_column]] %in% keywords, ]
    
    # 检查是否有符合条件的行
    if (nrow(filtered_data) > 0) {
      file_name <- basename(file)
      
      matched_string <- str_split(file_name,"_G",simplify = T)[,2]
      
      matched_string <- str_split(matched_string,"_",simplify = T)[,1]
      
      matched_string <- paste0("G",matched_string)
      
      print(matched_string)
      
      # 将特定字符串添加为新列
      filtered_data$id_outcome <- matched_string
      
      # 将筛选后的结果添加到数据框中
      filtered_df <- rbind(filtered_df, filtered_data)  
    }
   
}
write.csv(filtered_df,file = "../newout/filtered_steiger.csv",row.names = F)


