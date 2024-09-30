source("step1_lib.R")


# 创建一个新文件夹来存放筛选后的文件
# dir.create("./filterered_mr_results")

# 列出文件夹中的所有 CSV 文件
dir <- "./mr_outputs/"

# exposure <- str_split(csv_files,"_",simplify = T)[,1]
# exposure_n <- unique(exposure)
# dat <- data.frame(matrix(ncol = 2,nrow = length(exposure_n)))
# dat$id <- exposure_n
# write.csv(dat,file = "./outputs/exposure_list.csv",row.names = F)


# 获取文件夹中的所有CSV文件名
csv_files <- list.files(dir,pattern = "*.csv")


# 设置文件夹路径
setwd("./mr_outputs")


# 创建一个空列表来存储结果
result_list <- list()

# 循环遍历每个文件

# 循环遍历每个文件
for (file in csv_files) {
  df <- read.csv(file)
  
  # 筛选特定行
  keywords <- c("Inverse variance weighted", "Wald ratio") # 修改为你要筛选的特定值
  filtered_rows <- df %>% filter(method %in% keywords)
  
  # 检查数据框是否有多行且包含 "MR Egger" 方法
  if (nrow(df) > 1 && "MR Egger" %in% df$method) {
    
    # 提取 "MR Egger" 的 p 值
    mr_egger_row <- df %>% filter(method == "MR Egger")
    if (nrow(mr_egger_row) > 0) {
      mr_egger_pval <- mr_egger_row %>% dplyr::select(pval) %>% pull()
    } else {
      mr_egger_pval <- NA
    }
  } else {
    # 如果没有 "MR Egger" 方法，将 p 值设置为 NA
    mr_egger_pval <- NA
  }
  
  # 添加 "MR Egger" 的 p 值作为新列
  filtered_rows <- filtered_rows %>% mutate(mr_egger_pval = mr_egger_pval)
  
  # 检查数据框是否有多行且包含 "Weighted median" 方法
  if (nrow(df) > 1 && "Weighted median" %in% df$method) {
    
    # 提取 "Weighted median" 的 p 值
    mr_WME_row <- df %>% filter(method == "Weighted median")
    if (nrow(mr_WME_row) > 0) {
      mr_WME_pval <- mr_WME_row %>% dplyr::select(pval) %>% pull()
    } else {
      mr_WME_pval <- NA
    }
  } else {
    # 如果没有 "Weighted median" 方法，将 p 值设置为 NA
    mr_WME_pval <- NA
  }
  
  # 添加 "Weighted median" 的 p 值作为新列
  filtered_rows <- filtered_rows %>% mutate(mr_WME_pval = mr_WME_pval)
  
  # 检查数据框是否有多行且包含 "Weighted mode" 方法
  if (nrow(df) > 1 && "Weighted mode" %in% df$method) {
    
    # 提取 "Weighted mode" 的 p 值
    mr_MBE_row <- df %>% filter(method == "Weighted mode")
    if (nrow(mr_MBE_row) > 0) {
      mr_MBE_pval <- mr_MBE_row %>% dplyr::select(pval) %>% pull()
    } else {
      mr_MBE_pval <- NA
    }
  } else {
    # 如果没有 "Weighted mode" 方法，将 p 值设置为 NA
    mr_MBE_pval <- NA
  }
  
  # 添加 "Weighted mode" 的 p 值作为新列
  filtered_rows <- filtered_rows %>% mutate(mr_MBE_pval = mr_MBE_pval)
  
  # 提取文件名中的特定字符串
  file_name <- basename(file)
  matched_string <- str_split(file_name, "_GCST", simplify = TRUE)[, 2]
  matched_string <- str_split(matched_string, "_", simplify = TRUE)[, 1]
  matched_string <- paste0("GCST", matched_string)
  
  # 将特定字符串添加为新列
  filtered_rows <- filtered_rows %>% mutate(id_outcome = matched_string)
  
  # 将结果添加到结果列表中
  result_list[[length(result_list) + 1]] <- filtered_rows
}

# 将结果列表转换为一个数据框
final_result <- bind_rows(result_list)

# 打印结果
print(final_result[1:10,])


# 如果需要，可以将结果保存为一个新的CSV文件
write.csv(final_result, file = "../newout/final_ALL_result.csv", row.names = FALSE)
