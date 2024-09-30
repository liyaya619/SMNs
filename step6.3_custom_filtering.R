# Custom code for calculating PVE, F-statistic and power
## Load libraries and data
# Input data should be output of harmonise_data function in TwoSampleMR package. ncase and ncontrol is number of cases and controls in output data

rm(list = ls())
library(tidyverse)
library(dplyr)

### 需要根据自己的数据改动
final_dat <- read.csv("./newout/final_result.csv")
## 读取数据以后发现ukb-d-C_MALE_GENITAL和ukb-d-M13_GANGLION这两个暴露因为“_G”的问题，id_outcome未能被正确命名
# 使用 split 函数将数据框拆分为多个小数据框
split_exp <- split(final_dat, final_dat$id_outcome)

# 检查拆分后的数据框列表
print(names(split_exp))

# [1] "GCST90041873" "GCST90041874" "GCST90041875" "GCST90041876" "GCST90041877"
# [6] "GCST90041878" "GCST90041880"

overall_result_list <- readRDS("./Europe_overall_phegwas_result_list.rds")

# GCST90041873 ------------------------------------------------------------
GCST90041873<- split_exp[["GCST90041873"]]

# 创建一个空列表来存储结果
result_list <- list()

for (i in 1:nrow(GCST90041873)) {
  # i=1
  expo <- GCST90041873$id.exposure[i]
  print(expo)
  
  exposure <- overall_result_list[[expo]]
  
  outcome<-tryCatch(
      outcome <- read_outcome_data("./second_malignant_neoplasm/GCST90041873_buildGRCh37.tsv.gz", 
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
     dat<-harmonise_data(exposure,outcome)
     
     # 将结果添加到结果列表中
     result_list[[length(result_list) + 1]] <- dat
}

# 将结果列表转换为一个数据框
GCST90041873_harmonised <- bind_rows(result_list)    

write_tsv(GCST90041873_harmonised,file = "./newout/GCST90041873_harmonised.tsv")

## Add Minor allele and calculate proportion of variance explained
# Use PVE equation from http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0120758.s001

harmonised<- GCST90041873_harmonised
harmonised <- harmonised %>% mutate(maf=if_else(eaf.exposure>0.5,1-eaf.exposure,eaf.exposure))
harmonised <-  harmonised %>% mutate(pve= (2 * beta.exposure ^ 2 * maf * (1 - maf)) / ((2 * beta.exposure ^ 2 * maf * (1 - maf)) + (se.exposure ^ 2 * 2 * samplesize.exposure * maf * (1 - maf))))
head(harmonised)


## Compute F-statistic
# uses F statistic equations from https://academic.oup.com/ije/article/40/3/740/741448/Power-and-instrument-strength-requirements-for

exposure_statistics <- harmonised %>% group_by(id.exposure) %>% summarise(exposure=first(exposure), n_snp=n(), r2=sum(pve), samplesize=first(samplesize.exposure)) %>% mutate(fstat=(r2 * (samplesize - 1 - n_snp)) / ((1 - r2) * n_snp))
head(exposure_statistics)


## Calculate power for a range of odds ratios
# Uses power equations from https://cnsgenomics.shinyapps.io/mRnd/
  
ors <- c(1.05,1.1,1.25,1.33,1.5)

ncase <- 801  
ncontrol <- 455547 
n <- ncase + ncontrol 
k <- ncase/n
alpha <- 0.05
z.alpha <- qnorm(1 - alpha / 2, 0, 1)


calculate_power <-  function(r2, n, k, ors){b01 <- k * (ors / (1 + k * (ors - 1)) - 1)
ncp <- (n * r2) * (k * (ors / (1 + k * (ors - 1)) - 1)) ^ 2 / (k * (1 - k) - b01 ^ 2)
power <- 1 + pnorm(-z.alpha - sqrt(ncp), 0, 1) - pnorm(z.alpha - sqrt(ncp), 0, 1)
power}

power_table <- map_dfc(ors, ~mutate(exposure_statistics, !!sym(paste0("power",.)) := calculate_power(r2,n,k,.)) %>% dplyr::select(starts_with("power"))) # apply power function to table. 
head(power_table)

## Combine power, PVE and F-statistic

stats_dat <- bind_cols(exposure_statistics, power_table)
head(stats_dat)


## Filter data
# Filter 1) missing data e.g. sample size not available 2) F-stat >10 (rule of thumb for MR analysis) 3) duplicated exposures, taking the one with highest power. Write out tables. Filtered harmonised data can then be fed into mr function in TwoSampleMR.

stats_dat %>% filter(!is.na(power1.05)) -> stats_dat_no_NA #1 

stats_dat_no_NA %>% filter(fstat>=10) -> filter_data_1 #2

filter_data_1 %>%
  group_by(exposure) %>%              
  arrange(desc(power1.05)) %>%
  slice(1) -> test_data_2 #3

write_tsv(test_data_2, "./newout/GCST90041873_power_table_filter.txt")

test_data_2 %>% ungroup %>%
  dplyr::select(1) -> power_id

pruned_harmonised_data <- left_join(power_id, harmonised, by="id.exposure")
head(pruned_harmonised_data)
write_tsv(pruned_harmonised_data, "./newout/GCST90041873_filtered_harmonised_data.txt")



# GCST90041874 ------------------------------------------------------------
GCST90041874<- split_exp[["GCST90041874"]]

# 创建一个空列表来存储结果
result_list <- list()

for (i in 1:nrow(GCST90041874)) {
  # i=1
  expo <- GCST90041874$id.exposure[i]
  print(expo)
  
  exposure <- overall_result_list[[expo]]
  
  outcome<-tryCatch(
    outcome <- read_outcome_data("./second_malignant_neoplasm/GCST90041874_buildGRCh37.tsv.gz", 
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
  dat<-harmonise_data(exposure,outcome)
  
  # 将结果添加到结果列表中
  result_list[[length(result_list) + 1]] <- dat
}

# 将结果列表转换为一个数据框
GCST90041874_harmonised <- bind_rows(result_list)    

write_tsv(GCST90041874_harmonised,file = "./newout/GCST90041874_harmonised.tsv")

## Add Minor allele and calculate proportion of variance explained
# Use PVE equation from http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0120758.s001

harmonised<- GCST90041874_harmonised
harmonised <- harmonised %>% mutate(maf=if_else(eaf.exposure>0.5,1-eaf.exposure,eaf.exposure))
harmonised <-  harmonised %>% mutate(pve= (2 * beta.exposure ^ 2 * maf * (1 - maf)) / ((2 * beta.exposure ^ 2 * maf * (1 - maf)) + (se.exposure ^ 2 * 2 * samplesize.exposure * maf * (1 - maf))))
head(harmonised)


## Compute F-statistic
# uses F statistic equations from https://academic.oup.com/ije/article/40/3/740/741448/Power-and-instrument-strength-requirements-for

exposure_statistics <- harmonised %>% group_by(id.exposure) %>% summarise(exposure=first(exposure), n_snp=n(), r2=sum(pve), samplesize=first(samplesize.exposure)) %>% mutate(fstat=(r2 * (samplesize - 1 - n_snp)) / ((1 - r2) * n_snp))
head(exposure_statistics)


## Calculate power for a range of odds ratios
# Uses power equations from https://cnsgenomics.shinyapps.io/mRnd/

ors <- c(1.05,1.1,1.25,1.33,1.5)

ncase <- 582   
ncontrol <- 455766  
n <- ncase + ncontrol 
k <- ncase/n
alpha <- 0.05
z.alpha <- qnorm(1 - alpha / 2, 0, 1)


calculate_power <-  function(r2, n, k, ors){b01 <- k * (ors / (1 + k * (ors - 1)) - 1)
ncp <- (n * r2) * (k * (ors / (1 + k * (ors - 1)) - 1)) ^ 2 / (k * (1 - k) - b01 ^ 2)
power <- 1 + pnorm(-z.alpha - sqrt(ncp), 0, 1) - pnorm(z.alpha - sqrt(ncp), 0, 1)
power}

power_table <- map_dfc(ors, ~mutate(exposure_statistics, !!sym(paste0("power",.)) := calculate_power(r2,n,k,.)) %>% dplyr::select(starts_with("power"))) # apply power function to table. 
head(power_table)

## Combine power, PVE and F-statistic

stats_dat <- bind_cols(exposure_statistics, power_table)
head(stats_dat)


## Filter data
# Filter 1) missing data e.g. sample size not available 2) F-stat >10 (rule of thumb for MR analysis) 3) duplicated exposures, taking the one with highest power. Write out tables. Filtered harmonised data can then be fed into mr function in TwoSampleMR.

stats_dat %>% filter(!is.na(power1.05)) -> stats_dat_no_NA #1 

stats_dat_no_NA %>% filter(fstat>=10) -> filter_data_1 #2

filter_data_1 %>%
  group_by(exposure) %>%              
  arrange(desc(power1.05)) %>%
  slice(1) -> test_data_2 #3

write_tsv(test_data_2, "./newout/GCST90041874_power_table_filter.txt")

test_data_2 %>% ungroup %>%
  dplyr::select(1) -> power_id

pruned_harmonised_data <- left_join(power_id, harmonised, by="id.exposure")
head(pruned_harmonised_data)
write_tsv(pruned_harmonised_data, "./newout/GCST90041874_filtered_harmonised_data.txt")



# GCST90041875 -------------------------------------------------------------------------
GCST90041875<- split_exp[["GCST90041875"]]

# 创建一个空列表来存储结果
result_list <- list()

for (i in 1:nrow(GCST90041875)) {
  # i=1
  expo <- GCST90041875$id.exposure[i]
  print(expo)
  
  exposure <- overall_result_list[[expo]]
  
  outcome<-tryCatch(
    outcome <- read_outcome_data("./second_malignant_neoplasm/GCST90041875_buildGRCh37.tsv.gz", 
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
  dat<-harmonise_data(exposure,outcome)
  
  # 将结果添加到结果列表中
  result_list[[length(result_list) + 1]] <- dat
}

# 将结果列表转换为一个数据框
GCST90041875_harmonised <- bind_rows(result_list)    

write_tsv(GCST90041875_harmonised,file = "./newout/GCST90041875_harmonised.tsv")

## Add Minor allele and calculate proportion of variance explained
# Use PVE equation from http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0120758.s001

harmonised<- GCST90041875_harmonised
harmonised <- harmonised %>% mutate(maf=if_else(eaf.exposure>0.5,1-eaf.exposure,eaf.exposure))
harmonised <-  harmonised %>% mutate(pve= (2 * beta.exposure ^ 2 * maf * (1 - maf)) / ((2 * beta.exposure ^ 2 * maf * (1 - maf)) + (se.exposure ^ 2 * 2 * samplesize.exposure * maf * (1 - maf))))
head(harmonised)


## Compute F-statistic
# uses F statistic equations from https://academic.oup.com/ije/article/40/3/740/741448/Power-and-instrument-strength-requirements-for

exposure_statistics <- harmonised %>% group_by(id.exposure) %>% summarise(exposure=first(exposure), n_snp=n(), r2=sum(pve), samplesize=first(samplesize.exposure)) %>% mutate(fstat=(r2 * (samplesize - 1 - n_snp)) / ((1 - r2) * n_snp))
head(exposure_statistics)


## Calculate power for a range of odds ratios
# Uses power equations from https://cnsgenomics.shinyapps.io/mRnd/

ors <- c(1.05,1.1,1.25,1.33,1.5)

ncase <- 835    
ncontrol <- 455513  
n <- ncase + ncontrol 
k <- ncase/n
alpha <- 0.05
z.alpha <- qnorm(1 - alpha / 2, 0, 1)


calculate_power <-  function(r2, n, k, ors){b01 <- k * (ors / (1 + k * (ors - 1)) - 1)
ncp <- (n * r2) * (k * (ors / (1 + k * (ors - 1)) - 1)) ^ 2 / (k * (1 - k) - b01 ^ 2)
power <- 1 + pnorm(-z.alpha - sqrt(ncp), 0, 1) - pnorm(z.alpha - sqrt(ncp), 0, 1)
power}

power_table <- map_dfc(ors, ~mutate(exposure_statistics, !!sym(paste0("power",.)) := calculate_power(r2,n,k,.)) %>% dplyr::select(starts_with("power"))) # apply power function to table. 
head(power_table)

## Combine power, PVE and F-statistic

stats_dat <- bind_cols(exposure_statistics, power_table)
head(stats_dat)


## Filter data
# Filter 1) missing data e.g. sample size not available 2) F-stat >10 (rule of thumb for MR analysis) 3) duplicated exposures, taking the one with highest power. Write out tables. Filtered harmonised data can then be fed into mr function in TwoSampleMR.

stats_dat %>% filter(!is.na(power1.05)) -> stats_dat_no_NA #1 

stats_dat_no_NA %>% filter(fstat>=10) -> filter_data_1 #2

filter_data_1 %>%
  group_by(exposure) %>%              
  arrange(desc(power1.05)) %>%
  slice(1) -> test_data_2 #3

write_tsv(test_data_2, "./newout/GCST90041875_power_table_filter.txt")

test_data_2 %>% ungroup %>%
  dplyr::select(1) -> power_id

pruned_harmonised_data <- left_join(power_id, harmonised, by="id.exposure")
head(pruned_harmonised_data)
write_tsv(pruned_harmonised_data, "./newout/GCST90041875_filtered_harmonised_data.txt")



# GCST90041876 -------------------------------------------------------------------------
GCST90041876<- split_exp[["GCST90041876"]]

# 创建一个空列表来存储结果
result_list <- list()

for (i in 1:nrow(GCST90041876)) {
  # i=1
  expo <- GCST90041876$id.exposure[i]
  print(expo)
  
  exposure <- overall_result_list[[expo]]
  
  outcome<-tryCatch(
    outcome <- read_outcome_data("./second_malignant_neoplasm/GCST90041876_buildGRCh37.tsv.gz", 
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
  dat<-harmonise_data(exposure,outcome)
  
  # 将结果添加到结果列表中
  result_list[[length(result_list) + 1]] <- dat
}

# 将结果列表转换为一个数据框
GCST90041876_harmonised <- bind_rows(result_list)    

write_tsv(GCST90041876_harmonised,file = "./newout/GCST90041876_harmonised.tsv")

## Add Minor allele and calculate proportion of variance explained
# Use PVE equation from http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0120758.s001

harmonised<- GCST90041876_harmonised
harmonised <- harmonised %>% mutate(maf=if_else(eaf.exposure>0.5,1-eaf.exposure,eaf.exposure))
harmonised <-  harmonised %>% mutate(pve= (2 * beta.exposure ^ 2 * maf * (1 - maf)) / ((2 * beta.exposure ^ 2 * maf * (1 - maf)) + (se.exposure ^ 2 * 2 * samplesize.exposure * maf * (1 - maf))))
head(harmonised)


## Compute F-statistic
# uses F statistic equations from https://academic.oup.com/ije/article/40/3/740/741448/Power-and-instrument-strength-requirements-for

exposure_statistics <- harmonised %>% group_by(id.exposure) %>% summarise(exposure=first(exposure), n_snp=n(), r2=sum(pve), samplesize=first(samplesize.exposure)) %>% mutate(fstat=(r2 * (samplesize - 1 - n_snp)) / ((1 - r2) * n_snp))
head(exposure_statistics)


## Calculate power for a range of odds ratios
# Uses power equations from https://cnsgenomics.shinyapps.io/mRnd/

ors <- c(1.05,1.1,1.25,1.33,1.5)

ncase <- 464    
ncontrol <- 455884   
n <- ncase + ncontrol 
k <- ncase/n
alpha <- 0.05
z.alpha <- qnorm(1 - alpha / 2, 0, 1)


calculate_power <-  function(r2, n, k, ors){b01 <- k * (ors / (1 + k * (ors - 1)) - 1)
ncp <- (n * r2) * (k * (ors / (1 + k * (ors - 1)) - 1)) ^ 2 / (k * (1 - k) - b01 ^ 2)
power <- 1 + pnorm(-z.alpha - sqrt(ncp), 0, 1) - pnorm(z.alpha - sqrt(ncp), 0, 1)
power}

power_table <- map_dfc(ors, ~mutate(exposure_statistics, !!sym(paste0("power",.)) := calculate_power(r2,n,k,.)) %>% dplyr::select(starts_with("power"))) # apply power function to table. 
head(power_table)

## Combine power, PVE and F-statistic

stats_dat <- bind_cols(exposure_statistics, power_table)
head(stats_dat)


## Filter data
# Filter 1) missing data e.g. sample size not available 2) F-stat >10 (rule of thumb for MR analysis) 3) duplicated exposures, taking the one with highest power. Write out tables. Filtered harmonised data can then be fed into mr function in TwoSampleMR.

stats_dat %>% filter(!is.na(power1.05)) -> stats_dat_no_NA #1 

stats_dat_no_NA %>% filter(fstat>=10) -> filter_data_1 #2

filter_data_1 %>%
  group_by(exposure) %>%              
  arrange(desc(power1.05)) %>%
  slice(1) -> test_data_2 #3

write_tsv(test_data_2, "./newout/GCST90041876_power_table_filter.txt")

test_data_2 %>% ungroup %>%
  dplyr::select(1) -> power_id

pruned_harmonised_data <- left_join(power_id, harmonised, by="id.exposure")
head(pruned_harmonised_data)
write_tsv(pruned_harmonised_data, "./newout/GCST90041876_filtered_harmonised_data.txt")



# GCST90041877 -------------------------------------------------------------------------
GCST90041877<- split_exp[["GCST90041877"]]

# 创建一个空列表来存储结果
result_list <- list()

for (i in 1:nrow(GCST90041877)) {
  # i=1
  expo <- GCST90041877$id.exposure[i]
  print(expo)
  
  exposure <- overall_result_list[[expo]]
  
  outcome<-tryCatch(
    outcome <- read_outcome_data("./second_malignant_neoplasm/GCST90041877_buildGRCh37.tsv.gz", 
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
  dat<-harmonise_data(exposure,outcome)
  
  # 将结果添加到结果列表中
  result_list[[length(result_list) + 1]] <- dat
}

# 将结果列表转换为一个数据框
GCST90041877_harmonised <- bind_rows(result_list)    

write_tsv(GCST90041877_harmonised,file = "./newout/GCST90041877_harmonised.tsv")

## Add Minor allele and calculate proportion of variance explained
# Use PVE equation from http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0120758.s001

harmonised<- GCST90041877_harmonised
harmonised <- harmonised %>% mutate(maf=if_else(eaf.exposure>0.5,1-eaf.exposure,eaf.exposure))
harmonised <-  harmonised %>% mutate(pve= (2 * beta.exposure ^ 2 * maf * (1 - maf)) / ((2 * beta.exposure ^ 2 * maf * (1 - maf)) + (se.exposure ^ 2 * 2 * samplesize.exposure * maf * (1 - maf))))
head(harmonised)


## Compute F-statistic
# uses F statistic equations from https://academic.oup.com/ije/article/40/3/740/741448/Power-and-instrument-strength-requirements-for

exposure_statistics <- harmonised %>% group_by(id.exposure) %>% summarise(exposure=first(exposure), n_snp=n(), r2=sum(pve), samplesize=first(samplesize.exposure)) %>% mutate(fstat=(r2 * (samplesize - 1 - n_snp)) / ((1 - r2) * n_snp))
head(exposure_statistics)


## Calculate power for a range of odds ratios
# Uses power equations from https://cnsgenomics.shinyapps.io/mRnd/

ors <- c(1.05,1.1,1.25,1.33,1.5)

ncase <- 927    
ncontrol <- 455421   
n <- ncase + ncontrol 
k <- ncase/n
alpha <- 0.05
z.alpha <- qnorm(1 - alpha / 2, 0, 1)


calculate_power <-  function(r2, n, k, ors){b01 <- k * (ors / (1 + k * (ors - 1)) - 1)
ncp <- (n * r2) * (k * (ors / (1 + k * (ors - 1)) - 1)) ^ 2 / (k * (1 - k) - b01 ^ 2)
power <- 1 + pnorm(-z.alpha - sqrt(ncp), 0, 1) - pnorm(z.alpha - sqrt(ncp), 0, 1)
power}

power_table <- map_dfc(ors, ~mutate(exposure_statistics, !!sym(paste0("power",.)) := calculate_power(r2,n,k,.)) %>% dplyr::select(starts_with("power"))) # apply power function to table. 
head(power_table)

## Combine power, PVE and F-statistic

stats_dat <- bind_cols(exposure_statistics, power_table)
head(stats_dat)


## Filter data
# Filter 1) missing data e.g. sample size not available 2) F-stat >10 (rule of thumb for MR analysis) 3) duplicated exposures, taking the one with highest power. Write out tables. Filtered harmonised data can then be fed into mr function in TwoSampleMR.

stats_dat %>% filter(!is.na(power1.05)) -> stats_dat_no_NA #1 

stats_dat_no_NA %>% filter(fstat>=10) -> filter_data_1 #2

filter_data_1 %>%
  group_by(exposure) %>%              
  arrange(desc(power1.05)) %>%
  slice(1) -> test_data_2 #3

write_tsv(test_data_2, "./newout/GCST90041877_power_table_filter.txt")

test_data_2 %>% ungroup %>%
  dplyr::select(1) -> power_id

pruned_harmonised_data <- left_join(power_id, harmonised, by="id.exposure")
head(pruned_harmonised_data)
write_tsv(pruned_harmonised_data, "./newout/GCST90041877_filtered_harmonised_data.txt")



# GCST90041878 -------------------------------------------------------------------------
GCST90041878<- split_exp[["GCST90041878"]]

# 创建一个空列表来存储结果
result_list <- list()

for (i in 1:nrow(GCST90041878)) {
  # i=1
  expo <- GCST90041878$id.exposure[i]
  print(expo)
  
  exposure <- overall_result_list[[expo]]
  
  outcome<-tryCatch(
    outcome <- read_outcome_data("./second_malignant_neoplasm/GCST90041878_buildGRCh37.tsv.gz", 
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
  dat<-harmonise_data(exposure,outcome)
  
  # 将结果添加到结果列表中
  result_list[[length(result_list) + 1]] <- dat
}

# 将结果列表转换为一个数据框
GCST90041878_harmonised <- bind_rows(result_list)    

write_tsv(GCST90041878_harmonised,file = "./newout/GCST90041878_harmonised.tsv")

## Add Minor allele and calculate proportion of variance explained
# Use PVE equation from http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0120758.s001

harmonised<- GCST90041878_harmonised
harmonised <- harmonised %>% mutate(maf=if_else(eaf.exposure>0.5,1-eaf.exposure,eaf.exposure))
harmonised <-  harmonised %>% mutate(pve= (2 * beta.exposure ^ 2 * maf * (1 - maf)) / ((2 * beta.exposure ^ 2 * maf * (1 - maf)) + (se.exposure ^ 2 * 2 * samplesize.exposure * maf * (1 - maf))))
head(harmonised)


## Compute F-statistic
# uses F statistic equations from https://academic.oup.com/ije/article/40/3/740/741448/Power-and-instrument-strength-requirements-for

exposure_statistics <- harmonised %>% group_by(id.exposure) %>% summarise(exposure=first(exposure), n_snp=n(), r2=sum(pve), samplesize=first(samplesize.exposure)) %>% mutate(fstat=(r2 * (samplesize - 1 - n_snp)) / ((1 - r2) * n_snp))
head(exposure_statistics)


## Calculate power for a range of odds ratios
# Uses power equations from https://cnsgenomics.shinyapps.io/mRnd/

ors <- c(1.05,1.1,1.25,1.33,1.5)

ncase <- 107    
ncontrol <- 456241  
n <- ncase + ncontrol 
k <- ncase/n
alpha <- 0.05
z.alpha <- qnorm(1 - alpha / 2, 0, 1)


calculate_power <-  function(r2, n, k, ors){b01 <- k * (ors / (1 + k * (ors - 1)) - 1)
ncp <- (n * r2) * (k * (ors / (1 + k * (ors - 1)) - 1)) ^ 2 / (k * (1 - k) - b01 ^ 2)
power <- 1 + pnorm(-z.alpha - sqrt(ncp), 0, 1) - pnorm(z.alpha - sqrt(ncp), 0, 1)
power}

power_table <- map_dfc(ors, ~mutate(exposure_statistics, !!sym(paste0("power",.)) := calculate_power(r2,n,k,.)) %>% dplyr::select(starts_with("power"))) # apply power function to table. 
head(power_table)

## Combine power, PVE and F-statistic

stats_dat <- bind_cols(exposure_statistics, power_table)
head(stats_dat)


## Filter data
# Filter 1) missing data e.g. sample size not available 2) F-stat >10 (rule of thumb for MR analysis) 3) duplicated exposures, taking the one with highest power. Write out tables. Filtered harmonised data can then be fed into mr function in TwoSampleMR.

stats_dat %>% filter(!is.na(power1.05)) -> stats_dat_no_NA #1 

stats_dat_no_NA %>% filter(fstat>=10) -> filter_data_1 #2

filter_data_1 %>%
  group_by(exposure) %>%              
  arrange(desc(power1.05)) %>%
  slice(1) -> test_data_2 #3

write_tsv(test_data_2, "./newout/GCST90041878_power_table_filter.txt")

test_data_2 %>% ungroup %>%
  dplyr::select(1) -> power_id

pruned_harmonised_data <- left_join(power_id, harmonised, by="id.exposure")
head(pruned_harmonised_data)
write_tsv(pruned_harmonised_data, "./newout/GCST90041878_filtered_harmonised_data.txt")



# GCST90041880 -------------------------------------------------------------------------
GCST90041880<- split_exp[["GCST90041880"]]

# 创建一个空列表来存储结果
result_list <- list()

for (i in 1:nrow(GCST90041880)) {
  # i=1
  expo <- GCST90041880$id.exposure[i]
  print(expo)
  
  exposure <- overall_result_list[[expo]]
  
  outcome<-tryCatch(
    outcome <- read_outcome_data("./second_malignant_neoplasm/GCST90041880_buildGRCh37.tsv.gz", 
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
  dat<-harmonise_data(exposure,outcome)
  
  # 将结果添加到结果列表中
  result_list[[length(result_list) + 1]] <- dat
}

# 将结果列表转换为一个数据框
GCST90041880_harmonised <- bind_rows(result_list)    

write_tsv(GCST90041880_harmonised,file = "./newout/GCST90041880_harmonised.tsv")

## Add Minor allele and calculate proportion of variance explained
# Use PVE equation from http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0120758.s001

harmonised<- GCST90041880_harmonised
harmonised <- harmonised %>% mutate(maf=if_else(eaf.exposure>0.5,1-eaf.exposure,eaf.exposure))
harmonised <-  harmonised %>% mutate(pve= (2 * beta.exposure ^ 2 * maf * (1 - maf)) / ((2 * beta.exposure ^ 2 * maf * (1 - maf)) + (se.exposure ^ 2 * 2 * samplesize.exposure * maf * (1 - maf))))
head(harmonised)


## Compute F-statistic
# uses F statistic equations from https://academic.oup.com/ije/article/40/3/740/741448/Power-and-instrument-strength-requirements-for

exposure_statistics <- harmonised %>% group_by(id.exposure) %>% summarise(exposure=first(exposure), n_snp=n(), r2=sum(pve), samplesize=first(samplesize.exposure)) %>% mutate(fstat=(r2 * (samplesize - 1 - n_snp)) / ((1 - r2) * n_snp))
head(exposure_statistics)


## Calculate power for a range of odds ratios
# Uses power equations from https://cnsgenomics.shinyapps.io/mRnd/

ors <- c(1.05,1.1,1.25,1.33,1.5)

ncase <- 274    
ncontrol <- 456074   
n <- ncase + ncontrol 
k <- ncase/n
alpha <- 0.05
z.alpha <- qnorm(1 - alpha / 2, 0, 1)


calculate_power <-  function(r2, n, k, ors){b01 <- k * (ors / (1 + k * (ors - 1)) - 1)
ncp <- (n * r2) * (k * (ors / (1 + k * (ors - 1)) - 1)) ^ 2 / (k * (1 - k) - b01 ^ 2)
power <- 1 + pnorm(-z.alpha - sqrt(ncp), 0, 1) - pnorm(z.alpha - sqrt(ncp), 0, 1)
power}

power_table <- map_dfc(ors, ~mutate(exposure_statistics, !!sym(paste0("power",.)) := calculate_power(r2,n,k,.)) %>% dplyr::select(starts_with("power"))) # apply power function to table. 
head(power_table)

## Combine power, PVE and F-statistic

stats_dat <- bind_cols(exposure_statistics, power_table)
head(stats_dat)


## Filter data
# Filter 1) missing data e.g. sample size not available 2) F-stat >10 (rule of thumb for MR analysis) 3) duplicated exposures, taking the one with highest power. Write out tables. Filtered harmonised data can then be fed into mr function in TwoSampleMR.

stats_dat %>% filter(!is.na(power1.05)) -> stats_dat_no_NA #1 

stats_dat_no_NA %>% filter(fstat>=10) -> filter_data_1 #2

filter_data_1 %>%
  group_by(exposure) %>%              
  arrange(desc(power1.05)) %>%
  slice(1) -> test_data_2 #3

write_tsv(test_data_2, "./newout/GCST90041880_power_table_filter.txt")

test_data_2 %>% ungroup %>%
  dplyr::select(1) -> power_id

pruned_harmonised_data <- left_join(power_id, harmonised, by="id.exposure")
head(pruned_harmonised_data)
write_tsv(pruned_harmonised_data, "./newout/GCST90041880_filtered_harmonised_data.txt")



# combined ----------------------------------------------------------------
file_list <- list.files("./newout/",pattern = "filtered_harmonised_data.txt",full.names = T)

# 读取文件并合并为一个数据框
merged_df <- do.call(rbind, lapply(file_list, read_delim))
