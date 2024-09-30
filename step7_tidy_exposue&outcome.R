rm(list = ls())
library(tidyverse)
library(dplyr)

### 导入最终结果
final_dat <- read.csv("./newout/final_ALL_result.csv")

##导入暴露数据对应的trait
load("Europe_available_outcomes.Rdata")

Final_dat <- final_dat %>%
  mutate(Exposure = explist$trait[match(final_dat$id.exposure,explist$id)],
         Outcome = recode(id_outcome, 
                          'GCST90041873' = 'Respiratory Organs', 
                          'GCST90041874' = 'Digestive Systems',
                          'GCST90041875' = 'Liver',
                          'GCST90041876' = 'Brain/Spine',
                          'GCST90041877' = 'Bone',
                          'GCST90041878' = 'Skin',
                          'GCST90041880' = 'Unspceified')) %>% 
  dplyr::select(-c(1,3:5))

## add catagory
library(openxlsx)
catdat <- read.xlsx("./Systematic review of Mendelian randomization studies on risk of cancer.xlsx",sheet = "File S2")

if (T) {
  Final_dat$type <- ifelse(grepl("cancer|Cancer|alignant neoplasm|carcinoma|adenocarcinoma|tumor", Final_dat$Exposure), "First primary cancer",
                           ifelse(grepl("stool|Gut microbiota|abundance", Final_dat$Exposure), "Gut microbiota",
                                  ifelse(grepl("cholesterol|HDL|transferase|Apolipoprotein|riacylglycerol|LDL|transpeptidase|triglycerides|Triglycerides", Final_dat$Exposure), "Lipid metabolism biomarkers",
                                         ifelse(grepl("Height|fat|height|Weight|Body mass index|Heel|Hair|Hip|Waist|Ankle|Arm|colour", Final_dat$Exposure), "Anthropometrics",
                                                ifelse(grepl("CD|Chemokine|count|Count|cells|Interleukin|Fc|antigen|antibody|Complement|monocyte|Eosinophil|IgD|Dendritic|HLA|Macrophage|platelet", Final_dat$Exposure), "Inflammatory biomarkers", 
                                                       ifelse(grepl("ine|X-|3-|subunit", Final_dat$Exposure),"Amino acids and derivatives",
                                                              ifelse(grepl("telomere",Final_dat$Exposure),"Circulating leukocyte telomere length",
                                                                     ifelse(grepl("Pulse|Bone mineral density| blood pressure|function|FEV1|ECG|FVC|metabolic rate|Spleen|left|right",Final_dat$Exposure),"Clinical measurements",
                                                                            ifelse(grepl("Adiponectin|Fasting|IGF|HbA1C|diabetes|insulin|leptin",Final_dat$Exposure),"Diabetes and related biomarkers",
                                                                                   ifelse(grepl("Bilirubin|zinc|Vitamin D|consumption|intake|Intake|use|Bread|eat",Final_dat$Exposure),"Dietary intake and micronutrient concentrations",
                                                                                          ifelse(grepl("fatty acids",Final_dat$Exposure),"Fatty acids and derivatives",
                                                                                                 ifelse(grepl("growth factor",Final_dat$Exposure),"Growth factors",
                                                                                                        ifelse(grepl("smok|leep|frequency|Frequency|Duration|duration|Sun|sun|degree|college|activity|walking|cigarette|Time|hours|Day|morning",Final_dat$Exposure),"Lifestyle, education and behaviour",
                                                                                                               ifelse(grepl("Methylation",Final_dat$Exposure),"Methylations",
                                                                                                                      ifelse(grepl("Age|well being",Final_dat$Exposure),"Reproductive factors",
                                                                                                                             ifelse(grepl("one",Final_dat$Exposure),"Steroids",
                                                                                                                                    ifelse(grepl("Protein|protein|molecule|necrosis factor|Urate|reticulocyte|family member|NADPH",Final_dat$Exposure),"Other metabolites/biomarkers",
                                                                                                                                           ifelse(grepl("Operat|operat",Final_dat$Exposure),"Operation History",
                                                                                                                                               ifelse(grepl("Treatment|medication|isorder|pain|asthma|pancreatitis|syndrome|angina|Osteoarthritis|disease|eye|problems|llness|njury|COVID|biliary|Cirrhosis|keratosis|Medication|sclerosis|Cholelithiasis|renia|stroke|Depression|Neuroticism|Neoplasms|besity|Hyper",Final_dat$Exposure),"Other diseases and traits","unknown")))))))))))))))))))
table(Final_dat$type)
}

# 使用 split 函数将数据框拆分为多个小数据框
split_exp <- split(Final_dat, Final_dat$id_outcome)

library(openxlsx)
sheets = list("tidy_exposure_outcome" = Final_dat,
              "GCST90041873" = split_exp[[1]],
              "GCST90041874" = split_exp[[2]],
              "GCST90041875" = split_exp[[3]],
              "GCST90041876" = split_exp[[4]],
              "GCST90041877" = split_exp[[5]],
              "GCST90041878" = split_exp[[6]],
              "GCST90041880" = split_exp[[7]])
write.xlsx(sheets,"./newout/tidy_exposure_outcome.xlsx")

dat <- readxl::read_xlsx("./newout/tidy_exposure_outcome_manual_added.xlsx",sheet = 1)
table(dat$type) ## 检查一下~
