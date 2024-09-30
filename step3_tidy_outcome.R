source("step1_lib.R")
options(scipen = 100) # 小数点后100位不使用科学计数法

if (F) {
  dir <- "./second_malignant_neoplasm"
  dat_38 <- list.files(dir,pattern = "38.tsv.gz") 
  
  newfile <- str_split(dat_38,"_",simplify = T)[,1]
  
  # hg19
  # [1] "chromosome"              "variant_id"             
  # [3] "base_pair_location"      "effect_allele"          
  # [5] "other_allele"            "N"                      
  # [7] "effect_allele_frequency" "T"                      
  # [9] "SE_T"                    "P_noSPA"                
  # [11] "beta"                    "standard_error"         
  # [13] "p_value"                 "CONVERGE"
  
  # hg38
  # [1] "Name"                    "chromosome"             
  # [3] "base_pair_location"      "other_allele"           
  # [5] "effect_allele"           "Trait"                  
  # [7] "Cohort"                  "Model"                  
  # [9] "odds_ratio"              "ci_lower"               
  # [11] "ci_upper"                "p_value"                
  # [13] "effect_allele_frequency" "standard_error" 
  
  # dir.create("./second_malignant_neoplasm/dat_38/")
  
  for (i in dat_38) {
    f= file.path(dir,i)
    print(f)
    gwas <- fread(f,data.table = F)
    gwas[1:4,1:4]
    gwas$Effect <- log(gwas$odds_ratio)
    setnames(gwas,old = c("chromosome","base_pair_location","other_allele","effect_allele","effect_allele_frequency"),new = c("CHR","BP","A1","A2","EAF")) 
    
    newfile <- str_split(i,"_",simplify = T)[,1]
    # 判断数据框中A1或A2列中是否包含指定字符串
    if (any(apply(gwas[, c("A1", "A2")], 2, function(col) any(col %in% c("A", "G", "C", "T"))))) {
      # 如果包含指定字符串，则执行相应的操作
      reformatted2 <- format_sumstats(gwas,
                                      ref_genome = "GRCh38",
                                      convert_ref_genome = "GRCh37",
                                      nThread = 2,
                                      save_path =paste0("./second_malignant_neoplasm/dat_38/",newfile,".tsv.gz"))      # 在这里你可以执行其他的操作
    } else {
      # 如果不包含指定字符串，则打印"条件不符合"
      print(paste("uncorrect A1&A2:",i))
    }
    
  }
}



# FUMA_gwas ready-------------------------------------------------------------------------
dat <- fread("./second_malignant_neoplasm/GCST90041877_buildGRCh37.tsv.gz",data.table = F)
colnames(dat)
View(head(dat))

gwas <- dat %>% rename(.,c(rsID = variant_id,
                           position = base_pair_location,
                           allele1 = effect_allele,
                           allele2 = other_allele,
                           `pvalue` = p_value,
                           Beta = beta,
                           SE = SE_T)) %>% 
  dplyr::select(c("chromosome","rsID","position","allele1","allele2","pvalue","Beta","SE","N"))

View(head(gwas))
colnames(gwas) <- trimws(colnames(gwas))
gwas$pvalue <-as.numeric(gwas$pvalue) 

write.table(gwas,file = "./second_malignant_neoplasm/bonem_gwas.txt",quote = FALSE,row.names = F)




# SMR_plot ----------------------------------------------------------------
gwas_smr <- dat %>% rename(.,c(SNP = variant_id ,
                           A1 = effect_allele ,
                           A2 = other_allele,
                           p = p_value,
                           b = beta,
                           se = SE_T,
                           freq = effect_allele_frequency,
                           n = N)) %>% 
  dplyr::select(c("SNP","A1","A2","p","b","se","freq","n"))


write.table(gwas_smr,file = "./second_malignant_neoplasm/bonem_gwas_smr.txt",quote = FALSE,row.names = F,sep = "\t")
