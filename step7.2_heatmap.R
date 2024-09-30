rm(list = ls())

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(pheatmap)

dat <- readxl::read_xlsx("./newout/tidy_exposure_outcome.xlsx", sheet = 1) %>%
  filter(pval < 0.05)

# 将 mr_egger_pval 列中的 NA 值替换为 0.99
dat$mr_egger_pval[is.na(dat$mr_egger_pval)] <- 0.99

# 先分配 "Robust"
dat$Evidence <- NA
dat$Evidence[dat$pval < 1.359434e-05 & dat$mr_egger_pval > 0.05] <- "Robust"

# 再分配 "Probable"
dat$Evidence[is.na(dat$Evidence) & dat$pval < 0.05 & dat$mr_egger_pval > 0.05&
               (dat$mr_WME_pval < 0.05 | dat$mr_MBE_pval < 0.05) & 
               dat$nsnp > 1] <- "Probable"

dat$Evidence[is.na(dat$Evidence) & dat$pval < 0.05 & dat$nsnp == 1] <- "Probable"

# 分配 "Suggestive"
dat$Evidence[is.na(dat$Evidence) & dat$pval < 0.05 & dat$nsnp > 1] <- "Suggestive"
dat$Evidence[is.na(dat$Evidence) & dat$pval < 0.05 & dat$nsnp == 1] <- "Probable"

table(dat$Evidence)
# Probable     Robust Suggestive 
# 845          4        898  

table(dat$Type_manual_added[dat$Evidence == "Probable"])


####  绘图############
dat <- dat %>%  filter(mr_egger_pval>0.05) %>% tidyr::unite(.,"Traits",`id.exposure`,`Exposure`,sep = ":",remove = TRUE)
colnames(dat)

pdat <- dat[,c("Traits","Outcome","b")] 
pdat$Outcome <- trimws(str_to_title(pdat$Outcome))

##这里会报错是因为有些traits虽然ID不一样，但是traits一致，因此需要用ID:traits标识来区分
wdat <- spread(pdat,key = "Outcome",value = "b",fill = 0) %>% as.data.frame()
rownames(wdat) <- wdat$Traits
wdat <- wdat[,-1]

range(wdat)

dat$Category <- dat$Type_manual_added

Group <- dat[,c("Traits","Category")]
library(dplyr)
Group_sorted <- arrange(Group, Category)

annotation_row <- as.data.frame(unique(Group_sorted))
rownames(annotation_row) <- annotation_row$Traits
annotation_row <- annotation_row[,"Category",drop=FALSE] ##选取一列仍然保持是data.frame

wdat <-wdat[match(rownames(annotation_row), rownames(wdat)), ]

pdf("./newout/pheatmap.pdf",width =17,height = 8 )
pheatmap(wdat,
         scale = "row",
         show_colnames = T,
         show_rownames = F,
         cluster_rows = F,
         cluster_cols = F,
         cutree_rows = T,
         annotation_row = annotation_row,
         annotation_legend = T,
         legend = T,
         legend_breaks = seq(-100,75,2),
         treeheight_row=0 ,
         treeheight_col=0,
         fontsize=10,angle_col = 45)  # 将列名旋转 45 度
dev.off()

