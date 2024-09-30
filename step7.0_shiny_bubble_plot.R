rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyverse)
library(openxlsx)

dat <- readxl::read_xlsx("./newout/tidy_exposure_outcome.xlsx",sheet = 1) %>% filter(pval<0.05) 
# 将 mr_egger_pval 列中的 NA 值替换为 0.99
dat$mr_egger_pval[is.na(dat$mr_egger_pval)] <- 0.99
dat <- dat %>% filter(mr_egger_pval>0.05) 


# dat$beta <- ifelse(dat$b>0,"Increasing risk","Decreasing risk")
dat$`-log10pval` <- -log10(dat$pval)
dat$Type_manual_added <- trimws(str_to_title(dat$Type_manual_added))
dat$Outcome <- trimws(dat$Outcome)
dat$beta <- dat$b

# Cancer ---------------------------------------------
pdat <- dat[,c("Exposure","Outcome","Type_manual_added","beta","-log10pval")]

table(pdat$Outcome)
table(pdat$Type_manual_added)
sdat <- pdat[pdat$Type_manual_added=="First Primary Cancer",]

# 创建气泡图，根据 p 值定义点的大小
p1 <- ggplot(sdat, aes(x = Outcome, y = Exposure, color = beta, size = `-log10pval`)) +
  geom_point(stroke = 0.5) +  # 设置描边粗细
  scale_color_gradient2(low = "#215284", mid = "white", high = "#C1151F", midpoint = 0) +  # 根据beta值进行颜色渐变
  scale_size_continuous(range = c(1, 10)) +  # 设置点大小的范围
  scale_x_discrete(expand = c(0, 0.4), position = "bottom") +  # 设置X轴数据的边距和间距，X轴位于底端
  labs(x = "", y = "", size = "-log10(p)") +  # 取消XY轴标题，并为点大小添加标签
  theme(
    axis.ticks = element_blank(),  # 删除轴刻度
    legend.position = "top",  # 删除图例
    axis.text.x = element_text(angle = 45, hjust = 1),  # 调整X轴标签角度为45度
    axis.text.y = element_text(size = 8),  
    panel.grid.major = element_blank(),  # 删除网格
    panel.grid.minor = element_blank(),  # 删除网格
    panel.background = element_blank(),  # 删除背景
    axis.line.y = element_line(color = "black", size = 0.4)  # 设置左侧边框线条
  )
p1
ggsave(p1,file="./newout//bubble.example.pdf",width = 12,height = 7)

# start shiny -------------------------------------------------------------
source("step7.1_deploy_shiny.R")
