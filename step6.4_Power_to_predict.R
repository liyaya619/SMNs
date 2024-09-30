rm(list = ls())
library(tidyverse)
library(dplyr)

# [1] "GCST90041873" "GCST90041874" "GCST90041875" "GCST90041876" "GCST90041877"
# [6] "GCST90041878" "GCST90041880"


# GCST90041873 ------------------------------------------------------------
dat <- read_delim("./newout/GCST90041873_power_table_filter.txt")
pdat <- dat %>% dplyr::select(c(2,6:11)) %>% 
  gather(key = "Odds Ratio", value = "Power %",  3:7) 


pdat$`Power %` <- 100*pdat$`Power %`
pdat$`Odds Ratio` <- as.numeric(gsub("power","",pdat$`Odds Ratio`))

# 创建图表
p <- ggplot(pdat, aes(x = `Odds Ratio`, y = `Power %`, color = log(fstat), group = exposure)) +
  geom_smooth(se = FALSE, method = "loess", size = 0.3, alpha = 0.5) +
  scale_color_viridis_c(name = "log(F-stat)") +
  scale_x_continuous(breaks = seq(1.0, 1.5, by = 0.1)) +
  labs(x = NULL, y = NULL,title = "Respiratory Organs") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    axis.text.x = element_blank(),
    plot.title = element_text(size=12,hjust=0.5),legend.position="none"
  )

p ##第一行第一列



# GCST90041874 ------------------------------------------------------------
dat1 <- read_delim("./newout/GCST90041874_power_table_filter.txt")
pdat1 <- dat1 %>% dplyr::select(c(2,6:11)) %>% 
  gather(key = "Odds Ratio", value = "Power %",  3:7) 


pdat1$`Power %` <- 100*pdat1$`Power %`
pdat1$`Odds Ratio` <- as.numeric(gsub("power","",pdat1$`Odds Ratio`))

# 创建图表
p1 <- ggplot(pdat1, aes(x = `Odds Ratio`, y = `Power %`, color = log(fstat), group = exposure)) +
  geom_smooth(se = FALSE, method = "loess", size = 0.3, alpha = 0.5) +
  scale_color_viridis_c(name = "log(F-stat)") +
  scale_x_continuous(breaks = seq(1.0, 1.5, by = 0.1)) +
  labs(x = NULL, y = NULL,title = "Digestive Systems") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    plot.title = element_text(size=12,hjust=0.5),legend.position="none"
  )
p1 ##第二行第一列




# GCST90041875 ------------------------------------------------------------
dat2 <- read_delim("./newout/GCST90041875_power_table_filter.txt")
pdat2 <- dat2 %>% dplyr::select(c(2,6:11)) %>% 
  gather(key = "Odds Ratio", value = "Power %",  3:7) 


pdat2$`Power %` <- 100*pdat2$`Power %`
pdat2$`Odds Ratio` <- as.numeric(gsub("power","",pdat2$`Odds Ratio`))

# 创建图表
p2 <- ggplot(pdat2, aes(x = `Odds Ratio`, y = `Power %`, color = log(fstat), group = exposure)) +
  geom_smooth(se = FALSE, method = "loess", size = 0.3, alpha = 0.5) +
  scale_color_viridis_c(name = "log(F-stat)") +
  scale_x_continuous(breaks = seq(1.0, 1.5, by = 0.1)) +
  labs(x = NULL, y = NULL,title = "Liver") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(size=12,hjust=0.5),legend.position="none"
  )
p2 ##第一行第二列



# GCST90041876 ------------------------------------------------------------
dat3 <- read_delim("./newout/GCST90041876_power_table_filter.txt")
pdat3 <- dat3 %>% dplyr::select(c(2,6:11)) %>% 
  gather(key = "Odds Ratio", value = "Power %",  3:7) 


pdat3$`Power %` <- 100*pdat3$`Power %`
pdat3$`Odds Ratio` <- as.numeric(gsub("power","",pdat3$`Odds Ratio`))

# 创建图表
p3 <- ggplot(pdat3, aes(x = `Odds Ratio`, y = `Power %`, color = log(fstat), group = exposure)) +
  geom_smooth(se = FALSE, method = "loess", size = 0.3, alpha = 0.5) +
  scale_color_viridis_c(name = "log(F-stat)") +
  scale_x_continuous(breaks = seq(1.0, 1.5, by = 0.1)) +
  labs(x = NULL, y = NULL,title = "Brain/Spine") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(size=12,hjust=0.5),legend.position="none"
  )
p3 ##第一行第三列



# GCST90041877 ------------------------------------------------------------
dat4 <- read_delim("./newout/GCST90041877_power_table_filter.txt")
pdat4 <- dat4 %>% dplyr::select(c(2,6:11)) %>% 
  gather(key = "Odds Ratio", value = "Power %",  3:7) 


pdat4$`Power %` <- 100*pdat4$`Power %`
pdat4$`Odds Ratio` <- as.numeric(gsub("power","",pdat4$`Odds Ratio`))

# 创建图表
p4 <- ggplot(pdat4, aes(x = `Odds Ratio`, y = `Power %`, color = log(fstat), group = exposure)) +
  geom_smooth(se = FALSE, method = "loess", size = 0.3, alpha = 0.5) +
  scale_color_viridis_c(name = "log(F-stat)") +
  scale_x_continuous(breaks = seq(1.0, 1.5, by = 0.1)) +
  labs(x = NULL, y = NULL,title = "Bone") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    axis.text.y = element_blank(),
    plot.title = element_text(size=12,hjust=0.5),legend.position="none"
  )
p4 ##第二行第二列




# GCST90041878 -------------------------------------------------------------------------
dat5 <- read_delim("./newout/GCST90041878_power_table_filter.txt")
pdat5 <- dat5 %>% dplyr::select(c(2,6:11)) %>% 
  gather(key = "Odds Ratio", value = "Power %",  3:7) 


pdat5$`Power %` <- 100*pdat5$`Power %`
pdat5$`Odds Ratio` <- as.numeric(gsub("power","",pdat5$`Odds Ratio`))

# 创建图表
p5 <- ggplot(pdat5, aes(x = `Odds Ratio`, y = `Power %`, color = log(fstat), group = exposure)) +
  geom_smooth(se = FALSE, method = "loess", size = 0.3, alpha = 0.5) +
  scale_color_viridis_c(name = "log(F-stat)") +
  scale_x_continuous(breaks = seq(1.0, 1.5, by = 0.1)) +
  labs(x = NULL, y = NULL,title = "Skin") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    axis.text.y = element_blank(),
    plot.title = element_text(size=12,hjust=0.5),legend.position="none"
  )
p5 ##第二行第三列



# GCST90041880 ------------------------------------------------------------
dat6 <- read_delim("./newout/GCST90041880_power_table_filter.txt")
pdat6 <- dat6 %>% dplyr::select(c(2,6:11)) %>% 
  gather(key = "Odds Ratio", value = "Power %",  3:7) 


pdat6$`Power %` <- 100*pdat6$`Power %`
pdat6$`Odds Ratio` <- as.numeric(gsub("power","",pdat6$`Odds Ratio`))

# 创建图表
p6 <- ggplot(pdat6, aes(x = `Odds Ratio`, y = `Power %`, color = log(fstat), group = exposure)) +
  geom_smooth(se = FALSE, method = "loess", size = 0.3, alpha = 0.5) +
  scale_color_viridis_c(name = "log(F-stat)") +
  scale_x_continuous(breaks = seq(1.0, 1.5, by = 0.1)) +
  labs(x = NULL, y = NULL,title = "Secondary Malignant Neoplasm") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    plot.title = element_text(size=12,hjust=0.5)
  )
p6 ##fig a

## 拼图
library(ggplot2)
library(cowplot)
pr <- plot_grid(p, p2, p3, p1,p4,p5, ncol=3)
pr

pl <- p6;pl

pf <- plot_grid(
  pl, pr,nrow = 1, rel_widths = c(1, 2)
)
pf

ggsave(pf,"./newout/Power_to_predict.pdf",width = 15,height = 8)
