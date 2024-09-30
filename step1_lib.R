rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(meta)
library(ggrepel)
library(DT)
# devtools::install_github("PheWAS/PheWAS")
library(PheWAS)
# library(vroom)
# install.packages("janitor")
# library(janitor)

# Note: please update your ComplexHeatmap to the latest version!
# install.packages("devtools")
# devtools::install_github("junjunlab/ClusterGVis")
# devtools::install_github('junjunlab/scRNAtoolVis',force = TRUE)
# devtools::install_github("sajuukLyu/ggunchull", type = "source")
library(ClusterGVis)
library(scRNAtoolVis)
library(ggunchull)

# https://wei-lab.app.vumc.org/phecode ##对应PheCode 
# new.annotate.phenotype.description <- vroom::vroom("./WeiLab_Resources/original_phecodes_pheinfo.csv", .name = janitor::make_clean_names, delim = ",", col_types = c(phecode = "c", description = "c", groupnum = "double", group = "c", color = "c"))
##secondary cancer
# https://www.ebi.ac.uk/gwas/efotraits/EFO_0009812

if (T) {
  library(data.table)
  library(coloc)
  library(ggplot2)
  library(ggrepel)
  library(gwasglue)
  library(gassocplot)
  library(MRSamePopTest)
  library(TwoStepCisMR)
  library(TwoSampleMR)
  library(dplyr)
  library(stringr)
  library(MungeSumstats)
}
