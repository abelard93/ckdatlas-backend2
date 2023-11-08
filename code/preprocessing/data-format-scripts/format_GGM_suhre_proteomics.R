library(tidyverse)
zap()


### ----------------------- LOAD DATA ----------------------------
load("/Users/mwoerheide/Documents/PHD/projects/AD_atlas_proteomics/suhre_et_al/data/proteomics_suhre.Rdata")



write.table(
  protein_df %>% mutate(a_uniprot = a, b_uniprot = b),
  file = "/Users/mwoerheide/Documents/PHD/projects/AD_atlas_proteomics/suhre_et_al/data/proteomics_2017_suhre.csv",
  sep = ",",
  col.names = T,
  row.names = F
)

info <- data.frame(matrix(ncol = 2,nrow = 12))
info[1, ] <- c("sid","proteomics_2017_suhre")
info[2, ] <- c("node_type","edgeCOABUN")
info[3, ] <- c("file_name","proteomics_2017_suhre.csv")
info[4, ] <- c("sample_type","plasma|fasting")
info[5, ] <- c("cohort","KORA")
info[6, ] <- c("tax","hsa")
info[7, ] <- c("taxID","9606")
info[8, ] <- c("organism","Homo sapiens")
info[9, ] <- c("publication","suhre_2017_28240269")
info[10, ] <- c("significance","Bonf(p)<=0.05")
info[11, ] <- c("status","healthy")
info[12, ] <- c("study_type","population-based")

write.table(
  info,
  file = "/Users/mwoerheide/Documents/PHD/projects/AD_atlas_proteomics/suhre_et_al/data/proteomics_2017_suhre_info.csv",
  sep = ",",
  col.names = F,
  row.names = F
)
