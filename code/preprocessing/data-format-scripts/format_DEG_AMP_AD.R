library(tidyverse)
zap()




########
## DEG 
########

# READ AND FORMAT
deg <- read.csv("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/AMP-AD differentially expressed genes (DEG)/deg_meta_analysis.csv")
names(deg) <- c("gene_id","te_random","pvalue_fdr","direction")
# ADD TRAIT
deg$trait <- "AD Diagnosis (males and females)"
deg$tissue <- "brain"
deg <- deg %>% mutate(direction=tolower(direction)) 
# WRITE
write.csv(deg,"DEG/logsdon_deg_meta_analysis.csv",row.names = F,quote = T,na = "")

#####################

# FILE INFO
data = c(
  "file_name","logsdon_deg_meta_analysis.csv",
  "sid", "logsdon_deg_meta_analysis",
  "organism", "Homo Sapiens",
  "node_type", "DEG",
  "cohort", "Mayo|MSBB|ROSMAP",
  "publication","logsdon_2019",
  "significance","FDR(p)<=0.05",
  "tax","hsa",
  "taxID","9606"
)
# MAKE DF
info <- data.frame(matrix(data = data, ncol = 2, byrow = T))
# WRITE
write.table(
  info,
  file = "DEG/logsdon_deg_meta_analysis_info.csv",
  sep = ",",
  col.names = F,
  row.names = F
)
