library(tidyverse)
zap()

######################################
## reformat eQTL data aquired from ..
## on 25-07-2019
######################################

# READ AND FORMAT
eQTL <- data.table::fread("../code/data_transformation/raw_data/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.csv")

# rename colnames to match with GTEx names
colnames(eQTL) <- c("gene_chr",
                    "variant_pos",
                    "rs_id_dbSNP151_GRCh38p7", # not ideal (better would be rsid) 
                    "variant_id",
                    "gene_id",
                    "gene_symbol",
                    "statistic",
                    "pvalue_nominal",
                    "fdr",
                    "beta",
                    "A1",
                    "A2",
                    "A2freq",
                    "expressionIncreasingAllele", 
                    "strand",
                    "geneBiotype",
                    "geneStartPosition",
                    "geneEndPosition") 

# filter for significant eQTLs
eQTL <- eQTL %>% filter(fdr<=0.05) %>% 
  filter(rs_id_dbSNP151_GRCh38p7!=".") # ?? empty rsids 

write.csv(eQTL,"../data/QTLs/eQTL/cis/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo.csv",row.names = F,quote = T,na = "")

# FILE INFO
data = c(
  "file_name","Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo.csv",
  "sid", "Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo",
  "organism", "Homo Sapiens",
  "sample_type", "Brain Cortex",
  "node_type", "eQTL",
  "cohort", "ROSMAP|CMC|HBCC|Mayo",
  "publication","sieberts_2020_33046718",
  "significance","FDR(p)<=0.05",
  "tax","hsa",
  "taxID","9606"
)
# MAKE DF
info <- data.frame(matrix(data = data, ncol = 2, byrow = T))
# WRITE
write.table(
  info,
  file = "../data/QTLs/eQTL/cis/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_info.csv",
  sep = ",",
  col.names = F,
  row.names = F
)
