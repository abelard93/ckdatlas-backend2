library(tidyverse)
zap()

##############
# GWAS IGAP
##############

# READ AND FORMAT
## stage 1
stage1 <- data.table::fread("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/GWAS/Kunkle_etal_Stage1_results.txt") %>%
  mutate(pvalue = as.numeric(pvalue)) %>% 
  na.omit %>% 
  filter(startsWith(SNP, "rs"))# includes NA SNPs
## stage 2
stage2 <- data.table::fread("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/GWAS/Kunkle_etal_Stage2_results.txt") %>%
  mutate(
    effect_allele = toupper(effect_allele),
    non_effect_allele = toupper(non_effect_allele)
  ) %>% 
  na.omit %>% 
  filter(startsWith(SNP, "rs"))
## further analysis 3A and 3B
stage3A_meta <- data.table::fread("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/GWAS/Kunkle_etal_Stage3A_results.csv") %>%
  na.omit %>% 
  filter(startsWith(SNP, "rs"))
stage3B_meta <- data.table::fread("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/GWAS/Kunkle_etal_Stage3B_results.csv") %>% 
  na.omit %>% 
  filter(startsWith(SNP, "rs"))
## join 3A + B (no overlap)
meta <- full_join(stage3A_meta, stage3B_meta)

# ADD
## SNPs measured in stage2 and stage1
replace <- stage2$SNP[stage2$SNP %in% stage1$SNP]
## SNPS only measured in stage2
add <- stage2[!stage2$SNP %in% stage1$SNP, ]
## SNPS only measured in 3A or B
add_meta <- meta[!meta$SNP %in% stage1$SNP, ]
## add SNPS only measured in stage 2 / 3A / 3B to stage 1
stage1 <- bind_rows(stage1, add, add_meta)
## SNPs measured in stage1+2 and 3A or B
replace_meta <- meta$SNP[meta$SNP %in% stage1$SNP]

# SUBSET
rownames(stage1) <- NULL # make sure no rownames
## rownames = SNPid - easier subsetting
stage1 <- stage1 %>% column_to_rownames("SNP")
stage2 <- stage2 %>% column_to_rownames("SNP")
meta <- meta %>% column_to_rownames("SNP")
## subset results for SNPs measured in stage 2
stage1[replace, colnames(stage2)] <- stage2[replace, colnames(stage2)]
## subset results for SNPs measured in stage 3A or B
stage1[replace_meta, colnames(meta)] <- meta[replace_meta, colnames(meta)]
## convert rownames back to column
stage1 <- stage1 %>% rownames_to_column("SNP") %>% mutate(trait="CN vs. AD")

# WRITE
write.table(
  stage1,
  file = "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/QTLs/traitQTL/IGAP/IGAP_19_AD.csv",
  sep = ",",
  col.names = T,
  row.names = F
)

#####################