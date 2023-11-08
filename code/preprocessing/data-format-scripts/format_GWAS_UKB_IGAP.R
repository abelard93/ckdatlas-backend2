library(tidyverse)
zap()



#################################
# GWAS UKB_IGAP meta analysis
#################################

# READ AND FORMAT
# parental meta
gwas <- data.table::fread("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/GWAS/3_UKB_AD_parental_meta_summary_output_June2019.txt") %>% 
  transmute(
    SNP = SNP,
    chromosome = CHR,
    bp = BP,
    effect_direction = DIR,
    pvalue = P,
    effect_allele = A1,
    trait = "AD_by_proxy"
  )
# WRITE
write.table(gwas,file = "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/QTLs/traitQTL/UKBiobank/UKB_AD_parental_meta.csv",sep = ",",col.names = T,row.names = F)

# READ AND FORMAT
# parental and AD meta
gwas2 <- data.table::fread("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/GWAS/4_UKB_IGAP_AD_meta_summary_output_June2019.txt") %>% 
  transmute(
    SNP = SNP,
    chromosome = CHR,
    bp = BP,
    effect_direction = DIR,
    pvalue = P,
    effect_allele = A1,
    trait = "AD_by_proxy"
  )
# WRITE 
write.table(gwas2,file = "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/QTLs/traitQTL/UKBiobank/UKB_IGAP_AD_meta.csv",sep = ",",col.names = T,row.names = F)

#####################