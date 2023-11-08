library(tidyverse)
zap()

###########################
## SNPS_per_gene annotation
###########################

# READ AND FORMAT
snps_per_gene <-
  read.delim("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/numsnps_per_gene.tsv") %>%
  transmute(
    sid = ID,
    NUMSNPS_SNiPA = NUMSNPS,
    cutoff = 0.05 / NUMSNPS,
    cutoff_info = "Gene-wise Bonferroni; 0.05/NUMSNPS_SNiPA"
  )
# WRITE
write.csv(
  snps_per_gene,
  "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/gene/numsnps_per_gene_cutoff.csv",
  row.names = F,
  quote = T,
  na = ""
)



#####################