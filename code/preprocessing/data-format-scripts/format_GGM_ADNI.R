library(tidyverse)
zap()

########
## GGM 
########

## ADNI (2 GGMs)
# READ AND FORMAT (1)
ggm <- read.csv2(
  "/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/MWAS outcomes & GGMs/ADNI-1GO2-GGM.edges.bile-acids.tsv",
  sep = "\t",
  stringsAsFactors = F
) %>%
  transmute(
    a = node1,
    b = node2,
    partialCorr = pcor,
    pvalue = map_dbl(X.log10.pval.,  ~ 10 ^ (-1 * as.numeric(.x))), # transform -log10(pval) back
    qvalue = as.numeric(qval),
    significant = ifelse(qvalue <= 0.05, TRUE, FALSE)
  )
# WRITE
write.table(
  ggm,
  file = "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/GGM/ADNI_1GO2_bile_acids.csv",
  sep = ",",
  col.names = T,
  row.names = F
)
# READ AND FORMAT (2)
ggm <- read.csv2(
  "/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/MWAS outcomes & GGMs/ADNI-1GO2-GGM.edges.P180.tsv",
  sep = "\t",
  stringsAsFactors = F
) %>% 
  na.omit() %>%
  transmute(
    a = node1,
    b = node2,
    partialCorr = pcor,
    pvalue = map_dbl(X.log10.pval.,  ~ 10 ^ (-1 * as.numeric(.x))), # transform -log10(pval) back
    qvalue = as.numeric(qval),
    significant = ifelse(qvalue <= 0.05, TRUE, FALSE)
  )
# WRITE
write.table(
  ggm,
  file = "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/GGM/ADNI_1GO2_p180.csv",
  sep = ",",
  col.names = T,
  row.names = F
)

#####################