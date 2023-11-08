library(tidyverse)
zap()

########
## GGM 
########

## GGM KRUMSIEK
# READ AND FORMAT
GGM <- read_tsv("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/Emory brain proteomics GGM network/GGM.edges.brain.qval_lt_0.05.tsv") %>% 
  transmute(a=node1,b=node2,partialCorr=pcor,pvalue=pval,qvalue=qval)


GGM <- GGM %>% 
  mutate(significant = ifelse(qvalue <= 0.05, TRUE, FALSE))
# WRITE
write.table(
  GGM,
  file = "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/COABUN/consensus_proteomics_emory.csv",
  sep = ",",
  col.names = T,
  row.names = F
)

GGM <- read.csv("COABUN/consensus_proteomics_emory.csv") %>% mutate(a_uniprot=str_replace(a,"[.].*",""),b_uniprot=str_replace(b,"[.].*",""))
write.table(
  GGM,
  file = "COABUN/consensus_proteomics_emory.csv",
  sep = ",",
  col.names = T,
  row.names = F
)

# ## FILE INFO
# info <- readxl::read_excel(
#   "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/GGM/GGM.xlsx",
#   sheet = "info",
#   col_names = F
# )
# info[1, 2] <- "krumsiek_2012_GGM.csv"
# info[4, 2] <- "krumsiek_2012_GGM"
# # WRITE
# write.table(
#   info,
#   file = "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/GGM/krumsiek_2012_GGM_info.csv",
#   sep = ",",
#   col.names = F,
#   row.names = F
# )