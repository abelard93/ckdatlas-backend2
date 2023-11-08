library(tidyverse)
zap()

########
## GGM 
########

## GGM KRUMSIEK
# READ AND FORMAT
GGM <- readxl::read_excel("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/GGM_krumsiek/GGM.xlsx",
                          sheet = "GGM list",
                          skip = 1)
colnames(GGM) <- c("a", "b", "partialCorr", "pvalue")
GGM <- GGM %>% 
  mutate(significant = ifelse(abs(partialCorr) >= 0.1603 & pvalue <= 7.96e-7, TRUE, FALSE))
# WRITE
write.table(
  GGM,
  file = "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/GGM/krumsiek_2012_GGM.csv",
  sep = ",",
  col.names = T,
  row.names = F
)
## FILE INFO
info <- readxl::read_excel(
  "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/GGM/GGM.xlsx",
  sheet = "info",
  col_names = F
)
info[1, 2] <- "krumsiek_2012_GGM.csv"
info[4, 2] <- "krumsiek_2012_GGM"
# WRITE
write.table(
  info,
  file = "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/GGM/krumsiek_2012_GGM_info.csv",
  sep = ",",
  col.names = F,
  row.names = F
)
