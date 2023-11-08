####
# Script to parse SniPa SNP 2 gene mapping 
# and bring into correct format for Neo4j
# last modified: 2019-05-21
####

zap()
require(dplyr)
require(tidyr)

### COMMANDLINE ARGUMENTS
args = commandArgs(trailingOnly = T)
file          <- args[1]
output        <- args[2]

# output     <- "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/SNP/snp_2_gene_map.csv"
# file       <- "/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/SNiPA SNP-2-Gene mapping/snp_2_gene_map.tabix.gz"

## READ ---------------------------
snp <- read.delim(file=file) 

## UNSERIALIZE PHP ARRAY ----------
snp$SNPS <- character(nrow(snp))
for(i in 1:nrow(snp)){
  snp$SNPS[i] <- str_extract_all(snp$SNPARRAY[i], "rs[0-9]*")
}

## FORMAT ------------------------
snp <- snp %>%  
  select(-SNPARRAY) %>% 
  transmute(sid = SNPS, gene = ID, gene_name = NAME, chr = map_chr(CHR, ~paste0("chr",.))) %>% 
  unnest 

## WRITE -------------------------
write.csv(file = output,snp,row.names = F,quote = T,na="")
