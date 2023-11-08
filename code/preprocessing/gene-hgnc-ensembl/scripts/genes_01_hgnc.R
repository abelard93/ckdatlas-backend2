require(tidyverse)
library(tidyr)
## requires development version of tidy r 
## function unite needs option na.rm
##
##
## PARAMS
args = commandArgs(trailingOnly = T)
file.in  <- args[1]
file.out <- args[2]
##
## file.in <- "../raw/genes/hgnc_complete_set.txt"
##
################################################################################
## READ
################################################################################
data.hgnc <- read_delim(file.in, delim = "\t")# colClasses = list(character = c(32,48,41,33,39,35)))

## REMOVE WITHDRAWN
data.hgnc <- data.hgnc %>%
  filter(name != "entry withdrawn")

## ONLY UPPER CASE
data.hgnc <- data.hgnc %>%
  mutate_at(vars(symbol, prev_symbol, alias_symbol),toupper)

## ONLY ENTRIES WITH ENSEMBLE OR ENTREZ
data.hgnc <- data.hgnc %>% 
  filter(!is.na(ensembl_gene_id)|!is.na(entrez_id))

## SYNONYMS
data.hgnc <- data.hgnc %>% 
  tidyr::unite(col = symbol_all,symbol,prev_symbol,alias_symbol,sep = "|", na.rm=T) %>% 
  transmute(symbol_all, ensembl = ensembl_gene_id, entrez = entrez_id) %>% 
  unnest(symbol = strsplit(symbol_all, "\\|")) %>% 
  select(-symbol_all) %>% gather(source,sid,ensembl:entrez) %>% unique() %>% na.omit

data.entrez <- data.hgnc %>% filter(source=="entrez")
data.ensembl <- data.hgnc %>% filter(source=="ensembl")

################################################################################
## WRITE
################################################################################
write.csv(x = data.entrez, file = paste0(file.out,"_entrez.csv"), quote = T,na="",row.names = F)
write.csv(x = data.ensembl, file = paste0(file.out,"_ensembl.csv"), quote = T,na="",row.names = F)
