####
# Script from Jonas Zierer to parse recon3D
# adapted by Maria Woerheide for integration into Neo4j
# last modified: 2019-04-16
####
zap()
require(tidyverse)
require(xml2)
require(stringdist)
## GRAPHS
require(igraph)
library(zoo)
## PARAMS
args = commandArgs(trailingOnly = T)
file.info          <- args[1]
out.path          <- args[2]

#file.info <- "~/Documents/PHD/neo4j/neo4jserver/code/BiGG/data/20190429_microbes_info.tsv"

## READ FILE 
microbe.info <- read_tsv(file.info)


################################################################################
## MAKE INFO FILES
################################################################################

# add filename
microbe.info <- microbe.info %>% mutate(organism.orig = organism) %>% 
  mutate(organism = map_chr(organism, ~unlist(str_split(.,"="))[[1]])) %>% 
  mutate(file_name = map_chr(organism, ~str_replace_all(.,"[[:punct:]]","_"))) %>% 
  mutate(file_name = map_chr(file_name, ~paste0(str_replace_all(.," ","_"),".xml"))) %>%
  mutate(file_name = map_chr(file_name, ~str_replace_all(.,"(_+)","_"))) %>% 
  mutate(file_name = map_chr(file_name, ~str_replace_all(.,"(_.xml)",".xml"))) %>%
  mutate(sid = organism)

# check if all files existant
for(i in microbe.info$file_name){
  if(!file.exists(paste0("~/Documents/PHD/krumsieklab/neo4jserver/code/BiGG/data/raw/agora/",i))){
    print(paste0("Cannot find file: ",i))
  }
}

## write result csv
for(i in 1:nrow(microbe.info)){
  out <- microbe.info %>% 
    filter(row_number()==i) %>% 
    t %>% na.omit %>% data.frame 
  out.name <- microbe.info %>% 
    filter(row_number()==i) %>% .$file_name %>% str_replace(.,".xml", "_info.csv")
  write.table(out,file = paste0(out.path,out.name),sep = ",",col.names = F,row.names = T)
}

