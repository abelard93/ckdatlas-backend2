# quality control


library(RNeo4j)
library(log4r)
library(plyr)
library(tidyverse)
library(tictoc)
library(ssh.utils)
library(readxl)

## -------------- setup ---------------

zap <- function(){lst<-ls(envir=.GlobalEnv); lst<-lst[!(lst %in% c("zap","codes.makepath","data.makepath","store","restore","debugstore"))]; rm(list=lst, envir=.GlobalEnv) }

codes.makepath <- function(p){file.path("/storage/metabostore/project/ADatlas/repo",p)}
data.makepath <- function(p){file.path("/storage/metabostore/project/ADatlas/repo",p)}


log_path <- codes.makepath("log/")
data <- data.makepath("data/") 

db_addr <- "http://metabo-interactive:7474/db/data"
db = startGraph(db_addr, username = "neo4j", password = "sysdiab")

## SET WORKING DIRECTORY
setwd(data)

## LOAD SCRIPTS
for (nm in list.files(codes.makepath("code/setup/toLoad/"))) {
  source(file.path(codes.makepath("code/setup/toLoad"), nm))
}


## -------------- log file ---------------

log_file <- paste0(log_path, format(Sys.time(), "%Y%m%d_%H%M"), "_quality_control.log")
file.create(log_file)

writeLog <- function(...){
 msg <- paste(...,sep = "\t")
 message(msg)
 cat(msg, file = log_file, sep = "\n", append = T)
}

writeHeader <- function(string){
  msg = paste0("#---------- ",string," ---------- \n")
  message(msg)
  cat(msg, file = log_file, sep = "\n", append = T)
}

## --------------  header (ToDo) -----------


## -------------- eQTL ---------------

# logfile
writeHeader("eQTL")
writeLog("Type", "#Input", "#DB")

for(i in c("cis","trans")){
  
  countQ <- 0 # eQTL 
  countS <- c() # SNPs
  countG <- c() # genes
  
  for (file in getFiles(paste0("QTLs/eQTL/",i))) {
    
    df <- data.table::fread(file)
    countQ <- countQ + nrow(df) # count number of eQTL (edges)
    countS <- unique(c(countS,df$rs_id_dbSNP151_GRCh38p7)) # collect unique SNPs
    countG <- unique(c(countG,df$gene_id)) # collect unique genes
    
    message(paste0("Reading done: ",file))
    
  }
  
  ## number of eQTL nodes
  true <- paste0("match (e:eQTL:",i,") return count(e)") %>% cypher(db,.)
  writeLog(paste0("eQTL_",i),countQ, true)
  
  ## number of tissue
  true <- paste0("match (e:eQTL:",i,")-[:FROM]-(s:source) return count(distinct s.sample_type)") %>% cypher(db,.)
  writeLog(paste0("tissue_",i), length(getFiles(paste0("QTLs/eQTL/",i))), true)
  
  
  ## number of SNPs
  true <- paste0("match (e:eQTL:",i,")<-[:HIT]-(s:SNP) return count(distinct s)") %>% cypher(db,.)
  writeLog(paste0("SNPs_",i), length(countS), true)
  
  ## number of genes
  true <- paste0("match (e:eQTL:",i,")-[:HIT]->(s:gene) return count(distinct s)") %>% cypher(db,.)
  writeLog(paste0("genes_",i), length(countG), true)
  
}

# ## number of sources
# true <- paste0("match (s:source {node_type:'eQTL'}) return count(distinct s)") %>% cypher(db,.)
# countF <- length(list.files("QTLs/eQTL",recursive = T,pattern = "_info.csv"))
# writeLog(paste0("source",i), countF, true)
writeLog("\n")
## add node type to source!!!

## -------------- mWAS ---------------

# mWAS
writeHeader("mWAS")
writeLog("Type", "#Input", "#DB")

countQ <- 0 # mWAS
countT <- c() # traits
countM <- c() # metabolites


for (file in getFiles("QTLs/mWAS/")) {
  
  df <- data.table::fread(file)
  countQ <- countQ + nrow(df) # count number of eQTL (edges)
  countT <- unique(c(countT,df$trait)) # collect unique SNPs
  countM <- unique(c(countM,df$met)) # collect unique genes
  
  message(paste0("Reading done: ",file))
  
}

## number of mWAS nodes
true <- "match (e:mWAS) return count(e)" %>% cypher(db,.)
writeLog("mWAS",countQ, true)

## number of SNPs
true <- paste0("match (e:mWAS)-[:HIT]->(s:trait) return count(distinct s)") %>% cypher(db,.)
writeLog("traits", length(countT), true)

## number of genes
true <- paste0("match (e:mWAS)<-[:HIT]-(s:metMeasured) return count(distinct s)") %>% cypher(db,.)
writeLog("metabolites", length(countM), true)

## number of genes
true <- paste0("match (e:mWAS)<-[:HIT]-(s:metMetabolon) return count(distinct s)") %>% cypher(db,.)
writeLog("metabolon", "-", true)

## number of genes
true <- paste0("match (e:mWAS)<-[:HIT]-(s:metBiocrates) return count(distinct s)") %>% cypher(db,.)
writeLog("biocrates", "-", true)

writeLog("\n")



## -------------- mQTL ---------------

writeHeader("mQTL")
writeLog("Type", "#Input", "#DB")

countQ <- 0 # mWAS
countS <- c() # SNPs
countM <- c() # metabolites


for (file in getFiles("QTLs/mQTL")) {
  
  df <- data.table::fread(file)
  countQ <- countQ + nrow(df) # count number of eQTL (edges)
  countS <- unique(c(countS,df$rsID)) # collect unique SNPs
  countM <- unique(c(countM,df$metabolite)) # collect unique genes
  
  message(paste0("Reading done: ",file))
  
}

## number of mWAS nodes
true <- "match (e:mQTL) return count(e)" %>% cypher(db,.)
writeLog("mQTL",countQ, true)

## number of mWAS nodes
true <- "match (e:mQTL) where not (e)-[:HIT]-(:metMeasured) return count(e)" %>% cypher(db,.)
writeLog("isolated mQTL","-", true)

## number of SNPs
true <- paste0("match (e:mQTL)<-[:HIT]-(s:SNP) return count(distinct s)") %>% cypher(db,.)
writeLog("SNPs", length(countS), true)

## number of metabolites
true <- paste0("match (e:mQTL)-[:HIT]->(s:metMeasured) return count(distinct s)") %>% cypher(db,.)
writeLog("metabolites", length(countM), true)

## number of metabolites
true <- paste0("match (e:mQTL)-[:HIT]->(s:metMetabolon) return count(distinct s)") %>% cypher(db,.)
writeLog("metabolon", "-", true)

## number of metabolites
true <- paste0("match (e:mQTL)-[:HIT]->(s:metBiocrates) return count(distinct s)") %>% cypher(db,.)
writeLog("biocrates", "-", true)

writeLog("\n")



## -------------- traitQTL ---------------

writeHeader("traitQTL")
writeLog("Type", "#Input", "#DB")

countQ <- 0 # traitQTL
countS <- c() # SNPs
countM <- c() # traits


for (file in getFiles("QTLs/traitQTL")) {
  
  df <- data.table::fread(file)
  countQ <- countQ + nrow(df) # count number of eQTL (edges)
  countS <- unique(c(countS,df$SNP)) # collect unique SNPs
  countM <- unique(c(countM,df$trait)) # collect unique traits
  
  message(paste0("Reading done: ",file))
  
}

## number of mWAS nodes
true <- "match (e:traitQTL) return count(e)" %>% cypher(db,.)
writeLog("traitQTL",countQ, true)

## number of SNPs
true <- paste0("match (e:traitQTL)<-[:HIT]-(s:SNP) return count(distinct s)") %>% cypher(db,.)
writeLog("SNPs", length(countS), true)

## number of traits
true <- paste0("match (e:traitQTL)-[:HIT]->(s:trait) return count(distinct s)") %>% cypher(db,.)
writeLog("traits", length(countM), true)



writeLog("\n")

## -------------- DEG ---------------

writeHeader("DEG")
writeLog("Type", "#Input", "#DB")

countQ <- 0 # DEG
countS <- c() # gene
countM <- c() # traits


for (file in getFiles("DEG/")) {
  
  df <- data.table::fread(file)
  countQ <- countQ + nrow(df) # count number of eQTL (edges)
  countS <- unique(c(countS,df$gene_id)) # collect unique SNPs
  countM <- unique(c(countM,df$trait)) # collect unique traits
  
  message(paste0("Reading done: ",file))
  
}

## number of mWAS nodes
true <- "match (e:edgeDEG) return count(e)" %>% cypher(db,.)
writeLog("DEG",countQ, true)

## number of SNPs
true <- paste0("match (e:edgeDEG)<-[:DEG]-(s:gene) return count(distinct s)") %>% cypher(db,.)
writeLog("genes", length(countS), true)

## number of traits
true <- paste0("match (e:edgeDEG)-[:DEG]->(s:trait) return count(distinct s)") %>% cypher(db,.)
writeLog("traits", length(countM), true)



writeLog("\n")

## -------------- DEP ---------------

writeHeader("DEP")
writeLog("Type", "#Input", "#DB")

countQ <- 0 # DEG
countS <- c() # proteins
countM <- c() # traits


for (file in getFiles("DEP")) {
  
  df <- data.table::fread(file)
  countQ <- countQ + nrow(df) # count number of eQTL (edges)
  countS <- unique(c(countS,df$uniprot_id)) # collect unique SNPs
  countM <- unique(c(countM,df$trait)) # collect unique traits
  
  message(paste0("Reading done: ",file))
  
}

## number of mWAS nodes
true <- "match (e:edgeDEP) return count(e)" %>% cypher(db,.)
writeLog("DEP",countQ, true)

## number of uniprots not unique => higher number
true <- paste0("match (e:edgeDEP)<-[:DEP]-(s:uniprot) return count(distinct s)") %>% cypher(db,.)
writeLog("proteins", length(countS), true)

## number of traits
true <- paste0("match (e:edgeDEP)-[:DEP]->(s:trait) return count(distinct s)") %>% cypher(db,.)
writeLog("traits", length(countM), true)



writeLog("\n")


## -------------- COEXPR ---------------

writeHeader("COEXPR")
writeLog("Type", "#Input", "#DB")

countQ <- 0 # COEXPR
countS <- c() # genes



for (file in getFiles("COEXPR")) {
  
  df <- data.table::fread(file)
  countQ <- countQ + nrow(df) # count number of eQTL (edges)
  countS <- unique(c(countS,df$gene_a,df$gene_b)) # collect unique SNPs
 # countM <- unique(c(countM,df$trait)) # collect unique traits
  
  message(paste0("Reading done: ",file))
  
}

## number of mWAS nodes
true <- "match (e:edgeCOEXPR) return count(e)" %>% cypher(db,.)
writeLog("COEXPR",countQ, true)

## genes
true <- paste0("match (e:edgeCOEXPR)-[:COEXPR]-(s:gene) return count(distinct s)") %>% cypher(db,.)
writeLog("genes", length(countS), true)

# ## number of traits
# true <- paste0("match (e:edgeDEP)-[:DEP]->(s:trait) return count(distinct s)") %>% cypher(db,.)
# writeLog("traits", length(countM), true)



writeLog("\n")

## -------------- COABUN ---------------

writeHeader("COABUN")
writeLog("Type", "#Input", "#DB")

countQ <- 0 # COABUN
countS <- c() # proteins



for (file in getFiles("COABUN")) {
  
  df <- data.table::fread(file)
  countQ <- countQ + nrow(df) # count number of eQTL (edges)
  countS <- unique(c(countS,df$a_uniprot,df$b_uniprot)) # collect unique SNPs
  # countM <- unique(c(countM,df$trait)) # collect unique traits
  
  message(paste0("Reading done: ",file))
  
}

## number of mWAS nodes
true <- "match (e:edgeCOABUN) return count(e)" %>% cypher(db,.)
writeLog("COABUN",countQ, true)

## uniprot not unique
true <- paste0("match (e:edgeCOABUN)-[:COABUN]-(s:uniprot) return count(distinct s)") %>% cypher(db,.)
writeLog("proteins", length(countS), true)

# ## number of traits
# true <- paste0("match (e:edgeDEP)-[:DEP]->(s:trait) return count(distinct s)") %>% cypher(db,.)
# writeLog("traits", length(countM), true)



writeLog("\n")
