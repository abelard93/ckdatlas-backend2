library(getPass)
library(log4r)
library(plyr)
library(tidyverse)
library(tictoc)
library(ssh.utils)
library(readxl)
library(data.table)
library(reticulate)
library(data.table)

zap <- function(){lst<-ls(envir=.GlobalEnv); lst<-lst[!(lst %in% c("zap","codes.makepath","data.makepath","store","restore","debugstore"))]; rm(list=lst, envir=.GlobalEnv) }

#codes.makepath <- function(p){file.path("/storage/metabostore/project/ADatlas/repo",p)}
#data.makepath <- function(p){file.path("/storage/metabostore/project/ADatlas/repo",p)}
codes.makepath <- function(p){file.path("~/Dokumente/ckdatlas-backend",p)}
data.makepath <- function(p){file.path("~/local/neo4j/ADatlas/import",p)}

tic("Duration time of ATLAS setup:")

## PATHS #####################################################################################

data <- data.makepath("data/") # path to data repository 
code <- codes.makepath("code/") # path to git "code" repository 
log_path <- codes.makepath("log/")

## SETUP PYTHON #################################

#use_python(python = "/home/icb/maria.woerheide/anaconda3/envs/atlas/bin/python",required = T)
use_python(python = "/usr/bin/python3",required = T)
#use_virtualenv("atlas")
#uri = "neo4j://metabo-interactive:7687"
uri = "bolt://10.0.0.114:7687" # address of neo4j instance

neo4j <- import(module = "neo4j",as = "neo4j")
#token <- neo4j$basic_auth("neo4j",getPass(msg ="Please enter password for the Neo4j DB:"))
driver <- neo4j$GraphDatabase$driver(uri, 
                                     auth=
                                       neo4j$basic_auth(
                                         "neo4j",
                                         getPass(msg ="Please enter password for the Neo4j DB:")
                                       )
)

## SETUP #####################################################################################

## SET WORKING DIRECTORY
setwd(data)


## LOAD SCRIPTS
for (nm in list.files(codes.makepath("code/setup/toLoad/"))) {
  source(file.path(codes.makepath("code/setup/toLoad"), nm))
}

"CALL apoc.warmup.run(true,true,true)" %>% runQuery(read = TRUE)

for (nm in list.files(codes.makepath("code/setup_CKDatlas/toLoad/"))) {
  source(file.path(codes.makepath("code/setup_CKDatlas/toLoad"), nm))
}

## CREATE LOG FILES
log <- create_log(paste0(log_path, format(Sys.time(), "%Y%m%d_%H%M"), "_adatlas.log"))
err <- paste0(log_path, format(Sys.time(), "%Y%m%d_%H%M"), "_adatlas_err.log")
file.create(err)
msg = paste("Timestamp",
            "File/Node/Type",
            "ERROR_code",
            "Info [free text]",
            sep = "\t")
cat(msg, file = err, sep = "\n", append = T)


##



## CREATE META METABOLITE NODES      ----------------------------

addMetaMetabolite()

## EXTRACT AND ADD EXTRACTION LAYER  ----------------------------

getSignRels()
addSignRels()




## EXTRACT AND ADD EXTRACTION LAYER  ----------------------------
("MATCH (a:metMeta) where (a)-[:GENETIC_ASSOCIATION]-() set a :atlas") %>% runQuery()

("MATCH (a:proteinCoding) where (a)-[:GENETIC_ASSOCIATION]-() set a :atlas") %>% runQuery()

("MATCH (a:metMeta) where (a)-[:METABOLIC_ASSOCIATION]-() set a :atlas") %>% runQuery()

("MATCH (a:metMeta) where (a)-[:PARTIAL_CORRELATION]-() set a :atlas") %>% runQuery()

#("MATCH (a:proteinCoding) where (a)-[:COEXPRESSION]-() set a :atlas") %>% runQuery()

("MATCH (a:proteinCoding) where (a)-[:COREGULATION]-() set a :atlas") %>% runQuery()

("MATCH (a:proteinCoding) where (a)-[:COABUNDANCE]-() set a :atlas") %>% runQuery()

("MATCH (a:proteinCoding) where (a)-[:DEG]-() set a :atlas") %>% runQuery()

("MATCH (a:proteinCoding) where (a)-[:CODES]-(:transcript)-[:CODES]-(:protein)-[:DEP]-() set a :atlas") %>% runQuery()

# final log entry
info(log,capture.output(toc()))
