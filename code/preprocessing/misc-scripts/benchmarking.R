## BENCHMARKING
zap()
library(RNeo4j)
library(log4r)
library(plyr)
library(tidyverse)
#library(dplyr)
library(tictoc)
library(ssh.utils)
library(readxl)
#library(neo4r)

data <- data.makepath("neo4jserver/data/") # path to data repository 
code <- codes.makepath("neo4jserver/code/") # path to git "code" repository 
log_path <- codes.makepath("neo4jserver/log/")

setwd(data)

## LOAD SCRIPTS
for (nm in list.files(codes.makepath("neo4jserver/code/setup/toLoad/"))) {
  source(file.path(codes.makepath("neo4jserver/code/setup/toLoad"), nm))
}

## logs
log <- create_log(paste0(log_path, format(Sys.time(), "%Y%m%d_%H%M"), "_neo4j.log"))
err <- paste0(log_path, format(Sys.time(), "%Y%m%d_%H%M"), ".log")

db_addr <- "http://vmetaneo4j:7474/db/data/"

db = startGraph(db_addr, username = "neo4j", password = "sysdiab")
clear(db,input = FALSE)

addConstraint(db, "SNP" , "sid")


tic()
for (file in getFiles("QTLs/mQTL")) {
  source <- getInfo(file, c("sid", "file_name"), err)
  checkSource(source, db, log)
  
  ############
  
  query <- paste0("CALL apoc.periodic.iterate(
                       \"LOAD CSV WITH HEADERS FROM {properties} AS props return distinct props.rsID as snp\",
                 \"MERGE (:SNP {sid:snp})\", {batchSize:2000, parallel:true, iterateList:true, concurrency:4, params: {properties:{properties}}})
                 ")
  system.time(cypher(db,query, properties=file))
  
  query <- paste0("CALL apoc.periodic.iterate(
                 \"unwind {properties} as prop LOAD CSV WITH HEADERS FROM 'file:///'+prop.file AS props return props, prop.source as su\",
                 \"CREATE (n:mQTL {pvalue:toFloat(props.pvalue), rsID:props.rsID, metabolite:apoc.text.replace(toLower(props.metabolite), '[^A-Za-z0-9]', '')}) with n,su MATCH (s:source {sid:su}) CREATE (n)-[:FROM]->(s) with n MATCH (snp:SNP {sid:n.rsID}) CREATE (snp)-[:HIT]->(n) with n MATCH (:metName {sid:n.metabolite})<-[:HAS_NAME]-(m:metMetabolon) CREATE (n)-[:HIT]->(m)\", {batchSize:2000, parallel:true, iterateList:true, concurrency:4, params: {properties:{properties}}})
                 ")
  system.time(cypher(db,query,properties = list(file=file,source=source["sid"]$sid)))
  
}
toc()


#############


library(neo4r)
con <- neo4j_api$new(
  url = "http://vmetaneo4j:7474", 
  user = "neo4j", 
  password = "sysdiab"
)
con$ping()

'CREATE CONSTRAINT ON (h:SNP) ASSERT h.sid IS UNIQUE' %>%
  call_neo4j(con)


system.time("CALL apoc.periodic.iterate('LOAD CSV WITH HEADERS FROM {props}.file AS props return distinct props.rsID as snp', 'MERGE (:SNP {sid:snp})', {batchSize:2000, parallel:true, iterateList:true, concurrency:4, params: {props:{props}}});" %>% neo4r::call_neo4j(con, params=list(props=list(file=file))))




