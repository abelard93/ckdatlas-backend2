## SETUP ----------------------------------------
# load packages
library(getPass)
library(log4r)
library(plyr)
library(tidyverse)
library(tictoc)
library(ssh.utils)
library(readxl)
library(reticulate)
library(data.table)

#I've installed the following packages because I can't install tidyverse!
#library(dplyr)
#library(tidyr)
#library(readr)
#library(purrr)
#library(tibble)
#library(stringr)
#library(ggplot2)

#install.packages('R.utils')
#install.packages("tidyverse")
#install.packages("devtools")
#library(devtools); devtools::install_github("collectivemedia/ssh.utils")
#install.packages("ssh.utils", repos = "https://github.com/collectivemedia/ssh.utils")



# helper functions
zap <- function(){lst<-ls(envir=.GlobalEnv); lst<-lst[!(lst %in% c("zap","codes.makepath","data.makepath","store","restore","debugstore"))]; rm(list=lst, envir=.GlobalEnv) }
#codes.makepath <- function(p){file.path("~/Dokumente/adatlas-backend",p)}
#data.makepath <- function(p){file.path("/local/neo4j/ADatlas/import",p)}
codes.makepath <- function(p){file.path("~/Dokumente/ckdatlas-backend",p)}
data.makepath <- function(p){file.path("~/local/neo4j/ADatlas/import",p)}


# clear environment for clean start
zap()

# log time 
tic("Duration time of setup:")

# paths
data <- data.makepath("data/") # path to data repository data/
code <- codes.makepath("code/") # path to git "code" repository code/
log_path <- codes.makepath("log/")

#Current working directory: 
getwd()
#setwd("~/Dokumente/adatlas-backend")
# set working directory
setwd(data)
#setwd("~/Dokumente/adatlas-backend")
# setup python
use_python(python = "/usr/bin/python3",required = T)# 
#use_virtualenv("atlas") # environment where neo4j module is installed
#uri = "neo4j://sysmet-int-03:7687" # address of neo4j instance # 
uri = "bolt://10.0.0.114:7687" # address of neo4j instance # 
neo4j <- import(module = "neo4j",as = "neo4j")


### Connect to DB ----------------------------------------
driver <- neo4j$GraphDatabase$driver(uri, 
                                     auth=
                                       neo4j$basic_auth(
                                         "neo4j",
                                         getPass(msg ="Please enter password for the Neo4j DB:") 
                                         )
                                     )


### Load scripts ----------------------------------------
for (nm in list.files(codes.makepath("code/setup/toLoad/"))) {
  source(file.path(codes.makepath("code/setup/toLoad"), nm))
}


### Create log ----------------------------------------
log <- create_log(paste0(log_path, format(Sys.time(), "%Y%m%d_%H%M"), 
                         "_neo4j.log"))
err <- paste0(log_path, format(Sys.time(), "%Y%m%d_%H%M"), ".log")
file.create(err)
msg = paste("Timestamp",
            "File/Node/Type",
            "ERROR_code",
            "Info [free text]",
            sep = "\t")
cat(msg, file = err, sep = "\n", append = T)

### Set constraints/index ----------------------------------------
## constraints of uniqueness for node labels
nodes = c(
  "source",
  "organism",
  "cohort",
  "publication",
  "sampleType",
  "SMPDB",
  "metHMDB",
  "accHMDB",
  "metChEBI",
  "accChEBI",
  "metMetabolon",
  "gene",
  "transcript",
  "protein",
  "metOntology",
  "SNP",
  "metBiGG",
  "geBiGG",
  "reBiGG",
  "compBiGG",
  "trait",
  "metaTrait",
  "proteinCoding",
  "metSmiles",
  "metINCHI",
  "metIUPAC",
  "metFormula",
  "metName",
  "metKEGG",
  "metPubChem",
  "metSubPathway",
  "metSuperPathway",
  "metClassification"
)

add_constr(nodes, "sid")

# show constraints
results <-runQuery("SHOW CONSTRAINTS", read = TRUE)
#results$type

# index without constraint for faster search in functions
# new syntax for creating index: CREATE INDEX FOR (n:Label) ON (n.property)`
runQuery("CREATE INDEX ON :metMetabolon(COMP_ID)")
runQuery("CREATE INDEX ON :geneSymbol(sid)")
runQuery("CREATE INDEX ON :uniprot(sid)")
runQuery("CREATE INDEX ON :proteinCoding(cutoff)")
runQuery("CREATE INDEX ON :proteinCoding(ensembl_id)")
runQuery("CREATE INDEX ON :traitQTL(pvalue)")
runQuery("CREATE INDEX ON :mQTL(pvalue)")

# show indexes
#results <-runQuery("call db.indexes", read = TRUE)

## FILL DB ----------------------------------------

### Sources/traits ----------------------------------------
# additional info on sources
addSources()

# additional info on traits
addTraits()

### Knowledge DBs ----------------------------------------
# SMPDB
addSMPDB()
# HMDB
addHMDB()
# metChEBI
addChEBI()

### Metabolites ----------------------------------------
# metMetabolon
addMetabolon()
# biocrates metabolites
addBiocrates()
# CHECK metMEASURED, METABOLITE, HMDB MAPPING
checkMapping(err)
# GGM
addGGM()


### SNP - Gene - transcript - Protein ----------------------------------------

# genes
addGene()
# gene-wise cutoff
addSNiPACutoff()
# transcripts
addTranscript()
# proteins
addProtein()
# gene2transcript2protein
mappingGenes()
# SNPs
addSNP()

## BiGG models ----------------------------------------

addBiGG()

## gene networks ----------------------------------------

# differential expressed genes
addDEG()
# differential expressed proteins
addDEP()
# coexpression networks
addCOEXPR()

## QTLs ----------------------------------------

# expressionQTL (GTEx)
addeQTL()
# add traitQTLs (GWAS)
addTraitQTL()
# add mQTLs (population studies)
addmQTL()
# add mWAS 
addmWAS()


## gene symbols ----------------------------------------

# if genesymbols important for integration execute this before
addGeneSymbol()


## merge ----------------------------------------

mergeGenes()
mergeProteins()

# coabundance networks
addCOABUN()

# biodomains
addBioDomains()

## protein-coding? ----------------------------------------

query <- "MATCH (g:gene)-[:CODES]->(:transcript)-[:CODES]->(:protein) where g.chromosome IN ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT'] SET g :proteinCoding"
runQuery(query)

"MATCH (n:proteinCoding) where apoc.meta.type(n.sid)='STRING' set n.sid=[n.sid]" %>% runQuery()

query <- "MATCH (t:trait) set t.sid=[t.sid]"
runQuery(query)

query <- "MATCH (t:metaTrait) set t.sid=[t.sid]"
runQuery(query)

## fixes
query <- "MATCH (g:gene) where apoc.meta.type(g.sid)='STRING' set g.sid = [g.sid]" 
runQuery(query)

query <- "MATCH (g:gene) where not exists(g.ensembl_id) set g.ensembl_id = g.sid[0]"
runQuery(query)

"MATCH (n:gene) where apoc.meta.type(n.symbol_id)='STRING' set n.symbol_id=[n.symbol_id]" %>% runQuery()

"MATCH (n:uniprot) where apoc.meta.type(n.uniprot_id)='STRING' set n.uniprot_id=[n.uniprot_id]" %>% runQuery()


# query <- "Call apoc.periodic.iterate(\"MATCH (g:proteinCoding)-[:HIT]-(e:eQTL)-[:HIT]-(s:SNP) return g,s\" , \"MERGE (g)-[:MAP]->(s)\", {parallel:true,batchSize:100,iterateList:true,concurrency:3,retries:4})"
# cypher(db,query)



# final log entry
info(log,capture.output(toc()))

########## end ###########################################################################


## LEGACY - check if needed?

# metWCM
# addWCM()#
# CHECK IF NEW METABOLON CREATED
#checkMappingGwas(db,err)


# ## OPTIONAL: BLACKLISTING + ADJUSTMENTS

# 
# #blacklisting
# doBlacklisting("blacklist/recon3D_maria.csv")
# 
# ## fixes
# for (fix in list.files("../code/setup/fixes")) {
#   source(file.path("../code/setup/fixes", fix))
# }
# 

# ## GENERATE  HASH KEY AND REPRESENTATION OF DATA SCHEMA

# 
#key <- hashAllFiles(c(data, code))
#date <- format(Sys.time(), "%Y-%m-%d")
# cypher(db,"CREATE (:hashkey {key:{key}, date:{date}})",date=date,key=key)
#makeSchema(db,key,date,git)
 
