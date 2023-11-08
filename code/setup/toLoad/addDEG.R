## Differential expressed genes (DEG)
addDEG <- function() {
  for (file in getFiles("DEG")) {
    #read info on dataset
    source <- getInfo(file, c("sid", "file_name", "node_type"), err)
    #skip file if no/incomplete info
    if (!all(is.na(source))) {
      tic(paste0("Added DEG file ", source["file_name"]))
      #check/add source node for dataset
      checkSource(source, log)
      
      ##
      query <- paste0("CALL apoc.periodic.iterate(
                       \"LOAD CSV WITH HEADERS FROM 'file:///'+$properties AS props return distinct props.gene_id as g\",
                      \"MERGE (ge:gene:ensembl {sid:g,ensembl_id:g}) ON CREATE SET ge.source = 'deg'\", {batchSize:1000, parallel:true, iterateList:true, concurrency:8, params: {properties:$properties}})
                      ")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = file)
      ##
      query <- paste0("CALL apoc.periodic.iterate(
                       \"LOAD CSV WITH HEADERS FROM 'file:///'+$properties AS props return distinct props.trait as t\",
                      \"MERGE (:trait {sid:t}) \", {batchSize:1000, parallel:true, iterateList:true, concurrency:8, params: {properties:$properties}})
                      ")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = file)
      ## ADD edgeDEG NODE
      # sid is source_id
      query <- paste("
                     USING PERIODIC COMMIT 1000
                     LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props 
                     MATCH (g:gene:ensembl {sid:props.gene_id})
                     MATCH (s:source {sid:$properties.source}) 
                     MATCH (t:trait {sid:props.trait})
                     CREATE (n:edgeDEG)
                     SET n = props 
                     CREATE (g)-[:DEG]->(n)-[:DEG]->(t)
                     CREATE (n)-[:FROM]->(s)
                     ",sep="")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=file,source=source["sid"]$sid))

      info(log, capture.output(toc()))
 
    }
  }
}
