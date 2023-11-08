## Differential expressed proteins (DEP)
addDEP <- function() {
  for (file in getFiles("DEP")) {
    #read info on dataset
    source <- getInfo(file, c("sid", "file_name", "node_type"), err)
    #skip file if no/incomplete info
    if (!all(is.na(source))) {
      tic(paste0("Added DEP file ", source["file_name"]))
      #check/add source node for dataset
      checkSource(source, log)
      
      ## ADD edgeDEG NODE
      # sid is source_id
      query <- paste("
                     USING PERIODIC COMMIT 100
                     LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props 
                     CREATE (n:edgeDEP)
                     SET n = props  with n, props
                     MATCH (s:source {sid:$properties.source}) 
                     MERGE (t:trait {sid:props.trait})
                     CREATE (n)-[:DEP]->(t)
                     CREATE (n)-[:FROM]->(s) with n, props
                     MATCH (g:uniprot {sid:props.uniprot_id}) 
                     CREATE (g)-[:DEP]->(n)
                     ",sep="")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=file,source=source["sid"]$sid))
      info(log, capture.output(toc()))
      
      
      
    }
  }
}