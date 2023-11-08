#' Add SMPDB - The Small Molecule Pathway Database
#'
#' @description
#' Reads .csv file containing SMPDB db with following format:
#' * sid
#' * Name
#' * Subject
#' First creates a node for every entry/row and then creates link to source for every node

addSMPDB <- function(){
  
  for (file in getFiles("SMPDB")) {
    # read info on dataset
    source <-getInfo(file, c("sid", "file_name", "node_type", "version"), err)
    # skip file if no/incomplete info
    if (!all(is.na(source))) {
      tic(paste0("Added SMPDB file ",source["file_name"]))
      #check/add source node for dataset
      checkSource(source, log)
      
      # create node if not exists
      query <- paste("USING PERIODIC COMMIT 1000 
                       LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props 
                       MERGE (n:SMPDB {sid:props.sid}) 
                       SET n += props 
                       ",sep="")
      runQuery(query, param = TRUE, periodic = TRUE, param_list=list(file=file))
     
      # create link to source
      query <- paste("USING PERIODIC COMMIT 1000 
                     LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props 
                     MATCH (n:SMPDB {sid:props.sid}) 
                     MATCH (s:source {sid:$properties.source}) 
                     MERGE (n)-[:FROM]->(s) 
                     ",sep="")
      runQuery(query, param = TRUE, periodic = TRUE, param_list=list(file=file,source=source["sid"]$sid))
      info(log, capture.output(toc())) # log
    }
  }
}
