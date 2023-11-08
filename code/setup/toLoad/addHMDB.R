#' Add HMDB - The Human Metabolome Database 
#'
#' @description
#' Reads .csv file in 'data > HMDB' containing HMDB db and do:
#' * create hmdb node
#' * connect node to source
#' * create SA and PA mapping/nodes (PA: primary accession, SA: secondary accession)
#' * add hmdb to smpdb mapping

addHMDB <- function(){

    for (file in getFiles("metHMDB")) {
      # read info on dataset
      source <- getInfo(file, c("sid", "file_name", "node_type", "version"), err)
      # skip file if no/incomplete info
      if (!all(is.na(source))) {
        # check/add source node for dataset
        checkSource(source, log)
        tic(paste0("Added HMDB file ", source["file_name"]))
        
        # create node if not exists 
        query <- paste("
                       USING PERIODIC COMMIT 1000
                       LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS props 
                       MERGE (n:metHMDB {sid:props.sid}) 
                       SET n += props 
                       ",sep="")
        runQuery(query, periodic=TRUE) # send query
        
        # connect nodes to source
        query <- paste("
                       USING PERIODIC COMMIT 1000
                       LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props 
                       MATCH (n:metHMDB {sid:props.sid}) 
                       MATCH (s:source {sid:$properties.source}) 
                       MERGE (n)-[:FROM]->(s)
                       ",sep="")
        runQuery(query, param = TRUE, periodic = TRUE, param_list=list(file=file,source=source["sid"]$sid)) # send query
        info(log, capture.output(toc()))
        
        # create SA and PA mapping/nodes
        tic(paste0("Added HMDB sa and pa for file ", source["file_name"]))
        query <- paste("
                       MATCH (n:metHMDB) 
                       with n.sid as pa
                       MERGE (p:accHMDB {sid:pa})
                       ",sep="")
        runQuery(query) # send query
        
        query <- paste("
                       MATCH (n:metHMDB) 
                       with n as hmdb, n.sid as pa
                       MATCH (p:accHMDB {sid:pa})
                       MERGE (hmdb)-[:PA]->(p)
                       ",sep="")
        runQuery(query) # send query
        
        query <- paste("
                       MATCH (n:metHMDB) 
                       with n as hmdb, split(n.sa,'|') as sa
                       UNWIND sa as secAcc
                       MERGE (s:accHMDB {sid:secAcc})
                       CREATE (hmdb)-[:SA] -> (s)
                       ",sep="")
        runQuery(query) # send query
        info(log, capture.output(toc()))
        
        # check which secondary accessions (SA) are mapped to multiple primary accessions (PA)
        checkMappingHmdb(err) 
        
        # add HMDB to SMPDB links
        tic(paste0("Added HMDB-SMPDB mapping for file ", source["file_name"]))
        query <- paste("
                       MATCH (h:metHMDB) 
                       with h as hmdb, split(h.pathways_smpdb,'|') as smpdb
                       UNWIND smpdb as id
                       MERGE (s:SMPDB {sid:id})
                       CREATE (hmdb)-[:ANNOTATE]->(s)
                       ",sep="")
        runQuery(query) # send query
        info(log, capture.output(toc()))
        }
    }
}
