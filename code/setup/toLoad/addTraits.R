#' Add traits as nodes to neo4j DB
#'
#' @description
#' Loop through files in 'data > trait' and pre-create trait nodes 
#' Files are .csv tables which must contain `sid`,unique source id, and `metaTrait` as row name

addTraits<- function(){
  
  for (file in getFiles("trait")){
    #skip file if no/incomplete info
    tic(paste0("Added trait file: ", basename(file)))
    # add as node to db (sid)
    query <- paste("USING PERIODIC COMMIT 1000
                   LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props
                   MERGE (n:trait {sid:props.sid})
                   SET n += props
                   with n, split(n.metaTrait,'|') as mt
                   WHERE NOT '' IN mt
                   UNWIND mt as trait
                   MERGE (m:metaTrait {sid:trait})
                   MERGE (n)-[:PART_OF]->(m)
                   ",sep="")
    
    runQuery(query, param = TRUE, periodic = TRUE, param_list=list(file=file))
    info(log,capture.output(toc()))
  }
}