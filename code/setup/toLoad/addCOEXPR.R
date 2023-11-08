## Co-expression Network
addCOEXPR <- function() {
  for (file in getFiles("COEXPR")) {
    #read info on dataset
    source <- getInfo(file, c("sid", "file_name", "node_type"), err)
    #skip file if no/incomplete info
    if (!all(is.na(source))) {
      tic(paste0("Added COEXPR file ", source["file_name"]))
      #check/add source node for dataset
      checkSource(source, log)
      
      ##
      query <- paste0("CALL apoc.periodic.iterate(
                       \"LOAD CSV WITH HEADERS FROM 'file:///'+$properties AS props return distinct props.gene_a as g\",
                      \"MERGE (ge:gene:ensembl {sid:g,ensembl_id:g}) ON CREATE SET ge.source = 'coexpr' \", {batchSize:1000, parallel:true, iterateList:true, concurrency:8, params: {properties:$properties}})
                      ")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = file)
      
      ##
      query <- paste0("CALL apoc.periodic.iterate(
                       \"LOAD CSV WITH HEADERS FROM 'file:///'+$properties AS props return distinct props.gene_b as g\",
                      \"MERGE (ge:gene:ensembl {sid:g,ensembl_id:g}) ON CREATE SET ge.source = 'coexpr' \", {batchSize:1000, parallel:true, iterateList:true, concurrency:8, params: {properties:$properties}})
                      ")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = file)
      ## ADD edgeCOEXPR NODE
      # sid is source_id
      query <- paste("
                     USING PERIODIC COMMIT 1000
                     LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props 
                     MATCH (g1:gene:ensembl {sid:props.gene_a})
                     MATCH (g2:gene:ensembl {sid:props.gene_b})
                     MATCH (s:source {sid:$properties.source})
                     CREATE (n:edgeCOEXPR {gene_a:props.gene_a, gene_b:props.gene_b, sample_type:COALESCE(props.sampleType,null)})
                     CREATE (g1)-[:COEXPR]->(n)<-[:COEXPR]-(g2)
                     CREATE (n)-[:FROM]->(s)
                     ",sep="")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=file,source=source["sid"]$sid))
      info(log, capture.output(toc()))
      
      
      }
    }
}
