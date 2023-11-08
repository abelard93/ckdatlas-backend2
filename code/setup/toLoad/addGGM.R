#' Add GGMs - Gaussian Graphical Models   
#'
#' @description
#' Reads GGM edge list from 'data > GGM' were following columns are needed:
#' * a: metabolite name
#' * b: metabolite name
#' * pvalue
#' * partialCorr
#' * significant: TRUE/FALSE
#' * qvalue: optional
#' CAUTION: only implemented for biochemical metabolite names so far. 

addGGM <- function(){
  
  for(file in getFiles("GGM")){
    # read info on dataset
    source <- getInfo(file,c("sid","file_name","node_type","metaboliteID"),err)
    # skip file if no/incomplete info
    if (!all(is.na(source))) {
      tic(paste0("Added GGM ", source["file_name"]))
      # check if node for source is present, if not create
      checkSource(source, log)
      
      query <- paste("
                     LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS prop 
                     MATCH (:metName {sid:apoc.text.replace(toLower(prop.a), '[^A-Za-z0-9]', '')})-[:HAS_NAME]-(a:",source["metaboliteID"]$metaboliteID,")
                     MATCH (:metName {sid:apoc.text.replace(toLower(prop.b), '[^A-Za-z0-9]', '')})-[:HAS_NAME]-(b:",source["metaboliteID"]$metaboliteID,")
                     with distinct prop, apoc.coll.zip(collect(a),apoc.coll.fill('a', size(collect(a))))+apoc.coll.zip(collect(b),apoc.coll.fill('b', size(collect(b)))) as zip
                     MATCH (s:source {sid:$properties.source}) 
                     CREATE (g:edgeGGM {pvalue:toFloat(prop.pvalue),
                                        partialCorr:toFloat(prop.partialCorr),
                                        significant:prop.significant,
                                        sampleType:COALESCE(prop.sampleType, null),
                                        qvalue:COALESCE(toFloat(prop.qvalue),null)}) 
                     CREATE (g) -[:FROM]-> (s) with zip,g
                     unwind zip as pair
                     CALL apoc.get.nodes(pair[0]) yield node as n
                     MERGE (n)-[:GGM {metabo:pair[1]}]->(g)",sep="")
    
      runQuery(query, param=T, param_list=list(file=file,source=source["sid"]$sid))
      info(log, capture.output(toc()))
    }
  }
}
