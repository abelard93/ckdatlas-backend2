
## Co-expression Network
addCOABUN <- function() {
  for (file in getFiles("COABUN")) {
    #read info on dataset
    source <- getInfo(file, c("sid", "file_name", "node_type"), err)
    #skip file if no/incomplete info
    if (!all(is.na(source))) {
      
      tic(paste0("Added COABUN file ", source["file_name"]))
      #check/add source node for dataset
      checkSource(source, log)
  
      query <- paste0("CALL apoc.periodic.iterate(\"
                      LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS prop RETURN prop\",
                       \"with split(prop.a_uniprot,'|') as uni,split(prop.b_uniprot,'|') as uni2, prop
                       MATCH (a:uniprot) where any(x IN a.sid where x IN uni)
                       MATCH (b:uniprot) where any(x IN b.sid where x IN uni2)
                       with distinct prop, apoc.coll.zip(collect(a),apoc.coll.fill('a', size(collect(a))))+apoc.coll.zip(collect(b),apoc.coll.fill('b', size(collect(b)))) as zip
                       MATCH (s:source {sid:'",source["sid"]$sid,"'})
                       CREATE (g:edgeCOABUN {
                                 pvalue:COALESCE(toFloat(prop.pvalue),null),
                                 partialCorr:toFloat(prop.partialCorr),
                                 significant:prop.significant,
                                 sample_type:COALESCE(prop.sample_type, null),
                                 qvalue:COALESCE(toFloat(prop.qvalue),null)})
                       CREATE (g) -[:FROM]-> (s) with zip,g
                       unwind zip as pair
                       CALL apoc.get.nodes(pair[0]) yield node as n
                       MERGE (n)-[:COABUN {protein:pair[1]}]->(g)\",{batchSize:10, parallel:false})")

      
      runQuery(query, periodic = TRUE)
    
      info(log, capture.output(toc()))
      
    }
  }
}
