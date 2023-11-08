## metabolite trait associations form (mWAS)
addmWAS <- function() {
  for (file in getFiles("QTLs/mWAS")) {
    #read info on dataset
    source <- getInfo(file, c("sid", "file_name"), err)
    #skip file if no/incomplete info
    if (!all(is.na(source))) {
      tic(paste0("Added mWAS file ", source["file_name"]))
      #check/add source node for dataset
      checkSource(source, log)
      
      if(source["metaboliteID"]$metaboliteID =="COMP_ID"){
        
        query <- paste0("
                        USING PERIODIC COMMIT 500
                        LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS prop 
                        MATCH (m:metMetabolon {COMP_ID:prop.COMP_ID}) 
                        with collect(m) as metabo, prop
                        MATCH (s:source {sid:$properties.source}) 
                        MERGE (t:trait {sid:prop.trait})
                        CREATE (e:mWAS {pvalue:toFloat(prop.pvalue),
                                trait:COALESCE(prop.trait,null),
                                beta:COALESCE(prop.beta,null),
                                SE:COALESCE(prop.SE,null),
                                qvalue:COALESCE(toFloat(prop.qvalue),null),
                                sampleType:COALESCE(prop.sampleType,null)}) 
                        CREATE (e) -[:FROM]-> (s) 
                        CREATE (e)-[:HIT]->(t) with metabo, e
                        unwind metabo as met
                        CALL apoc.get.nodes(met) yield node 
                        CREATE (node)-[:HIT]->(e)")
        
        runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=file,source=source["sid"]$sid))

      }else{ ##biocrates
        
        query <- paste0("
                        USING PERIODIC COMMIT 500
                        LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS prop 
                        MATCH (:metName {sid:apoc.text.replace(toLower(prop.Metabolite), '[^A-Za-z0-9]', '')})<-[:HAS_NAME]-(m:",source["metaboliteID"]$metaboliteID,")
                        with collect(m) as metabo, prop
                        MATCH (s:source {sid:$properties.source}) 
                        MERGE (t:trait {sid:prop.trait})
                        CREATE (e:mWAS {pvalue:toFloat(prop.pvalue),
                                        trait:COALESCE(prop.trait,null),
                                        beta:COALESCE(prop.beta,null),
                                        SE:COALESCE(prop.SE,null),
                                        qvalue:COALESCE(toFloat(prop.qvalue),null),
                                        sampleType:COALESCE(prop.sampleType,null)}) 
                        CREATE (e) -[:FROM]-> (s)  
                        CREATE (e)-[:HIT]->(t) with metabo, e
                        unwind metabo as met
                        CALL apoc.get.nodes(met) yield node 
                        CREATE (node)-[:HIT]->(e)")
        runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=file,source=source["sid"]$sid))
      }
     
      info(log, capture.output(toc()))
     
    }
  }
}