## metabolite quatitative trait loci (mQTL)
addmQTL <- function() {
  for (file in getFiles("QTLs/mQTL/")) {
    #read info on dataset
    source <- getInfo(file, c("sid", "file_name"), err)
    #skip file if no/incomplete info
    if (!all(is.na(source))) {
      tic(paste0("Added mQTL file ", source["file_name"]))
      #check/add source node for dataset
      checkSource(source, log)
      
      if(source["metaboliteID"]$metaboliteID =="COMP_ID"){
        
        query <- paste("
                       USING PERIODIC COMMIT 500
                       LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props with props
                       where exists(props.metabolite)
                       MERGE (n:metName {sid:apoc.text.replace(toLower(props.metabolite), '[^A-Za-z0-9]', '')})
                       MERGE (m:metMetabolon:metMeasured {COMP_ID:props.COMP_ID})
                       ON CREATE SET m.sid = props.COMP_ID + '_unknown'
                       MERGE (n)<-[:HAS_NAME]-(m)  
                       ",sep="")
        runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=file,source=source["sid"]$sid))

        query <- paste("
                       USING PERIODIC COMMIT 500
                       LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props 
                       MATCH (m:metMetabolon {COMP_ID:props.COMP_ID})  
                       with collect(m) as metabo, props
                       MATCH (s:source {sid:$properties.source})
                       MATCH (snp:SNP {sid:props.rsID})
                       CREATE (n:mQTL {
                                       pvalue:toFloat(props.pvalue), 
                                       rsID:COALESCE(props.rsID,null), 
                                       fluid:COALESCE(props.fluid,null), 
                                       metabolite:COALESCE(props.metabolite,null),
                                       effect_allele:COALESCE(props.effect_allele,null),
                                       effect_direction:COALESCE(props.effect_direction,null)
                        }) 
                       CREATE (n)-[:FROM]->(s)
                       CREATE (snp)-[:HIT]->(n) with metabo, n
                       unwind metabo as met
                       CALL apoc.get.nodes(met) yield node 
                       CREATE (n)-[:HIT]->(node)
                       ",sep="")
        runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=file,source=source["sid"]$sid))
        
        
      }else{ ##BIOCHEMICAL biocrates
        
        query <- paste("
                       USING PERIODIC COMMIT 500
                       LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props 
                       MATCH (:metName {sid:apoc.text.replace(toLower(props.metabolite), '[^A-Za-z0-9]', '')})<-[:HAS_NAME]-(m:",source["metaboliteID"]$metaboliteID,")
                       with collect(m) as metabo, props
                       MATCH (s:source {sid:$properties.source})
                       MATCH (snp:SNP {sid:props.rsID})
                       CREATE (n:mQTL {
                         pvalue:toFloat(props.pvalue), 
                         rsID:COALESCE(props.rsID,null), 
                         fluid:COALESCE(props.fluid,null), 
                         metabolite:COALESCE(props.metabolite,null),
                         effect_allele:COALESCE(props.effect_allele,null),
                         effect_direction:COALESCE(props.effect_direction,null)
                       })  
                       CREATE (n)-[:FROM]->(s)
                       CREATE (snp)-[:HIT]->(n) with metabo, n
                       unwind metabo as met
                       CALL apoc.get.nodes(met) yield node 
                       CREATE (n)-[:HIT]->(node)
                       ",sep="")
        runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=file,source=source["sid"]$sid))
      }
      
      info(log, capture.output(toc()))
    
   }
  }
}