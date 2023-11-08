## metabolite quatitative trait loci (mQTL)
addTraitQTL <- function() {
  
  for (file in getFiles("QTLs/traitQTL")) {
    #read info on dataset
    source <- getInfo(file, c("sid", "file_name"), err)
    #skip file if no/incomplete info
    if (!all(is.na(source))) {
      tic(paste0("Added traitQTL file ", source["file_name"]))
      #check/add source node for dataset
      checkSource(source, log)
      
      query <- paste("CALL apoc.periodic.iterate(
                     \"LOAD CSV WITH HEADERS FROM 'file:///'+$properties AS props return distinct props.trait as trait\",
                     \"MERGE (:trait {sid:trait})\", {batchSize:500, parallel:true, iterateList:true, concurrency:8, params: {properties:$properties}})
                     ",sep="")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = file)
      
      query <- paste("
                     USING PERIODIC COMMIT 200
                     LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props 
                     MATCH (s:source {sid:$properties.source})
                     MATCH (snp:SNP {sid:props.SNP})
                     MATCH (m:trait {sid:props.trait})  
                     CREATE (n:traitQTL {
                                          pvalue:toFloat(props.pvalue),
                                          trait:COALESCE(props.trait,null), 
                                          bp:COALESCE(props.bp,null),
                                          chromosome:COALESCE(props.chromosome,null), 
                                          effect_allele:COALESCE(props.effect_allele,null),
                                          non_effect_allele:COALESCE(props.non_effect_allele,null),
                                          major_allele:COALESCE(props.major_allele,null),
                                          minor_allele:COALESCE(props.minor_allele,null),
                                          MAF:COALESCE(props.MAF,null),
                                          OR:COALESCE(props.OR,null),
                                          effect_direction:COALESCE(props.effect_direction,null),
                                          bp_38:COALESCE(props.bp_38,null),
                                          EAF:COALESCE(props.EAF,null),
                                          SE:COALESCE(props.SE,null),
                                          sample_size:COALESCE(props.sample_size,null), 
                                          effect_sign:COALESCE(props.effect_sign,null)}) 
                     CREATE (n)-[:FROM]->(s)
                     CREATE (snp)-[:HIT]->(n)-[:HIT]->(m)
                     ",sep="")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=file,source=source["sid"]$sid))
      
      info(log, capture.output(toc()))
    }
  }
}